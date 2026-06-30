import hashlib
import os
import glob
from datetime import date

import streamlit as st  # type: ignore
import requests         # type: ignore

import appV3
import constants


# Helpers

def compute_md5(file_bytes: bytes) -> str:
    return hashlib.md5(file_bytes).hexdigest()


def parse_md5_file(md5_content: str) -> str | None:
    for line in md5_content.strip().splitlines():
        line = line.strip()
        if line:
            return line.split()[0].lower()
    return None


def fetch_pgs_file_url(pgs_id: str, genome_build: str) -> str | None:
    """
    Return the harmonized file URL for a PGS score from the PGS Catalog REST API.

    Args:
        pgs_id:       e.g. "PGS000014"
        genome_build: "GRCh38" or "GRCh37"

    Returns:
        Direct .txt.gz URL string, or None if unavailable.
    """
    url = f"https://www.pgscatalog.org/rest/score/{pgs_id}"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code != 200:
            return None
        data = resp.json()
        harmonized = data.get("ftp_harmonized_scoring_files", {})
        build_key = "GRCh38" if "38" in genome_build else "GRCh37"
        entry = harmonized.get(build_key, {})
        return entry.get("positions", {}).get("url", None)
    except Exception:
        return None


def clear_ld_cache() -> int:
    """Delete all files in the LD cache directory. Returns count deleted."""
    pattern = os.path.join(constants.LD_CACHE_DIR, "*.txt")
    files = glob.glob(pattern)
    for f in files:
        try:
            os.remove(f)
        except OSError:
            pass
    return len(files)


def render_metadata_card(metadata: dict) -> None:
    """Render a compact info card for PGS file metadata in the Streamlit UI."""
    pgs_id   = metadata.get("pgs_id",      metadata.get("PGS ID",       "N/A"))
    name     = metadata.get("pgs_name",    metadata.get("trait_mapped",  "N/A"))
    trait    = metadata.get("trait_efo",   metadata.get("trait_efo_id",  "N/A"))
    build    = metadata.get("genome_build",metadata.get("HmPOS_build",   "N/A"))
    citation = metadata.get("citation",    metadata.get("Citation",       "N/A"))
    pgp_id   = metadata.get("pgp_id",      metadata.get("PGP ID",        "N/A"))
    license_ = metadata.get("license",     metadata.get("License",        "N/A"))
    accessed = metadata.get("date_accessed", date.today().isoformat())

    st.markdown(
        f"""
        <div style="background:#f0f4ff;border-left:4px solid #4a7cdc;
                    border-radius:6px;padding:12px 16px;margin-bottom:12px;
                    font-size:0.88rem;line-height:1.7;">
          <b>PGS ID:</b> {pgs_id} &nbsp;|&nbsp; <b>Publication:</b> {pgp_id}<br>
          <b>Name / Trait:</b> {name}<br>
          <b>Trait EFO:</b> {trait} &nbsp;|&nbsp; <b>Build:</b> {build}<br>
          <b>Author / Citation:</b> {citation}<br>
          <b>License:</b> {license_}<br>
          <b>Date Accessed:</b> {accessed}
        </div>
        """,
        unsafe_allow_html=True,
    )


# ─── Page config ─────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="PGS Grep",
    layout="wide",
    page_icon="🧬",
)

# Columns to compress title and description into smaller space
title, right_buffer = st.columns([5,7])

with title:
    st.title("🧬 PGS Grep")
    st.markdown(
        "PGS Grep is a utility for searching SNPs within Polygenic Score (PGS) files from the [PGS Catalog](https://www.pgscatalog.org/). " \
        "Locate variant information across PGS datasets and recover related results using linkage equilibrium (LD) data fetched from [LDlink API](https://ldlink.nih.gov/apiaccess)."
    )
st.divider()

# ─── Sidebar ─────────────────────────────────────────────────────────────────
st.sidebar.header("🛠️ Search Configuration")

token_input = st.sidebar.text_input(
    "LDLink API Token",
    value="",
    type="password",
    help="Get a free token at https://ldlink.nih.gov/apiaccess. Required for first-time queries, identical fetched results are cached for future use.",
)

# ── Target variant ────────────────────────────────────────────────────────────
st.sidebar.subheader("🎯 Target Variant", help = "This tool **searches by chromosomal position**, not by rsID. "
        "You can look up the position of a variant at https://www.ncbi.nlm.nih.gov/gdv.")
target_rsid_input = st.sidebar.text_input(
    "Target rsID (must match search position)",
    value="rs10305420",
    help="This must be within the genomic window you specify below. Used to query the LDLink API.",
)
genome_build = st.sidebar.selectbox("Genome Assembly", ["GRCh38", "GRCh37"], index=0, help="The genome assembly used in the PGS file.")
population   = st.sidebar.selectbox(
    "LD Population (1000G)", ["EUR", "AMR", "AFR", "EAS", "SAS"], index=0, help="The population group for linkage disequilibrium calculation. EUR = European, AMR = Admixed American, AFR = African, EAS = East Asian, SAS = South Asian."
)

st.sidebar.subheader("📍 Genomic Search Window", help = "The scan searches for variants and nearby LD Proxies **in this genomic window**. " \
"To find the position of a variant, you can look it up at https://www.ncbi.nlm.nih.gov/gdv."
)

if "resolved_chr" not in st.session_state:
    st.session_state["resolved_chr"] = 6
if "resolved_pos" not in st.session_state:
    st.session_state["resolved_pos"] = 39_048_860  # GRCh38 default for rs10305420

chromosome  = st.sidebar.number_input(
    "Chromosome #",
    min_value=1, max_value=25,
    value=st.session_state["resolved_chr"],
)
target_pos = st.sidebar.number_input(
    "Center Position (SNP position)",
    value=st.session_state["resolved_pos"],
    help="The scan searches for variants at this base-pair position and LD proxies within the specified window.",
)
window_size = st.sidebar.number_input("Flanking Size (+/-)", value=10_000, step=1_000, help="Creates a window of this many base pairs on either side of the center point.")

start_window = target_pos - window_size
end_window   = target_pos + window_size
st.sidebar.caption(f"Search window: Chr{chromosome}:{start_window:,} - {end_window:,}")

# ── LD Metric selection ───────────────────────────────────────────────────────
st.sidebar.subheader(
    "📐 LD Metric",
    help="Choose which linkage disequilibrium statistic(s) to use when filtering proxy variants. "
         "R² measures the correlation between alleles; D′ measures the maximum possible LD. ")
ld_metric = st.sidebar.radio(
    "Filter proxies by",
    options=["R² only", "D′ only", "R² and D′ (both must pass)"],
    index=0,
    horizontal=False,
    help="Select which LD metric(s) proxies must satisfy to be included.",
)

use_r2     = ld_metric in ("R² only",              "R² and D′ (both must pass)")
use_dprime = ld_metric in ("D′ only",              "R² and D′ (both must pass)")

r2_filter = dprime_filter = None
if use_r2:
    r2_filter = st.sidebar.slider("R² Threshold", 0.0, 1.0, 0.7, 0.05, key="r2_slider")
if use_dprime:
    dprime_filter = st.sidebar.slider("D′ Threshold", 0.0, 1.0, 0.8, 0.05, key="dprime_slider")

proxies_only = st.sidebar.checkbox(
    "Hide Unlinked Variants",
    value=False,
    help="Hide variants in genomic window that do not reach R² and D′ thresholds.",
)

# ─── PGS File Source ──────────────────────────────────────────────────────────
st.subheader("📁 Polygenic Score (PGS) File")

file_to_scan       = None
pgs_source_label   = ""

# ── Tab: Manual upload (backup) ─────────────────────────────────────────────
col_upload, col_md5 = st.columns([3, 2])
with col_upload:
    st.markdown("**Upload PGS File**", help="Download a harmonized PGS file from https://www.pgscatalog.org/ and attach here. The file must be in the `.txt.gz` format. No need to unzip!")
    uploaded_file = st.file_uploader(
        "Harmonized PGS Score File (.txt.gz)", type=["gz"]
    )
with col_md5:
    st.markdown("**MD5 Verification** *(optional)*", help="An MD5 file can be provided to verify the integrity of the uploaded PGS file. If the MD5 hash does not match, the scan will still proceed but results may be unreliable.")
    uploaded_md5 = st.file_uploader(
        "MD5 File (.md5) — leave blank to skip", type=["md5", "txt"], key="md5"
    )
local_file_fallback = "PGS Data/PGS000014_hmPOS_GRCh38.txt.gz"
if uploaded_file is not None:
    file_bytes = uploaded_file.read()
    uploaded_file.seek(0)
    if uploaded_md5 is not None:
        md5_content   = uploaded_md5.read().decode("utf-8", errors="replace")
        expected_hash = parse_md5_file(md5_content)
        actual_hash   = compute_md5(file_bytes)
        if expected_hash is None:
            st.warning("⚠️ Could not parse the MD5 file — skipping integrity check.")
        elif actual_hash.lower() == expected_hash.lower():
            st.success(f"✅ MD5 verified: `{actual_hash}`")
        else:
            st.error(
                f"❌ MD5 mismatch!\n"
                f"- Expected : `{expected_hash}`\n"
                f"- Got      : `{actual_hash}`\n\n"
                "The file may be corrupted or incomplete. Proceed with caution."
            )
    # Uploaded file takes precedence over API-fetched file
    file_to_scan     = uploaded_file
    pgs_source_label = uploaded_file.name
elif file_to_scan is None and os.path.exists(local_file_fallback):
    file_to_scan     = local_file_fallback
    pgs_source_label = local_file_fallback
    st.info(f"💡 Using local dev file: `{local_file_fallback}`")
elif file_to_scan is None:
    st.info(" Upload a file to begin! ")


# ─── Metadata preview (shown as soon as a file is available) ─────────────────
if file_to_scan is not None:
    try:
        preview_meta = appV3.extract_pgs_metadata(file_to_scan)
        # Reset seek position after metadata read
        if hasattr(file_to_scan, "seek"):
            file_to_scan.seek(0)
        render_metadata_card(preview_meta)
    except Exception:
        pass  # Don't crash if metadata extraction fails on a malformed file


# ─── Cache status indicator ───────────────────────────────────────────────────
cache_path = os.path.join(
    constants.LD_CACHE_DIR,
    f"LD{target_rsid_input}_{population}_{genome_build.replace(' ', '_')}.txt",
)
if os.path.exists(cache_path):
    st.info(f"🗄️ LD cache found for **{target_rsid_input}** ({population}, {genome_build})")
else:
    st.caption(f"No LD cache found at `{cache_path}`: will query LDlink API on scan.")


# ─── Scan execution ───────────────────────────────────────────────────────────
if file_to_scan and st.button("🚀 Execute Genomic Scan", type="primary"):

    engine = appV3.PGSScanEngine(token=token_input)

    # Step 1 — LD proxies
    with st.spinner("📡 Fetching / loading LD proxy map…"):
        engine.fetch_ld_proxies(
            target_rsid   = target_rsid_input,
            genome_build  = genome_build,
            r2_threshold  = r2_filter,
            dprime_threshold = dprime_filter,
            population    = population,
        )

    if not engine.ld_map:
        st.warning(
            "⚠️ No LD proxies returned above the threshold(s). "
            "Only the exact target position will be searched."
        )
    else:
        metric_desc = (
            f"R² ≥ {r2_filter}" if use_r2 and not use_dprime
            else f"D′ ≥ {dprime_filter}" if use_dprime and not use_r2
            else f"R² ≥ {r2_filter} and D′ ≥ {dprime_filter}"
        )
        st.success(f"✅ {len(engine.ld_map)} proxy position(s) loaded ({metric_desc})")

    # Step 2 — scan file with progress bar
    # TODO: Loading bar for file scanning, update progress every 500 lines or so.
    st.markdown("**🔍 Scanning PGS file…**")
    progress_bar      = st.progress(0, text="Starting scan…")
    progress_status   = st.empty()          # text slot updated alongside the bar

    def on_scan_progress(lines_read: int, lines_total: int | None, variants_found: int) -> None:
        """
        Callback passed into engine.execute_scan() and called every ~500 lines.

        Args:
            lines_read:    Number of data lines processed so far.
            lines_total:   Total data lines in the file, or None if unknown
                           (e.g. streaming from an uploaded file without a pre-count).
            variants_found: Number of matching variants found so far.
        """
        if lines_total:
            fraction = min(lines_read / lines_total, 1.0)
            pct      = int(fraction * 100)
            progress_bar.progress(fraction, text=f"Scanning… {pct}% ({lines_read:,} / {lines_total:,} lines)")
        else:
            # File size unknown — show a pulsing indeterminate-style counter
            pct_guess = min(lines_read / 500_000, 0.99)   # assume ~500k lines max for pulse cap
            progress_bar.progress(pct_guess, text=f"Scanning… {lines_read:,} lines read")
        progress_status.caption(f"Variants in window so far: **{variants_found}**")

    exact_match, proxy_matches, results_df = engine.execute_scan(
        file_object      = file_to_scan,
        chr_number       = chromosome,
        target_pos       = target_pos,
        start_window     = start_window,
        end_window       = end_window,
        target_rsid      = target_rsid_input,
        progress_callback= on_scan_progress,
    )

    progress_bar.progress(1.0, text="✅ Scan complete!")
    progress_status.empty()

    # ── Metadata card (from scanned file header) ──────────────────────────
    st.subheader("📄 Source File Data")
    render_metadata_card(engine.last_metadata)

    # Step 3 — results dashboard
    st.subheader("📊 Scan Results")

    if exact_match:
        st.balloons()
        st.success(f"✅ **Exact match:** `{target_rsid_input}` found directly at Chr{chromosome}:{target_pos:,}.")
    elif proxy_matches:
        st.info(
            f"🔮 **Proxy match:** Target absent at Chr{chromosome}:{target_pos:,}, but "
            f"**{len(proxy_matches)}** LD-linked proxy variant(s) found in the window."
        )
    else:
        st.error(
            f"❌ **No match:** Neither the target position (Chr{chromosome}:{target_pos:,}) "
            "nor any qualifying proxies were found in this window."
        )

    m1, m2, m3 = st.columns(3)
    m1.metric("Exact Position Match", "Yes" if exact_match else "No")
    m2.metric("LD Proxies Captured",  len(proxy_matches))
    m3.metric("Variants in Window",   len(results_df))

    if not results_df.empty:
        st.subheader("📋 Data Viewer")

        # Apply the proxies_only sidebar filter before rendering the table.
        display_df = results_df.copy()
        if proxies_only:
            display_df = display_df[display_df["Match_Status"].str.contains("PROXY|EXACT", na=False)]
            if display_df.empty:
                st.info("No proxy-matched variants in this window.")
            else:
                st.caption(
                    f"Showing {len(display_df)} proxy-matched variant(s) out of "
                    f"{len(results_df)} total in window."
                )

        def highlight_status(val):
            if "EXACT" in val:
                return "background-color: #d4edda; color: #155724;"
            if "PROXY" in val:
                return "background-color: #cce5ff; color: #004085;"
            if "Unlinked" in val:
                return "background-color: #fff3cd; color: #856404;"
            return ""

        if not display_df.empty:
            styled = display_df.style.map(highlight_status, subset=["Match_Status"])
            st.dataframe(styled, width='stretch', hide_index=True)

        # Build output CSV with metadata header comments
        # Note: export always uses the full results_df regardless of the proxy filter,
        # so the user gets a complete record. The filter is a view-only convenience.
        header_comments = appV3.build_output_header(
            metadata      = engine.last_metadata,
            target_rsid   = target_rsid_input,
            population    = population,
            genome_build  = genome_build,
            r2_threshold  = r2_filter,
            dprime_threshold = dprime_filter,
        )
        csv_data   = results_df.to_csv(index=False)
        csv_bytes  = (header_comments + csv_data).encode("utf-8")

        st.download_button(
            label     = "📥 Export Results as CSV",
            data      = csv_bytes,
            file_name = f"pgs_scan_{target_rsid_input}_{date.today().isoformat()}.csv",
            mime      = "text/csv",
        )
    else:
        st.warning("No variants fell within the specified genomic window.")
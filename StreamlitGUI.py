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
          <b>📄 PGS File Metadata</b><br>
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
st.title("🧬 PGS Grep")
st.markdown(
    "PGS Grep is a utility for searching SNPs within Polygenic Score (PGS) files from the [PGS Catalog](https://www.pgscatalog.org/). " \
    "It helps users locate variant information across PGS datasets and recover " 
    "related results using linkage equilibrium (LD) data from the "
    "[1000 Genome project](https://www.internationalgenome.org/). " \
    "This LD information is fetched from the public [LDlink API](https://ldlink.nih.gov/apiaccess)."
)

# ─── Sidebar ─────────────────────────────────────────────────────────────────
st.sidebar.header("🛠️ Configuration")

token_input = st.sidebar.text_input(
    "LDLink API Token",
    value="",
    type="password",
    help="Required for first-time queries. Cached results do not need a token.",
)


# ── Target variant ────────────────────────────────────────────────────────────
st.sidebar.subheader("🎯 Target Variant", help = "This tool **searches by chromosomal position** (not by rsID). "
        "Set the chromosome and position manually. You can look up the position of a variant at https://www.ncbi.nlm.nih.gov/gdv")
target_rsid_input = st.sidebar.text_input(
    "Target rsID Name (must match search position)",
    value="rs10305420",
    help="Used to caption your scan results and query the LDLink API",
)
genome_build = st.sidebar.selectbox("Genome Assembly", ["GRCh38", "GRCh37"], index=0, help="The genome assembly used in the PGS file.")
population   = st.sidebar.selectbox(
    "LD Population (1000G)", ["EUR", "AMR", "AFR", "EAS", "SAS"], index=0, help="The population group for linkage disequilibrium calculation. EUR = European, AMR = Admixed American, AFR = African, EAS = East Asian, SAS = South Asian."
)
r2_filter = st.sidebar.slider("R² Linkage Threshold", 0.0, 1.0, 0.7, 0.05)

st.sidebar.subheader("📍 Genomic Search Position", help = "The scan searches for variants **at this chromosomal position** (and nearby "
    "LD proxies). Use the button below to auto-fill from the rsID above."
)

# ── Auto-fill button ──────────────────────────────────────────────────────────
if "resolved_chr" not in st.session_state:
    st.session_state["resolved_chr"] = 6
if "resolved_pos" not in st.session_state:
    st.session_state["resolved_pos"] = 39_048_860  # GRCh38 default for rs10305420

chromosome  = st.sidebar.number_input(
    "Chromosome #",
    min_value=1, max_value=23,
    value=st.session_state["resolved_chr"],
)
target_pos = st.sidebar.number_input(
    "Center Position (BP) ← scan target",
    value=st.session_state["resolved_pos"],
    help="The scan looks for variants AT this base-pair position (exact match) and LD proxies within the window below.",
)
window_size = st.sidebar.number_input("Flanking Size (+/- BP)", value=10_000, step=1_000)

start_window = target_pos - window_size
end_window   = target_pos + window_size
st.sidebar.caption(f"Search window: Chr{chromosome}:{start_window:,} – {end_window:,}")

# ─── PGS File Source ──────────────────────────────────────────────────────────
st.subheader("📁 PGS Score File")

#source_tab_api, source_tab_upload = st.tabs(["🌐 Fetch from PGS Catalog", "📤 Upload file manually"])

file_to_scan       = None
pgs_source_label   = ""

# ── Tab: Manual upload (backup) ─────────────────────────────────────────────
col_upload, col_md5 = st.columns([3, 2])
with col_upload:
    uploaded_file = st.file_uploader(
        "**Harmonized PGS Score File (.txt.gz)**", type=["gz"]
    )
with col_md5:
    st.markdown("**MD5 Checksum Verification** *(optional)*")
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
    st.info("📥 Upload a harmonized `.txt.gz` PGS Catalog file, or use the Fetch tab above.")


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
            target_rsid  = target_rsid_input,
            genome_build = genome_build,
            r2_threshold = r2_filter,
            population   = population,
        )

    if not engine.ld_map:
        st.warning(
            "⚠️ No LD proxies returned above the R² threshold. "
            "Only the exact target position will be searched."
        )
    else:
        st.success(f"✅ {len(engine.ld_map)} proxy position(s) loaded (R² ≥ {r2_filter})")

    # Step 2 — scan file
    with st.spinner("🔍 Scanning PGS file…"):
        exact_match, proxy_matches, results_df = engine.execute_scan(
            file_object   = file_to_scan,
            chr_number    = chromosome,
            target_pos    = target_pos,
            start_window  = start_window,
            end_window    = end_window,
            target_rsid   = target_rsid_input,
        )

    # ── Metadata card (from scanned file header) ──────────────────────────
    st.subheader("📄 Source File Metadata")
    render_metadata_card(engine.last_metadata)

    # Step 3 — results dashboard
    st.subheader("📊 Scan Results")

    if exact_match:
        st.balloons()
        st.success(f"🎉 **Exact match:** `{target_rsid_input}` found directly at Chr{chromosome}:{target_pos:,}.")
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
        st.subheader("📋 Genomic Data Viewer")

        def highlight_status(val):
            if "EXACT" in val:
                return "background-color: #d4edda; color: #155724;"
            if "PROXY" in val:
                return "background-color: #cce5ff; color: #004085;"
            if "Unlinked" in val:
                return "background-color: #fff3cd; color: #856404;"
            return ""

        styled = results_df.style.map(highlight_status, subset=["Match_Status"])
        st.dataframe(styled, use_container_width=True, hide_index=True)

        # Build output CSV with metadata header comments
        header_comments = appV3.build_output_header(
            metadata      = engine.last_metadata,
            target_rsid   = target_rsid_input,
            population    = population,
            genome_build  = genome_build,
            r2_threshold  = r2_filter,
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
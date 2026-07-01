import hashlib
import io
import os
from datetime import date

import streamlit as st  # type: ignore

import appV3
import constants


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def compute_md5(file_bytes: bytes) -> str:
    """Return the MD5 hex digest of raw file bytes."""
    return hashlib.md5(file_bytes).hexdigest()


def parse_md5_file(md5_content: str) -> str | None:
    """Pull the first hash token out of an uploaded .md5 checksum file."""
    for line in md5_content.strip().splitlines():
        line = line.strip()
        if line:
            return line.split()[0].lower()
    return None


def render_metadata_card(metadata: dict) -> None:
    """Render a compact info card for PGS file metadata in the Streamlit UI."""
    pgs_id = metadata.get("pgs_id", "N/A")
    name = metadata.get("pgs_name", metadata.get("trait_mapped", "N/A"))
    trait = metadata.get("trait_efo", metadata.get("trait_efo_id", "N/A"))
    # Prefer the harmonized build (hmpos_build): it reflects the actual
    # coordinates in this file (hm_chr/hm_pos). The legacy 'genome_build'
    # field instead records the original study's assembly, which is often
    # GRCh37 even in a file harmonized to GRCh38, so it's only used as a
    # fallback when hmpos_build is absent (e.g. older V1 files).
    build = metadata.get("hmpos_build", metadata.get("genome_build", "N/A"))
    citation = metadata.get("citation", "N/A")
    pgp_id = metadata.get("pgp_id", "N/A")
    license_ = metadata.get("license", "N/A")
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


def normalize_build_label(raw_build: str | None) -> str:
    """
    Map a free-text genome-build string pulled from PGS file metadata
    (e.g. 'GRCh38', 'hg38', 'GRCh37') onto one of the two selectbox
    options used by this app: 'GRCh38' or 'GRCh37'.
    """
    if not raw_build:
        return "GRCh38"
    text = str(raw_build)
    if "37" in text:
        return "GRCh37"
    if "38" in text:
        return "GRCh38"
    return "GRCh38"


def get_display_build(metadata: dict | None) -> str:
    """Return the human-readable genome build string extracted from PGS file metadata, or 'Unknown'."""
    if not metadata:
        return "Unknown"
    raw_build = metadata.get("hmpos_build", metadata.get("genome_build"))
    return str(raw_build) if raw_build else "Unknown"


def reset_downstream_state() -> None:
    """
    Clear everything downstream of the source file: any cached scan
    results, match metrics, and the previous metadata preview. Called
    whenever a newly-uploaded file is detected, so stale results from a
    prior file can never be displayed alongside — or exported with — a
    new one.
    """
    st.session_state["scan_results"] = None
    st.session_state["preview_metadata"] = None


# ─────────────────────────────────────────────────────────────────────────────
# Session state
# ─────────────────────────────────────────────────────────────────────────────

STEP_ORDER = ["welcome", "upload", "rsid_guide", "ld_choice", "ld_auth", "config", "execute"]
STEP_DISPLAY_NUMBER = {
    "welcome": "1",
    "upload": "2",
    "rsid_guide": "3",
    "ld_choice": "4",
    "ld_auth": "4.5",
    "config": "5",
    "execute": "6",
}
STEP_TITLE = {
    "welcome": "Welcome",
    "upload": "Upload File",
    "rsid_guide": "Locate Your RSID",
    "ld_choice": "LD Proxies?",
    "ld_auth": "API Token",
    "config": "Search Setup",
    "execute": "Results",
}
TOTAL_STEPS = 6  # counts the conditional 4.5 step as part of step 4


def init_session_state() -> None:
    """Populate st.session_state with default values the first time the app runs."""
    defaults = {
        "wizard_step": "welcome",
        # Step 2 — file acquisition
        "pgs_file_bytes": None,
        "pgs_file_name": None,
        "pgs_file_signature": None,
        "preview_metadata": None,
        # Step 3 / 3.5 — LD proxies
        "want_ld_proxies": "No — scan target position only",
        "ldlink_token": "",
        # Step 4 — search configuration
        "target_rsid": "rs10305420",
        "genome_build": "GRCh38",
        "population": "EUR",
        "chromosome": 6,
        "target_pos": 39_048_860,
        "window_size": 10_000,
        "ld_metric": "R² only",
        "r2_filter": 0.7,
        "dprime_filter": 0.8,
        "proxies_only": True,
        # Step 5 — results cache
        "scan_results": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def wants_ld_proxies() -> bool:
    """True if the user opted into LD proxy searching in Step 3."""
    return st.session_state["want_ld_proxies"].startswith("Yes")


def active_steps() -> list[str]:
    """The ordered list of steps that actually apply, given the user's choices so far."""
    return [s for s in STEP_ORDER if s != "ld_auth" or wants_ld_proxies()]


def go_to_next_step() -> None:
    """Advance the wizard to the next applicable step (skips 4.5 when not needed)."""
    steps = active_steps()
    idx = steps.index(st.session_state["wizard_step"])
    if idx < len(steps) - 1:
        st.session_state["wizard_step"] = steps[idx + 1]


def go_to_prev_step() -> None:
    """Move the wizard back to the previous applicable step (skips 4.5 when not needed)."""
    steps = active_steps()
    idx = steps.index(st.session_state["wizard_step"])
    if idx > 0:
        st.session_state["wizard_step"] = steps[idx - 1]


def render_progress_indicator() -> None:
    """Draw a slim progress bar + 'Step X of Y' caption at the top of the wizard."""
    steps = active_steps()
    current = st.session_state["wizard_step"]
    idx = steps.index(current) + 1
    st.progress(idx / len(steps))
    st.caption(f"Step {STEP_DISPLAY_NUMBER[current]} of {TOTAL_STEPS} &nbsp;·&nbsp; **{STEP_TITLE[current]}**")


# ─────────────────────────────────────────────────────────────────────────────
# Page config
# ─────────────────────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="PGS Grep",
    layout="centered",
    page_icon="🧬",
)

# Custom typography: local Roboto .ttf files served from a /fonts directory
# alongside this script, plus minimal spacing / click-target CSS, plus a
# small fixed-position credit caption in the corner of the app.
st.markdown(
    """
    <style>
    @font-face {
        font-family: 'Roboto';
        src: url('/fonts/Roboto-Regular.ttf') format('truetype');
        font-weight: 400;
        font-style: normal;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('/fonts/Roboto-Medium.ttf') format('truetype');
        font-weight: 500;
        font-style: normal;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('/fonts/Roboto-Bold.ttf') format('truetype');
        font-weight: 700;
        font-style: normal;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('/fonts/Roboto-Italic.ttf') format('truetype');
        font-weight: 400;
        font-style: italic;
    }
    html, body, [class*="css"] {
        font-family: 'Roboto', sans-serif !important;
    }
    div.stButton > button {
        padding: 0.6em 1.4em;
        font-size: 1.05rem;
        border-radius: 8px;
    }
    div[data-testid="stRadio"] label {
        font-size: 1.02rem;
        padding: 4px 0;
    }
    #app-credit-caption {
        position: fixed;
        bottom: 6px;
        right: 10px;
        font-size: 0.7rem;
        color: #9a9a9a;
        z-index: 9999;
        opacity: 0.8;
        pointer-events: none;
    }
    #app-credit-caption a {
        color: #9a9a9a;
        pointer-events: auto;
    }
    </style>
    <div id="app-credit-caption">
        Charlie Buhanan, PGS Grep V1.1, 2026, Download for better results here:
        <a href="https://github.com/CharlieBuhanan/PGS-Grep-Public" target="_blank">
        https://github.com/CharlieBuhanan/PGS-Grep-Public</a>
    </div>
    """,
    unsafe_allow_html=True,
)

init_session_state()

st.title("🧬 PGS Grep")
st.caption(
    "A step-by-step application for locating a single SNP RSID (and its LD proxies) "
    "inside a Polygenic Score (PGS) file."
)
render_progress_indicator()
st.divider()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1: Welcome & Requirements (Start Page)
# ─────────────────────────────────────────────────────────────────────────────

def render_step1_welcome() -> None:
    """
    Introduce the app and list what the user needs before starting.
    Ends with a single 'Get Started' button that begins the wizard.
    """
    st.header("Welcome")
    st.markdown(
        "**PGS Grep** scans a Polygenic Score (PGS) file from the "
        "[PGS Catalog](https://www.pgscatalog.org/) to find a single target SNP "
        "(by chromosomal position), and can optionally recover related variants "
        "using linkage disequilibrium (LD) data from the "
        "[LDlink API](https://ldlink.nih.gov/apiaccess)."
    )

    st.subheader("Before you start, consider gathering:")
    st.markdown(
        """
- A **target SNP RSID** (ex: rs10305420) that you want to locate
- 📁 A Polygenic Score (PGS) file (or its PGS Catalog / publication ID)
- Optional: [LDLink API Token](https://ldlink.nih.gov/apiaccess) *(only if you want to consider Linkage Disequilibrium in your results)*
        """
    )

    st.info("The next steps will guide you through finding a suitable PGS file, configuring the search, and finding an LDLink Token if necessary.")

    st.write("")
    if st.button("Get Started", type="primary", width='stretch'):
        go_to_next_step()
        st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2: File Acquisition & Upload
# ─────────────────────────────────────────────────────────────────────────────

def render_step2_upload() -> None:
    """
    Guide the user to download a harmonized PGS file, let them upload it,
    validate it, and cache the raw bytes plus a metadata preview in
    session_state. Detects newly-uploaded files and clears any stale
    downstream results so they can never be mixed with a new source file.
    """
    st.header("📁 Get Your PGS File")

    with st.expander("How do I download a harmonized PGS file?", expanded=False):
        st.markdown(
            """
1. Go to the [PGS Catalog](https://www.pgscatalog.org/) and find your score (e.g. `PGS000014`).
2. Open its **Download Score** section.
3. Navigate to the **Harmonized** tab folder.
4. Download the **harmonized** scoring file. Filename ends in `_hmPOS_GRCh38.txt.gz`
   or `_hmPOS_GRCh37.txt.gz`, depending on build.
5. Download the optional corresponding MD5 file if you want to verify the download's integrity.
6. No need to unzip. Upload the `.txt.gz` file directly below and the MD5 file if you have it.
            """
        )

    col_upload, col_md5 = st.columns([3, 2])
    with col_upload:
        st.markdown("**Harmonized PGS File** *(`.txt.gz`)*")
        uploaded_file = st.file_uploader(
            "Upload file", type=["gz"], label_visibility="collapsed"
        )
    with col_md5:
        st.markdown("**MD5 Checksum** *(optional)*")
        uploaded_md5 = st.file_uploader(
            "Upload MD5", type=["md5", "txt"], key="md5_uploader",
            label_visibility="collapsed",
        )

    if uploaded_file is not None:
        file_bytes = uploaded_file.read()
        file_name = uploaded_file.name

        # ── File Upload Reset ────────────────────────────────────────────
        # Compare a signature (name + content hash) of this upload against
        # whatever was cached before. A mismatch means a genuinely new
        # file has arrived, so every downstream result / cached value
        # tied to the old file must be cleared out immediately.
        signature = f"{file_name}:{compute_md5(file_bytes)}"
        if signature != st.session_state["pgs_file_signature"]:
            reset_downstream_state()
            st.session_state["pgs_file_signature"] = signature
            st.session_state["pgs_file_bytes"] = file_bytes
            st.session_state["pgs_file_name"] = file_name

        if uploaded_md5 is not None:
            md5_content = uploaded_md5.read().decode("utf-8", errors="replace")
            expected_hash = parse_md5_file(md5_content)
            actual_hash = compute_md5(file_bytes)
            if expected_hash is None:
                st.warning("Could not parse the MD5 file — skipping integrity check.")
            elif actual_hash.lower() == expected_hash.lower():
                st.success(f"MD5 verified: `{actual_hash}`")
            else:
                st.error(
                    f"MD5 mismatch! Expected `{expected_hash}`, got `{actual_hash}`. "
                    "The file may be corrupted — proceed with caution."
                )

    # Metadata preview + validation, using whatever is currently cached.
    if st.session_state["pgs_file_bytes"] is not None:
        st.success(f"📄 File ready: `{st.session_state['pgs_file_name']}`")
        if st.session_state["preview_metadata"] is None:
            try:
                preview_file = io.BytesIO(st.session_state["pgs_file_bytes"])
                preview_meta = appV3.extract_pgs_metadata(preview_file)
                st.session_state["preview_metadata"] = preview_meta
                # Auto-map genome build from file metadata for later steps.
                raw_build = preview_meta.get("hmpos_build", preview_meta.get("genome_build"))
                if raw_build:
                    st.session_state["genome_build"] = normalize_build_label(raw_build)
            except Exception:
                st.session_state["preview_metadata"] = None
                st.warning("Could not read metadata from this file, but you can still proceed.")
        if st.session_state["preview_metadata"]:
            render_metadata_card(st.session_state["preview_metadata"])
    else:
        st.info("Upload a file to continue.")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = st.session_state["pgs_file_bytes"] is not None
        if st.button("Next", type="primary", width='stretch', disabled=not can_advance):
            go_to_next_step()
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 3: RSID Location Guide
# ─────────────────────────────────────────────────────────────────────────────

def render_step3_rsid_guide() -> None:
    """
    Help the user find the physical genomic coordinates (chromosome +
    base-pair position) of their target RSID, using the genome build
    parsed from the PGS file they just uploaded as the correct reference
    assembly to select in NCBI's Genome Data Viewer.
    """
    st.header("Locate Your Target RSID")
    st.markdown(
        "This tool searches by **chromosomal position**, not by rsID directly, "
        "so you'll need to look up your target SNP's coordinates first."
    )

    detected_build = get_display_build(st.session_state["preview_metadata"])

    st.subheader("How to find the coordinates")
    st.markdown(
        f"""
1. Open the [NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/gdv).
2. Search for your rsID (e.g. `rs10305420`).
3. **Set the reference assembly to match your uploaded PGS file: `{detected_build}`.**
   Coordinates differ between assemblies, so using the wrong one will point you
   to the wrong position.
4. Note down the **chromosome number** and **base-pair position** shown for
   that assembly. You'll need this when setting the search window.
        """
    )

    if detected_build == "Unknown":
        st.warning(
            "Could not detect a genome build from your file's metadata. "
            "Double-check the assembly manually before proceeding."
        )
    else:
        st.info(f"Detected genome build from your uploaded file: **{detected_build}**")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        if st.button("Next", type="primary", width='stretch'):
            go_to_next_step()
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4: LD Proxy Selection
# ─────────────────────────────────────────────────────────────────────────────

def render_step4_ld_choice() -> None:
    """
    Ask whether the user wants to search for LD proxies of the target SNP,
    with a short plain-language explanation to inform the choice.
    """
    st.header("Search for LD Proxies?")
    st.markdown(
        """
**LD proxies** are nearby SNPs that tend to be inherited together with your
target SNP (they're in *linkage disequilibrium*). If your exact target SNP
isn't in the PGS file, a proxy can act as a stand-in that captures a
similar genetic signal.
        """
    )

    st.session_state["want_ld_proxies"] = st.radio(
        "Would you like to search for LD Proxies?",
        ["No, scan target position only", "Yes, also search LD proxies (requires a free token)"],
        index=0 if st.session_state["want_ld_proxies"].startswith("No") else 1,
    )

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        if st.button("Next", type="primary", width='stretch'):
            go_to_next_step()
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 4.5: LDLink API Authentication (Conditional)
# ─────────────────────────────────────────────────────────────────────────────

def render_step4_5_ld_auth() -> None:
    """
    Shown only when the user opted into LD proxy search in Step 4.
    Guides them to a free LDLink token and collects it securely.
    """
    st.header("LDLink API Token")
    st.markdown(
        "LD proxy lookups are powered by the **LDlink API**, which requires a "
        "free personal access token."
    )
    st.markdown("[Get a free token here (ldlink.nih.gov)](https://ldlink.nih.gov/apiaccess)")

    st.session_state["ldlink_token"] = st.text_input(
        "LDLink API Token",
        value=st.session_state["ldlink_token"],
        type="password",
        help="Your token is only used for this session's LDlink API calls. "
             "Identical queries are cached locally for future runs.",
    )

    if not st.session_state["ldlink_token"].strip():
        st.warning("A token is required to fetch LD proxies. You can leave this blank and go Back to skip LD proxies instead.")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = bool(st.session_state["ldlink_token"].strip())
        if st.button("Next", type="primary", width='stretch', disabled=not can_advance):
            go_to_next_step()
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 5: Genomic Window & Search Configurations
# ─────────────────────────────────────────────────────────────────────────────

def render_step5_config() -> None:
    """
    Collect the target SNP, genomic search window, LD thresholds (if
    applicable), and output filtering preferences. Genome build is
    pre-filled from Step 2's metadata extraction when available.
    """
    st.header("Search Configuration")

    st.subheader("Target Variant")
    st.caption(
        "Enter the chromosome and position you looked up in the previous step."
    )
    st.session_state["target_rsid"] = st.text_input(
        "Target rsID (must match center position)",
        value=st.session_state["target_rsid"],
    )

    col_build, col_chr = st.columns(2)
    with col_build:
        build_options = ["GRCh38", "GRCh37"]
        st.session_state["genome_build"] = st.selectbox(
            "Genome Assembly",
            build_options,
            index=build_options.index(st.session_state["genome_build"]),
            help="Auto-filled from your uploaded PGS file when available.",
        )
    with col_chr:
        st.session_state["chromosome"] = st.number_input(
            "Chromosome #", min_value=1, max_value=25,
            value=st.session_state["chromosome"],
        )

    st.session_state["target_pos"] = st.number_input(
        "Center Position (target variant position)",
        value=st.session_state["target_pos"],
        help="The scan searches for variants at this base-pair position, plus the window below.",
    )

    st.subheader("Genomic Search Window")
    st.session_state["window_size"] = st.number_input(
        "Flanking Size (+/- base pairs)",
        value=st.session_state["window_size"], step=1_000,
        help="Creates a window this many base pairs on either side of the center position.",
    )
    start_window = st.session_state["target_pos"] - st.session_state["window_size"]
    end_window = st.session_state["target_pos"] + st.session_state["window_size"]
    st.caption(
        f"Search window: **Chr{st.session_state['chromosome']}:"
        f"{start_window:,} – {end_window:,}**"
    )

    if wants_ld_proxies():
        st.subheader("LD Proxy Settings")
        st.session_state["population"] = st.selectbox(
            "LD Population (1000 Genomes)",
            ["EUR", "AMR", "AFR", "EAS", "SAS"],
            index=["EUR", "AMR", "AFR", "EAS", "SAS"].index(st.session_state["population"]),
            help="EUR = European, AMR = Admixed American, AFR = African, EAS = East Asian, SAS = South Asian.",
        )
        st.session_state["ld_metric"] = st.radio(
            "Filter proxies by",
            ["R² only", "D′ only", "R² and D′ (both must pass)"],
            index=["R² only", "D′ only", "R² and D′ (both must pass)"].index(st.session_state["ld_metric"]),
            help="R² measures correlation between alleles; D′ measures maximum possible LD.",
        )
        use_r2 = st.session_state["ld_metric"] in ("R² only", "R² and D′ (both must pass)")
        use_dprime = st.session_state["ld_metric"] in ("D′ only", "R² and D′ (both must pass)")
        if use_r2:
            st.session_state["r2_filter"] = st.slider(
                "R² Threshold", 0.0, 1.0, st.session_state["r2_filter"], 0.05
            )
        if use_dprime:
            st.session_state["dprime_filter"] = st.slider(
                "D′ Threshold", 0.0, 1.0, st.session_state["dprime_filter"], 0.05
            )

    st.subheader("📋 Output Filtering")
    st.session_state["proxies_only"] = st.checkbox(
        "Hide unlinked variants",
        value=st.session_state["proxies_only"],
        help="Hide variants in the window that are neither the exact target nor a qualifying LD proxy. "
             "Applies to both the on-screen table and the exported CSV.",
    )

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        if st.button("Next", type="primary", width='stretch'):
            go_to_next_step()
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# STEP 6: Execution, Caching & Download
# ─────────────────────────────────────────────────────────────────────────────

def render_scan_summary() -> None:
    """Print a plain-language summary of every configuration choice before running the scan."""
    s = st.session_state
    start_window = s["target_pos"] - s["window_size"]
    end_window = s["target_pos"] + s["window_size"]

    st.subheader("Summary")
    lines = [
        f"- **File:** `{s['pgs_file_name']}`",
        f"- **Target rsID:** `{s['target_rsid']}` at Chr{s['chromosome']}:{s['target_pos']:,}",
        f"- **Search window:** Chr{s['chromosome']}:{start_window:,} – {end_window:,}",
        f"- **Genome build:** {s['genome_build']}",
        f"- **LD Proxies:** {'Enabled (' + s['population'] + ', ' + s['ld_metric'] + ')' if wants_ld_proxies() else 'Disabled'}",
        f"- **Hide unlinked variants:** {'Yes' if s['proxies_only'] else 'No'}",
    ]
    st.markdown("\n".join(lines))


def render_cache_status() -> None:
    """Show whether an LD cache file already exists for the current query, to set expectations."""
    s = st.session_state
    if not wants_ld_proxies():
        return
    cache_path = os.path.join(
        constants.LD_CACHE_DIR,
        f"LD{s['target_rsid']}_{s['population']}_{s['genome_build'].replace(' ', '_')}.txt",
    )
    if os.path.exists(cache_path):
        st.caption(f"LD cache found for **{s['target_rsid']}** — repeated LDlink API calls can be skipped.")
    else:
        st.caption("No LD cache found yet: will query the LDlink API on scan.")


def run_scan() -> None:
    """
    Execute the genomic scan using the cached PGS file and configuration,
    then store the results in session_state so they survive accidental
    reruns. The progress bar is deliberately capped below 100% for the
    duration of the scan and is only set to 1.0 once execute_scan() has
    actually returned, so it can never visually claim completion before
    the underlying work is done.
    """
    s = st.session_state
    start_window = s["target_pos"] - s["window_size"]
    end_window = s["target_pos"] + s["window_size"]

    engine = appV3.PGSScanEngine(token=s["ldlink_token"])
    status_placeholder = st.empty()

    with status_placeholder.container():
        if wants_ld_proxies():
            with st.spinner("Fetching / loading LD proxy map…"):
                engine.fetch_ld_proxies(
                    target_rsid=s["target_rsid"],
                    genome_build=s["genome_build"],
                    r2_threshold=s.get("r2_filter"),
                    dprime_threshold=s.get("dprime_filter"),
                    population=s["population"],
                )
            if not engine.ld_map:
                st.warning("No LD proxies returned above the threshold(s). Only the exact target position will be searched.")
            else:
                st.success(f"{len(engine.ld_map)} proxy position(s) loaded")

        st.markdown("**Scanning PGS file…**")
        progress_bar = st.progress(0, text="Starting scan…")
        progress_status = st.empty()

        # Cap displayed progress just below 100% during the scan itself —
        # execute_scan() may report a percent_complete that reaches (or
        # briefly exceeds, near the edge of the window) 1.0 well before
        # the underlying read/parse loop has actually finished, which
        # previously made the bar show "done" prematurely. The bar is
        # only ever set to a true 1.0 after execute_scan() returns below.
        PROGRESS_CAP = 0.99

        def on_scan_progress(current_pos: int, percent_complete: float, variants_found: int) -> None:
            """
            Callback passed into engine.execute_scan(), invoked periodically during the scan.

            Args:
                current_pos:       Base-pair position currently being read.
                percent_complete:  Fraction (0.0-1.0) of the way through the search window.
                variants_found:    Number of matching variants found so far.
            """
            clamped = max(0.0, min(PROGRESS_CAP, percent_complete))
            pct = int(clamped * 100)
            progress_bar.progress(
                clamped,
                text=f"Scanning… {pct}% (position {current_pos:,} of Chr{s['chromosome']}:{start_window:,}-{end_window:,})",
            )
            progress_status.caption(f"Variants in window so far: **{variants_found}**")

        file_object = io.BytesIO(s["pgs_file_bytes"])
        exact_match, proxy_matches, results_df = engine.execute_scan(
            file_object=file_object,
            chr_number=s["chromosome"],
            target_pos=s["target_pos"],
            start_window=start_window,
            end_window=end_window,
            target_rsid=s["target_rsid"],
            progress_callback=on_scan_progress,
        )
        # Only now, after execute_scan() has truly returned, do we mark the bar complete.
        progress_bar.progress(1.0, text="Scan complete!")

    status_placeholder.empty()

    # Cache results so a stray widget interaction doesn't lose them.
    s["scan_results"] = {
        "exact_match": exact_match,
        "proxy_matches": proxy_matches,
        "results_df": results_df,
        "last_metadata": engine.last_metadata,
    }


def render_results() -> None:
    """Render the cached scan results: match banner, metrics, data table, and CSV download."""
    s = st.session_state
    results = s["scan_results"]
    exact_match = results["exact_match"]
    proxy_matches = results["proxy_matches"]
    results_df = results["results_df"]
    last_metadata = results["last_metadata"]

    st.subheader("📄 Source File Data")
    render_metadata_card(last_metadata)

    st.subheader("Scan Results")
    if exact_match:
        st.balloons()
        st.markdown(
            f"""
            <div style="background:#d4edda;border:2px solid #28a745;border-radius:8px;
                        padding:16px 20px;margin-bottom:12px;font-size:1.15rem;
                        font-weight:700;color:#155724;text-align:center;">
              EXACT MATCH FOUND — {s['target_rsid']} at Chr{s['chromosome']}:{s['target_pos']:,}
            </div>
            """,
            unsafe_allow_html=True,
        )
    elif proxy_matches:
        st.markdown(
            f"""
            <div style="background:#cce5ff;border:2px solid #004085;border-radius:8px;
                        padding:16px 20px;margin-bottom:12px;font-size:1.1rem;
                        font-weight:700;color:#004085;text-align:center;">
              PROXY MATCH — target absent at Chr{s['chromosome']}:{s['target_pos']:,}, but
              {len(proxy_matches)} LD-linked proxy variant(s) found in the window.
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.error(
            f"No match: neither the target position (Chr{s['chromosome']}:{s['target_pos']:,}) "
            "nor any qualifying proxies were found in this window."
        )

    m1, m2, _ = st.columns([1, 1, 3])
    m1.metric("Exact Position Match", "Yes" if exact_match else "No")
    m2.metric("LD Proxies Captured", len(proxy_matches))

    if results_df.empty:
        st.warning("No variants fell within the specified genomic window.")
        return

    st.subheader("📋 Data Viewer")
    if s["proxies_only"]:
        export_df = results_df[results_df["Match_Status"].str.contains("PROXY|EXACT", na=False)]
    else:
        export_df = results_df
    display_df = export_df.copy()

    if s["proxies_only"]:
        if display_df.empty:
            st.info("No exact-match or proxy-matched variants in this window.")
        else:
            st.caption(
                f"Showing {len(display_df)} matched variant(s) (unlinked variants hidden) "
                f"out of {len(results_df)} total in window."
            )

    def highlight_status(val):
        """Row-color a Match_Status cell based on whether it's an exact, proxy, or unlinked call."""
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

    # ── CSV metadata injection ──────────────────────────────────────────
    # A greyed-out / commented credit line is prepended ahead of the
    # existing metadata header comments produced by appV3.build_output_header.
    # The GitHub download URL shown in the on-screen corner caption is
    # deliberately excluded from this exported line.
    csv_credit_line = "# Charlie Buhanan, PGS Grep V1.1, 2026\n"
    header_comments = appV3.build_output_header(
        metadata=last_metadata,
        target_rsid=s["target_rsid"],
        population=s.get("population"),
        genome_build=s["genome_build"],
        r2_threshold=s.get("r2_filter"),
        dprime_threshold=s.get("dprime_filter"),
    )
    csv_data = export_df.to_csv(index=False)
    csv_bytes = (csv_credit_line + header_comments + csv_data).encode("utf-8")

    st.write("")
    st.download_button(
        label="📥 Export Results as CSV",
        data=csv_bytes,
        file_name=f"pgs_scan_{s['target_rsid']}_{date.today().isoformat()}.csv",
        mime="text/csv",
        type="primary",
        width='stretch',
    )


def render_step6_execute() -> None:
    """
    Show a configuration summary, run (or re-display) the genomic scan,
    and offer the CSV download once results are available.
    """
    st.header("Run the Scan")
    render_scan_summary()
    render_cache_status()

    st.write("")
    if st.button("Execute Genomic Scan", type="primary", width='stretch'):
        run_scan()
        st.rerun()

    if st.session_state["scan_results"] is not None:
        st.divider()
        render_results()

    st.write("")
    col_back, col_restart = st.columns(2)
    with col_back:
        if st.button("Back to Config", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_restart:
        if st.button("Start Over", width='stretch'):
            for key in list(st.session_state.keys()):
                del st.session_state[key]
            st.rerun()


# ─────────────────────────────────────────────────────────────────────────────
# Router
# ─────────────────────────────────────────────────────────────────────────────

STEP_RENDERERS = {
    "welcome": render_step1_welcome,
    "upload": render_step2_upload,
    "rsid_guide": render_step3_rsid_guide,
    "ld_choice": render_step4_ld_choice,
    "ld_auth": render_step4_5_ld_auth,
    "config": render_step5_config,
    "execute": render_step6_execute,
}

STEP_RENDERERS[st.session_state["wizard_step"]]()
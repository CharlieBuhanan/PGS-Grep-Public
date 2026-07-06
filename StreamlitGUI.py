"""
R²
PGS Grep Streamlit wizard UI.

Architecture notes for developers:

- Widget/durable key pattern: any input widget that is only mounted on a
  single wizard step has its own "*_widget" session_state key. Streamlit
  drops a widget's session_state entry the moment that widget stops being
  drawn on a rerun, so binding a widget's key= directly to the value the
  rest of the app reads would wipe that value out as soon as the user
  navigated away from the step. Each widget's value is therefore copied,
  immediately after the widget is drawn, into a separate durable key
  (e.g. "target_rsid", "chromosome", "target_pos") that is a plain
  session_state entry and persists regardless of which step is currently
  rendering. It's the durable key that the rest of the app (scan summary,
  run_scan, CSV export, and the read-only restatement in Step 5) reads.
- Custom typography/CSS below relies on Roboto being served locally from
  ./static/ (requires enableStaticServing=true in .streamlit/config.toml).
"""

import hashlib
import io
import os
from datetime import date

import streamlit as st  # type: ignore

import appV3
import constants


def compute_md5(file_bytes: bytes) -> str:
    """Return the MD5 hex digest of raw file bytes."""
    return hashlib.md5(file_bytes).hexdigest()


def _sync_target_pos_from_text() -> None:
    """
    Callback for the Center Position text field. Strips thousands
    separators and whitespace (positions copied from NCBI's Genome Data
    Viewer are formatted like "39,048,860", which a native numeric input
    silently rejects on paste) and, if what's left is a valid integer,
    writes it into the numeric 'target_pos' session-state key that the
    rest of the app relies on. Invalid input is left in the text box
    as-is and simply doesn't update 'target_pos', so a bad paste can't
    silently corrupt the search window.

    Reads from the widget's own key ('target_pos_text_widget') rather
    than the durable 'target_pos_text' key, and immediately mirrors the
    raw text into 'target_pos_text' too, per the durable/widget key
    split described in the module docstring.
    """
    raw = st.session_state.get("target_pos_text_widget", "")
    st.session_state["target_pos_text"] = raw
    cleaned = raw.replace(",", "").replace(" ", "").strip()
    if cleaned.lstrip("-").isdigit():
        st.session_state["target_pos"] = int(cleaned)


def parse_md5_file(md5_content: str) -> str | None:
    """Pull the first hash token out of an uploaded .md5 checksum file."""
    for line in md5_content.strip().splitlines():
        line = line.strip()
        if line:
            return line.split()[0].lower()
    return None


def render_metadata_card(metadata: dict) -> None:
    """
    Render a compact info card for PGS file metadata in the Streamlit UI.

    Prefers the harmonized build ('hmpos_build'), which reflects the
    actual coordinates in this file (hm_chr/hm_pos). The legacy
    'genome_build' field records the original study's assembly instead,
    which is often GRCh37 even in a file harmonized to GRCh38, so it's
    only used as a fallback when hmpos_build is absent (e.g. older V1
    files).
    """
    pgs_id = metadata.get("pgs_id", "N/A")
    name = metadata.get("pgs_name", metadata.get("trait_mapped", "N/A"))
    trait = metadata.get("trait_efo", metadata.get("trait_efo_id", "N/A"))
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
    "rsid_guide": "Locate & Enter Target SNP",
    "ld_choice": "LD Proxies?",
    "ld_auth": "API Token",
    "config": "Search Setup",
    "execute": "Results",
}
TOTAL_STEPS = 6


def init_session_state() -> None:
    """Populate st.session_state with default values the first time the app runs."""
    defaults = {
        "wizard_step": "welcome",
        "pgs_file_bytes": None,
        "pgs_file_name": None,
        "pgs_file_signature": None,
        "pgs_file_is_harmonized": None,
        "preview_metadata": None,
        "want_ld_proxies": "No, scan target position only",
        "ldlink_token": "",
        "target_rsid": "rs10305420",
        "genome_build": "GRCh38",
        "population": "Choose...",
        "chromosome": 6,
        "target_pos": 39_048_860,
        "target_pos_text": "39,048,860",
        "window_size": 5_000,
        "ld_metric": "R² only",
        "r2_filter": 0.7,
        "dprime_filter": 0.8,
        "proxies_only": True,
        "scan_results": None,
    }
    for key, value in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = value


def wants_ld_proxies() -> bool:
    """True if the user opted into LD proxy searching in Step 3."""
    return st.session_state["want_ld_proxies"].startswith("Yes")


def go_to_next_step() -> None:
    """
    Advance the wizard to the next step in the fixed STEP_ORDER sequence,
    automatically skipping over 'ld_auth' when LD proxies are disabled.

    STEP_ORDER never changes shape (unlike a dynamically-filtered list),
    so st.session_state["wizard_step"] is always a valid member of it and
    STEP_ORDER.index(...) can never raise ValueError.
    """
    idx = STEP_ORDER.index(st.session_state["wizard_step"]) + 1
    while idx < len(STEP_ORDER) and STEP_ORDER[idx] == "ld_auth" and not wants_ld_proxies():
        idx += 1
    if idx < len(STEP_ORDER):
        st.session_state["wizard_step"] = STEP_ORDER[idx]


def go_to_prev_step() -> None:
    """
    Move the wizard back to the previous step in the fixed STEP_ORDER
    sequence, automatically skipping over 'ld_auth' when LD proxies are
    disabled. See go_to_next_step() for why STEP_ORDER (rather than a
    dynamically-filtered list) is used for the index lookup.
    """
    idx = STEP_ORDER.index(st.session_state["wizard_step"]) - 1
    while idx >= 0 and STEP_ORDER[idx] == "ld_auth" and not wants_ld_proxies():
        idx -= 1
    if idx >= 0:
        st.session_state["wizard_step"] = STEP_ORDER[idx]


def render_progress_indicator() -> None:
    """
    Draw a slim progress bar + 'Step X of Y' caption at the top of the
    wizard. Uses the fixed step numbering (STEP_DISPLAY_NUMBER /
    TOTAL_STEPS) rather than deriving a fraction from position-in-sequence,
    since whether 'ld_auth' is visited or skipped depends on the LD-proxies
    Yes/No choice and would otherwise make the fraction jump around while
    sitting still on the same step.
    """
    current = st.session_state["wizard_step"]
    fraction = float(STEP_DISPLAY_NUMBER[current]) / TOTAL_STEPS
    st.progress(min(fraction, 1.0))
    st.caption(f"Step {STEP_DISPLAY_NUMBER[current]} of {TOTAL_STEPS} &nbsp;·&nbsp; **{STEP_TITLE[current]}**")


st.set_page_config(
    page_title="PGS Grep",
    layout="centered",
    page_icon="🧬",
)

st.markdown(
    """
    <style>
    @font-face {
        font-family: 'Roboto';
        src: url('app/static/Roboto-Regular.ttf') format('truetype');
        font-weight: 400;
        font-style: normal;
        font-display: swap;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('app/static/Roboto-Medium.ttf') format('truetype');
        font-weight: 500;
        font-style: normal;
        font-display: swap;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('app/static/Roboto-Bold.ttf') format('truetype');
        font-weight: 700;
        font-style: normal;
        font-display: swap;
    }
    @font-face {
        font-family: 'Roboto';
        src: url('app/static/Roboto-Italic.ttf') format('truetype');
        font-weight: 400;
        font-style: italic;
        font-display: swap;
    }

    html, body, [class*="css"], .stApp, .stMarkdown, .stApp p,
    .stApp span:not([data-testid="stIconMaterial"]),
    .stApp label, .stApp div {
        font-family: 'Roboto', 'Helvetica Neue', Arial, sans-serif !important;
    }
    code, pre, .stCode, [data-testid="stCodeBlock"] {
        font-family: 'Roboto Mono', 'Courier New', monospace !important;
    }
    h1, h2, h3, h4 {
        font-family: 'Roboto', 'Helvetica Neue', Arial, sans-serif !important;
        font-weight: 700 !important;
        letter-spacing: -0.01em;
        color: #1a2b4c;
    }

    .block-container {
        max-width: 840px;
        padding-top: 2.2rem;
        padding-bottom: 3rem;
    }

    [data-testid="stAppViewContainer"] h1 {
        font-weight: 800 !important;
        color: #000000;
    }

    div.stButton > button {
        padding: 0.62em 1.5em;
        font-size: 1.02rem;
        font-weight: 600;
        border-radius: 10px;
        border: 1px solid #d5dbe6;
        transition: box-shadow 0.15s ease-in-out;
    }
    div.stButton > button:hover {
        box-shadow: 0 2px 6px rgba(16, 24, 40, 0.12);
    }
    div.stDownloadButton > button {
        border-radius: 10px;
        font-weight: 600;
        transition: box-shadow 0.15s ease-in-out;
    }
    div.stDownloadButton > button:hover {
        box-shadow: 0 2px 6px rgba(16, 24, 40, 0.12);
    }

    div[data-testid="stFileUploaderDropzone"] {
        min-height: 112px;
    }

    div[data-testid="stTextInput"] input,
    div[data-testid="stNumberInput"] input,
    div[data-testid="stSelectbox"] div[data-baseweb="select"] {
        border-radius: 8px !important;
    }
    div[data-testid="stRadio"] label {
        font-size: 1.02rem;
        padding: 4px 0;
    }
    div[data-testid="stRadio"] > div {
        gap: 0.35rem;
    }

    div[data-testid="stExpander"] {
        border-radius: 12px;
        border: 1px solid #e3e8f0;
        overflow: hidden;
    }
    div[data-testid="stAlert"] {
        border-radius: 10px;
        padding: 0.9rem 1.1rem;
    }

    div[data-testid="stMetric"] {
        background: #f7f9fc;
        border: 1px solid #e3e8f0;
        border-radius: 12px;
        padding: 0.8rem 1rem;
    }

    h2, h3 {
        margin-top: 1.4rem;
    }
    hr {
        margin: 1.6rem 0;
    }

    div[data-testid="stProgress"] > div > div {
        border-radius: 8px;
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
        Charlie Buhanan, PGS Grep V1.1, 2026, Download for better results at:
        <a href="https://github.com/CharlieBuhanan/PGS-Grep-Public" target="_blank">
        https://github.com/CharlieBuhanan/PGS-Grep-Public</a>
    </div>
    """,
    unsafe_allow_html=True,
)

init_session_state()

st.title("🧬 PGS Grep")
st.caption(
    "An application for locating an SNP (and related variants) inside a Polygenic Score (PGS) file."
)
render_progress_indicator()
st.divider()


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

    with st.expander("ℹ️ Information on Key Terms", expanded=False):
        st.markdown("""**Polygenic Score (PGS)**: a numeric score that estimates an individual's genetic predisposition to a trait or disease, based on the combined effect of many genetic variants.
PGS files are published in the PGS Catalog (https://www.pgscatalog.org/) and detail the genetic variants and how much each contributes to a target disease.

**SNP (Single Nucleotide Polymorphism)**: a single base pair location in the DNA
sequence where individuals commonly differ. This app searches for SNPs by their position in the human genome. 

**LD (Linkage Disequilibrium)**: a measure of how often two nearby genetic
variants are inherited together. Variants close together on a chromosome tend to have a higher LD.

**LD Proxies**: SNPs that are in strong linkage disequilibrium with a target SNP, meaning
they're usually inherited alongside it. They are useful because if an SNP is missing from a PGS file, a proxy that
tracks closely with it can act as a stand-in and capture a similar
genetic signal.

**LDLink API Access Token**: this optional token can be obtained for free at [LDLink API](https://ldlink.nih.gov/apiaccess), hosted and developed by the National Cancer Institute (NCI), 
part of the U.S. National Institutes of Health (NIH). This token authenticates a request for LD proxy data from the LDLink database. Your token is never stored or shared. If you choose to provide one, it will only be used for LDLink's API calls. You can see how your token is
used by inspecting the source code at https://github.com/CharlieBuhanan/PGS-Grep-Public.""")

    st.subheader("Before you start, consider:")
    st.markdown(
        """
- A **target SNP RSID** (ex: rs10305420) that you want to locate
- A Polygenic Score (PGS) file (or a publication on the PGS Catalog)
- Whether or not you want to consider Linkage Disequilibrium in your results (this requires a free LDLink Token from the National Cancer Institute, see the [LDLink API](https://ldlink.nih.gov/apiaccess))
        """
    )

    st.info("The next steps will guide you through finding a suitable PGS file, configuring the search, and obtaining an LDLink Token if necessary.")

    st.write("")
    if st.button("Get Started", type="primary", width='stretch'):
        go_to_next_step()
        st.rerun()

    st.warning(
        "This application is hosted on the free Streamlit "
        "Community Cloud for demonstration purposes. Heavy usage may cause "
        "temporary slowdowns or outages. For the best experience, download and run the application locally "
        "from: https://github.com/CharlieBuhanan/PGS-Grep-Public"
    )

    with st.expander("License", expanded=False):
        st.markdown(
            "PGS Grep is released under the **GNU General Public License v3.0 (GPL-3.0).** You are free to reuse and modify this code under the terms of the GPL-3.0 license, but any derivative work must also be released under the same license."
            " See details at https://www.gnu.org/licenses/gpl-3.0.html"
        )


MAX_PGS_FILE_SIZE_BYTES = 200 * 1024 * 1024


def render_step2_upload() -> None:
    """
    Guide the user to download a harmonized PGS file, let them upload it,
    validate it, and cache the raw bytes plus a metadata preview in
    session_state. Detects newly-uploaded files by comparing a signature
    (name + content hash) against whatever was cached before; a mismatch
    clears every downstream result tied to the old file so it can never
    be mixed with a new source file.
    """
    st.header("📁 Get Your PGS File from the PGS Catalog")

    with st.expander("ℹ️ How do I download a harmonized PGS file?", expanded=False):
        st.markdown(
            """
1. Go to the **[PGS Catalog](https://www.pgscatalog.org/)** and find your score (e.g. `PGS000014`). 
You can also search the Catalog by publication (author name, journal, PGP ID, or PubMed ID).\n
2. Open its **Download Score** section.\n
3. Navigate to the **Harmonized** tab folder.
""")    
        st.markdown("""
4. Download the **harmonized** scoring file. Filename ends in `_hmPOS_GRCh38.txt.gz`
   or `_hmPOS_GRCh37.txt.gz`, depending on the genome build.""", help = "The genome build GRCh38 (hg38) and GRCh37 (hg19) specifies the reference human genome assembly used to define the chromosomal coordinates of the SNPs in the file. Be sure to download the harmonized version of the score  file for your desired genome build.")

        st.markdown("5. Download the optional corresponding MD5 file if you want to verify the download's integrity.", help = "Optional. The .md5 / .txt checksum file from the PGS Catalog is "
                 "used to confirm your download wasn't corrupted or truncated.")
        st.markdown("6. No need to unzip. Upload the `.txt.gz` file directly below and the MD5 file if you have it. 200 megabytes is the maximum size for file uploads.", help="PGS Grep can read the .gz zipped file type directly. No zipped files on the PGS Catalogue exceed 200 MB as of June 2026.")

    col_upload, col_md5 = st.columns([1,1])
    with col_upload:
        st.markdown("**Harmonized PGS File** *(`.txt.gz`)*")
        uploaded_file = st.file_uploader(
            "Upload file", type=["gz"], label_visibility="collapsed"
        )
    with col_md5:
        st.markdown("**MD5 Checksum** *(optional)*", help = "Optional. The .md5 / .txt checksum file from the PGS Catalog is used to confirm your download wasn't corrupted or truncated.")
        uploaded_md5 = st.file_uploader(
            "Upload MD5", type=["md5"], key="md5_uploader",
            label_visibility="collapsed",
        )
    st.caption("200 MB maximum file size. Must be in harmonized format.")

    if uploaded_file is not None:
        if uploaded_file.size > MAX_PGS_FILE_SIZE_BYTES:
            st.error(
                f"**'{uploaded_file.name}' is "
                f"{uploaded_file.size / (1024 * 1024):.1f} MB, which exceeds the 200 MB limit** "
                "for this hosted demo. Please use a smaller/compressed file, or download and run "
                "PGS Grep locally (see the note on the Welcome page) to remove this limit."
            )
        else:
            file_bytes = uploaded_file.read()
            file_name = uploaded_file.name

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
                    st.warning("Could not parse the MD5 file, skipping integrity check.")
                elif actual_hash.lower() == expected_hash.lower():
                    st.success(f"MD5 verified: `{actual_hash}`")
                else:
                    st.error(
                        f"MD5 mismatch! Expected `{expected_hash}`, got `{actual_hash}`. "
                        "The file may be corrupted, proceed with caution."
                    )

    if st.session_state["pgs_file_bytes"] is not None:
        st.success(f"📄 File ready: `{st.session_state['pgs_file_name']}`")
        if st.session_state["preview_metadata"] is None:
            try:
                preview_file = io.BytesIO(st.session_state["pgs_file_bytes"])
                preview_meta = appV3.extract_pgs_metadata(preview_file)
                st.session_state["preview_metadata"] = preview_meta
                raw_build = preview_meta.get("hmpos_build", preview_meta.get("genome_build"))
                if raw_build:
                    st.session_state["genome_build"] = normalize_build_label(raw_build)
                file_name_lower = (st.session_state["pgs_file_name"] or "").lower()
                st.session_state["pgs_file_is_harmonized"] = bool(
                    preview_meta.get("hmpos_build")
                ) or "_hmpos_" in file_name_lower
            except Exception:
                st.session_state["preview_metadata"] = None
                st.session_state["pgs_file_is_harmonized"] = None
                st.warning("Could not read metadata from this file, but you can still proceed.")
        if st.session_state["preview_metadata"]:
            render_metadata_card(st.session_state["preview_metadata"])
        if st.session_state.get("pgs_file_is_harmonized") is False:
            st.error(
                "**This doesn't look like a harmonized PGS file.** No harmonized-build "
                "metadata (`hm_pos`) was found, and the filename doesn't include `_hmPOS_`. "
                "Please re-download the **harmonized** version of this score from the PGS Catalog's "
                "'Harmonized' folder (see the instructions above) before continuing."
            )
    else:
        st.info("Upload a file to continue.")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = (
            st.session_state["pgs_file_bytes"] is not None
            and st.session_state.get("pgs_file_is_harmonized") is not False
        )
        if st.button("Next", type="primary", width='stretch', disabled=not can_advance):
            go_to_next_step()
            st.rerun()


def render_step3_rsid_guide() -> None:
    """
    Help the user find the physical genomic coordinates (chromosome +
    base-pair position) of their target RSID, using the genome build
    parsed from the PGS file they just uploaded as the correct reference
    assembly to select in NCBI's Genome Data Viewer, then collect the
    target rsID, chromosome, and center position directly on this step.

    These three fields are bound to durable session_state keys
    ("target_rsid", "chromosome", "target_pos"/"target_pos_text") via the
    widget/durable key split described in the module docstring, so their
    values survive navigating to other steps or triggering a script
    rerun, and are later restated read-only in Step 5.
    """
    st.header("Locate Your Target RSID")
    st.markdown(
        "An **rsID** (e.g. `rs10305420`) is a unique reference ID assigned to a "
        "specific SNP in NCBI's dbSNP database. This tool searches by "
        "**chromosomal position**, not by rsID directly, so you'll need to look "
        "up your target SNP's coordinates first."
    )

    detected_build = get_display_build(st.session_state["preview_metadata"])

    st.subheader("How to find the coordinates")
    st.markdown(
        f"""
1. Open the [NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/gdv). This resource is managed by the National Library of Medicine and provides a visual interface to explore genomic data.
2. Search for your rsID (e.g. `rs10305420`).
3. **Set the reference assembly to match your uploaded PGS file: `{detected_build}`.**
   Coordinates differ between assemblies, so using the wrong one will point you
   to the wrong position.
4. Note down the **chromosome number** and **base-pair position** shown for
   that assembly, or keep the tab open. You'll need this to fill in the fields below.
        """
    )

    if detected_build == "Unknown":
        st.warning(
            "Could not detect a genome build from your file's metadata. "
            "Double-check the assembly manually before proceeding."
        )
    else:
        st.info(f"Detected genome build from your uploaded file: **{detected_build}**")

    st.subheader("Target Variant")
    st.caption("Enter the rsID, chromosome, and position you looked up above. Default values are shown for rs10305420 as an example.")

    st.text_input(
        "Target rsID (must match center position)",
        value=st.session_state["target_rsid"],
        key="target_rsid_widget",
        help="The rsID (e.g. rs10305420) you looked up above. Used to label "
             "results and, if LD proxies are enabled, as the query variant "
             "for the LDlink lookup.",
    )
    st.session_state["target_rsid"] = st.session_state["target_rsid_widget"]

    st.number_input(
        "Chromosome #", min_value=1, max_value=25,
        value=st.session_state["chromosome"],
        key="chromosome_widget",
        help="The chromosome your target variant is on (1-22, 23=X, 24=Y, 25=MT).",
    )
    st.session_state["chromosome"] = st.session_state["chromosome_widget"]

    st.text_input(
        "Center Position (target variant position)",
        value=st.session_state["target_pos_text"],
        key="target_pos_text_widget",
        on_change=_sync_target_pos_from_text,
        help="The base-pair position of your target variant, e.g. from NCBI's "
             "Genome Data Viewer. You can paste it with or without comma "
             "separators (e.g. 39048860 or 39,048,860). The scan searches "
             "for variants at this position, plus a window configured later.",
    )
    st.session_state["target_pos_text"] = st.session_state["target_pos_text_widget"]
    position_valid = st.session_state["target_pos_text"].replace(",", "").replace(" ", "").strip().lstrip("-").isdigit()
    if not position_valid:
        st.caption("⚠️ Enter digits only (commas are fine)")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = bool(st.session_state["target_rsid"].strip()) and position_valid
        if st.button("Next", type="primary", width='stretch', disabled=not can_advance):
            go_to_next_step()
            st.rerun()


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

    ld_proxy_options = [
        "No, scan target position only",
        "Yes, also search LD proxies (requires a free token)",
    ]
    selection = st.radio(
        "Would you like to search for LD Proxies?",
        ld_proxy_options,
        index=ld_proxy_options.index(st.session_state["want_ld_proxies"]),
        key="want_ld_proxies_widget",
    )
    st.session_state["want_ld_proxies"] = selection

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


def render_step4_5_ld_auth() -> None:
    """
    Shown only when the user opted into LD proxy search in Step 4.
    Guides them to a free LDLink token and collects it securely.
    """
    st.header("LDLink API Token")
    st.markdown(
        "LD proxy lookups are powered by the **LDlink API**, a free service run "
        "by the National Cancer Institute. It requires a personal access "
        "token to authenticate requests. This "
        "app only uses it to make LD queries on your behalf during this session."
    )
    st.caption("LDLink API Access Token: this optional token can be obtained for free at [LDLink API](https://ldlink.nih.gov/apiaccess), hosted and developed by the National Cancer Institute (NCI), "
        "part of the U.S. National Institutes of Health (NIH). This token authenticates a request for LD proxy data from the LDLink database. "
        "Your token is never stored or shared. If you choose to provide one, it will only be used for LDLink's API calls. You can see how your token is "
        "used by inspecting the source code at https://github.com/CharlieBuhanan/PGS-Grep-Public.")
    st.markdown("[Get a free token here (ldlink.nih.gov)](https://ldlink.nih.gov/apiaccess)")
    

    st.text_input(
        "LDLink API Token",
        value=st.session_state["ldlink_token"],
        type="password",
        key="ldlink_token_widget",
        help="Your token is only used for this session's LDlink API calls. "
             "Identical queries are cached locally for future runs.",
    )
    st.session_state["ldlink_token"] = st.session_state["ldlink_token_widget"]

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


def render_step5_config() -> None:
    """
    Restate the target rsID, chromosome, position, and detected genome
    build collected in Step 3 as read-only fields, then collect the
    genomic search window, LD thresholds (if applicable), and output
    filtering preferences. The genome build used for the scan is
    pre-filled from Step 2's metadata extraction when available.
    """
    st.header("Search Configuration")

    st.subheader("Target Variant")
    st.caption("Set in the previous step. Go Back to change these values.")

    st.text_input(
        "Target rsID",
        value=st.session_state["target_rsid"],
        disabled=True,
    )

    col_build, col_chr = st.columns(2)
    with col_build:
        st.text_input(
            "Detected Genome Build (from harmonized file)",
            value=get_display_build(st.session_state["preview_metadata"]),
            disabled=True,
        )
    with col_chr:
        st.number_input(
            "Chromosome #",
            value=st.session_state["chromosome"],
            disabled=True,
        )

    st.text_input(
        "Center Position (target variant position)",
        value=f"{st.session_state['target_pos']:,}",
        disabled=True,
    )

    st.subheader("Genome Assembly")
    build_options = ["GRCh38", "GRCh37"]
    st.selectbox(
        "Genome Assembly Used for Search",
        build_options,
        index=build_options.index(st.session_state["genome_build"]),
        key="genome_build_widget",
        help="The reference genome assembly your coordinates are in. "
             "Auto-filled from your uploaded PGS file when available — "
             "change it only if you looked up coordinates in a different build.",
    )
    st.session_state["genome_build"] = st.session_state["genome_build_widget"]

    st.subheader("Genomic Search Window")
    st.number_input(
        "Flanking Size (+/- base pairs)",
        step=1_000,
        value=st.session_state["window_size"],
        key="window_size_widget",
        help="Creates a window this many base pairs on either side of the center position.",
    )
    st.session_state["window_size"] = st.session_state["window_size_widget"]
    start_window = st.session_state["target_pos"] - st.session_state["window_size"]
    end_window = st.session_state["target_pos"] + st.session_state["window_size"]
    st.caption(
        f"Search window: **Chr{st.session_state['chromosome']}:"
        f"{start_window:,} - {end_window:,}**"
    )

    population_selected = True
    if wants_ld_proxies():
        st.subheader("LD Proxy Settings")
        population_options = ["Choose...", "EUR", "AMR", "AFR", "EAS", "SAS"]
        st.selectbox(
            "LD Population (1000 Genomes)",
            population_options,
            index=population_options.index(st.session_state["population"]),
            key="population_widget",
            help="The 1000 Genomes super-population used to calculate LD. This "
                 "matters because LD patterns differ by ancestry — proxies "
                 "strongly linked in one population may not be linked in "
                 "another. EUR = European, AMR = Admixed American, "
                 "AFR = African, EAS = East Asian, SAS = South Asian. "
                 "Choose the population that best matches your study cohort.",
        )
        st.session_state["population"] = st.session_state["population_widget"]
        population_selected = st.session_state["population"] != "Choose..."
        if not population_selected:
            st.warning("Please select an LD population before continuing.")

        st.radio(
            "Filter proxies by",
            ["R² only", "D′ only", "R² and D′ (both must pass)"],
            index=["R² only", "D′ only", "R² and D′ (both must pass)"].index(st.session_state["ld_metric"]),
            key="ld_metric_widget",
            help="Which LD statistic(s) a candidate proxy must pass its threshold "
                 "on to be kept. R² measures correlation between alleles; "
                 "D′ measures maximum possible LD.",
        )
        st.session_state["ld_metric"] = st.session_state["ld_metric_widget"]
        use_r2 = st.session_state["ld_metric"] in ("R² only", "R² and D′ (both must pass)")
        use_dprime = st.session_state["ld_metric"] in ("D′ only", "R² and D′ (both must pass)")
        if use_r2:
            st.slider(
                "R² Threshold", 0.0, 1.0, step=0.05,
                value=st.session_state["r2_filter"], key="r2_filter_widget",
                help="Minimum R² (allele correlation) a candidate proxy must have "
                     "with the target SNP to be kept. Higher = stricter, fewer proxies.",
            )
            st.session_state["r2_filter"] = st.session_state["r2_filter_widget"]
        if use_dprime:
            st.slider(
                "D′ Threshold", 0.0, 1.0, step=0.05,
                value=st.session_state["dprime_filter"], key="dprime_filter_widget",
                help="Minimum D′ (normalized LD coefficient) a candidate proxy must "
                     "have with the target SNP to be kept. Higher = stricter, fewer proxies.",
            )
            st.session_state["dprime_filter"] = st.session_state["dprime_filter_widget"]

    st.checkbox(
        "Hide unlinked variants",
        value=st.session_state["proxies_only"],
        key="proxies_only_widget",
        help="Hide variants in the window that are neither the exact target nor a qualifying LD proxy. "
             "Applies to both the on-screen table and the exported CSV.",
    )
    st.session_state["proxies_only"] = st.session_state["proxies_only_widget"]

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        if st.button("Next", type="primary", width='stretch', disabled=not population_selected):
            go_to_next_step()
            st.rerun()


def render_scan_summary() -> None:
    """Print a plain-language summary of every configuration choice before running the scan."""
    s = st.session_state
    start_window = s["target_pos"] - s["window_size"]
    end_window = s["target_pos"] + s["window_size"]

    st.subheader("Summary")
    lines = [
        f"- **File:** `{s['pgs_file_name']}`",
        f"- **Target rsID:** `{s['target_rsid']}` at Chr{s['chromosome']}:{s['target_pos']:,}",
        f"- **Search window:** Chr{s['chromosome']}:{start_window:,} - {end_window:,}",
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

        PROGRESS_CAP = 0.99

        def on_scan_progress(current_pos: int, percent_complete: float, variants_found: int) -> None:
            """
            Callback passed into engine.execute_scan(), invoked periodically during the scan.

            Args:
                current_pos:       Base-pair position currently being read.
                percent_complete:  Fraction (0.0-1.0) of the way through the search window.
                variants_found:    Number of matching variants found so far.

            Displayed position is clamped to the search window because the
            backend's read/parse loop can occasionally report a position
            past end_window (e.g. a final buffered read overshooting the
            window edge).
            """
            clamped = max(0.0, min(PROGRESS_CAP, percent_complete))
            pct = int(clamped * 100)
            display_pos = max(start_window, min(current_pos, end_window))
            progress_bar.progress(
                clamped,
                text=f"Scanning… {pct}% (position {display_pos:,} of Chr{s['chromosome']}:{start_window:,}-{end_window:,})",
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
        progress_bar.progress(1.0, text="Scan complete!")

    status_placeholder.empty()

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
              EXACT MATCH FOUND: {s['target_rsid']} at Chr{s['chromosome']}:{s['target_pos']:,}
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
              PROXY MATCH: target absent at Chr{s['chromosome']}:{s['target_pos']:,}, but
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
    if st.button("Back to Config", width='stretch'):
        go_to_prev_step()
        st.rerun()


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
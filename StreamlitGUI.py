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
  rendering. It is the durable key that the rest of the app (scan summary,
  run_scan, CSV export, and the read-only restatement in Step 5) reads.
- The Roboto typeface is applied entirely through Streamlit's native theming
  ([theme] font = "Roboto" plus [[theme.fontFaces]] in .streamlit/config.toml,
  which serves the files in ./static/ via enableStaticServing=true). Do not
  redeclare @font-face or force font-family via injected CSS below; a
  duplicate declaration previously fought the native theme for control on
  every rerun and made a pre-existing flicker worse.
  Known issue: a brief flash back to the browser's default font on every
  widget interaction is a Streamlit 1.58 frontend behavior, not a bug in
  this file. Streamlit resends the full theme (including fontFaces) as
  part of a new "new_session" message on every script rerun, not just on
  initial page load (confirmed in streamlit/runtime/app_session.py,
  _create_new_session_message, triggered by ScriptRunnerEvent.SCRIPT_STARTED),
  and the frontend re-registers the custom @font-face rules from a freshly
  parsed (but content-identical) array each time, forcing a brief style
  recompute. This is outside this app's control short of patching
  Streamlit's bundled frontend, which is unsupported; it has been accepted
  as a known trade-off of using a self-hosted custom theme font.
"""

import hashlib
import os
import re
import tempfile
import time
import uuid
from datetime import date

import streamlit as st  # type: ignore

import appV3
import constants

RSID_PATTERN = re.compile(r"^rs\d+$", re.IGNORECASE)


def is_valid_rsid(value: str) -> bool:
    """True if value matches the standard NCBI dbSNP rsID format ('rs' + digits)."""
    return bool(RSID_PATTERN.match(value.strip()))


def compute_md5(file_bytes: bytes) -> str:
    """Return the MD5 hex digest of raw file bytes."""
    return hashlib.md5(file_bytes).hexdigest()


def _sync_target_pos_from_text() -> None:
    """
    Callback for the Center Position text field. Strips thousands
    separators and whitespace (positions copied from NCBI's Genome Data
    Viewer are formatted like "27,508,073", which a native numeric input
    silently rejects on paste) and, if what is left is a positive integer,
    writes it into the numeric 'target_pos' session-state key that the
    rest of the app relies on. Invalid input (non-digits, zero, or a
    negative value) is left in the text box as-is and simply does not
    update 'target_pos', so a bad paste cannot silently corrupt the
    search window.

    Reads from the widget's own key ('target_pos_text_widget') rather
    than the durable 'target_pos_text' key, and immediately mirrors the
    raw text into 'target_pos_text' too, per the durable/widget key
    split described in the module docstring.
    """
    raw = st.session_state.get("target_pos_text_widget", "")
    st.session_state["target_pos_text"] = raw
    cleaned = raw.replace(",", "").replace(" ", "").strip()
    if cleaned.isdigit() and int(cleaned) > 0:
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
    which is often GRCh37 even in a file harmonized to GRCh38, so it is
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
    prior file can never be displayed alongside, or exported with, a
    new one.
    """
    st.session_state["scan_results"] = None
    st.session_state["preview_metadata"] = None


# Uploaded score files live on ephemeral disk rather than in session_state.
# Session state is per-user and persists for the whole session, so keeping a
# multi-hundred-MB upload there multiplies with concurrency against the 1 GB a
# free Streamlit Community Cloud app gets; disk is the cheaper resource.
UPLOAD_TMP_DIR = os.path.join(tempfile.gettempdir(), "pgs_grep_uploads")

# How long an untouched upload survives before the startup sweep reclaims it.
UPLOAD_MAX_AGE_SECONDS = 6 * 60 * 60


def sweep_stale_uploads() -> None:
    """
    Delete orphaned uploads left behind by sessions that ended.

    Streamlit exposes no session-end hook, so a closed browser tab strands its
    file. Sweeping on startup bounds how much disk those can occupy. Only files
    untouched for UPLOAD_MAX_AGE_SECONDS are removed, and stored_upload_path()
    refreshes the timestamp on every use, so a file a live session still needs
    is never taken out from under it.
    """
    try:
        entries = list(os.scandir(UPLOAD_TMP_DIR))
    except OSError:
        return  # directory not created yet, or unreadable; nothing to sweep

    cutoff = time.time() - UPLOAD_MAX_AGE_SECONDS
    for entry in entries:
        try:
            if entry.is_file() and entry.stat().st_mtime < cutoff:
                os.remove(entry.path)
        except OSError:
            continue  # already gone, or held open elsewhere; skip it


def store_uploaded_file(file_bytes: bytes) -> str:
    """Write the upload to ephemeral disk and return the path to it."""
    os.makedirs(UPLOAD_TMP_DIR, exist_ok=True)
    path = os.path.join(UPLOAD_TMP_DIR, f"{uuid.uuid4().hex}.txt.gz")
    with open(path, "wb") as fh:
        fh.write(file_bytes)
    return path


def discard_stored_upload() -> None:
    """Delete this session's stored upload, if it still has one."""
    path = st.session_state.get("pgs_file_path")
    if path:
        try:
            os.remove(path)
        except OSError:
            pass
    st.session_state["pgs_file_path"] = None


def stored_upload_path() -> str | None:
    """
    Path to this session's uploaded file, or None if it was never stored or
    has since been swept away.

    Touching the file's mtime on each access marks it as still in use, which
    is what keeps sweep_stale_uploads() from reclaiming a file belonging to a
    session that is simply taking a long time.
    """
    path = st.session_state.get("pgs_file_path")
    if not path or not os.path.exists(path):
        return None
    try:
        os.utime(path, None)
    except OSError:
        pass
    return path


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
        "pgs_file_path": None,
        "pgs_file_name": None,
        "pgs_file_signature": None,
        "pgs_file_is_harmonized": None,
        "preview_metadata": None,
        "want_ld_proxies": "No, scan target position only",
        "ldlink_token": "",
        "target_rsid": "rs1260326",
        "genome_build": "GRCh38",
        "population": "Choose...",
        "chromosome": 2,
        "target_pos": 27_508_073,
        "target_pos_text": "27,508,073",
        "window_size": 0,
        "ld_window_size": 5_000,
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


def active_ld_thresholds() -> tuple[float | None, float | None]:
    """
    Return the (r2_threshold, dprime_threshold) actually in effect, based on
    the "Filter proxies by" choice (Step 5). Only the threshold(s) matching
    that choice are non-None; the other is forced to None regardless of its
    stored slider value, so a stale/unused threshold (e.g. the R² slider's
    value when "D' only" is selected) can never silently get applied to the
    scan or misreported in the CSV output header.
    """
    if not wants_ld_proxies():
        return None, None
    ld_metric = st.session_state["ld_metric"]
    use_r2 = ld_metric in ("R² only", "R² and D′ (both must pass)")
    use_dprime = ld_metric in ("D′ only", "R² and D′ (both must pass)")
    r2_threshold = st.session_state["r2_filter"] if use_r2 else None
    dprime_threshold = st.session_state["dprime_filter"] if use_dprime else None
    return r2_threshold, dprime_threshold


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
    code, pre, .stCode, [data-testid="stCodeBlock"] {
        font-family: 'Roboto Mono', 'Courier New', monospace !important;
    }
    h1, h2, h3, h4 {
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
sweep_stale_uploads()

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
        st.markdown("""**Polygenic Score (PGS)**: A single number summarizing how much
someone's genetics affect their likelihood of having a trait or disease, found from adding up small contributions
from many genetic variants. A PGS file lists those variants and the weight each one carries towards the disease outcome.

**SNP (Single Nucleotide Polymorphism)**: One position in the genome where the DNA
letter (A, T, C, or G) varies between people. PGS Grep locates SNPs by their chromosome and position.

**RSID**: The "rs"-prefixed label assigns to a SNP, such as rs1260326. It is determined by dbSNP, a free public archive.

**LD (Linkage Disequilibrium)**: The tendency of two nearby variants to appear together rather than being inherited independently.

**LD Proxy**: A variant closely linked to your target. When the target itself is absent
from a PGS file, a proxy carries much of the same information.

**LDLink API Access Token**: A free credential from the National Cancer Institute's
[LDLink service](https://ldlink.nih.gov/apiaccess), needed only for proxy searches. It
identifies your requests to LDLink and nothing else. PGS Grep never stores or shares it;
the [source code](https://github.com/CharlieBuhanan/PGS-Grep-Public) shows exactly how it
is used.""")

    st.subheader("Before you start, consider:")
    st.markdown(
        """
- A **target SNP RSID** (ex: rs1260326) that you want to locate
- A Polygenic Score (PGS) file (or a publication on the PGS Catalog)
- Whether or not you want to consider Linkage Disequilibrium in your results (this requires a free [LDLink Token](https://ldlink.nih.gov/apiaccess) from the National Cancer Institute)
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


DEFAULT_MAX_PGS_FILE_SIZE_MB = 200


def max_pgs_file_size_mb() -> int:
    """
    The upload cap, in MB, read back from Streamlit's own server.maxUploadSize.

    Sourcing it here rather than hardcoding it means the in-app check can never
    drift from what the transport layer actually enforces, including when a local
    run raises it with `--server.maxUploadSize=N`.
    """
    try:
        return int(st.config.get_option("server.maxUploadSize"))
    except Exception:
        return DEFAULT_MAX_PGS_FILE_SIZE_MB


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

    max_mb = max_pgs_file_size_mb()

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
   or `_hmPOS_GRCh37.txt.gz`, depending on the genome build.""", help = "GRCh38 (hg38) and GRCh37 (hg19) are two reference assemblies, and a variant sits at a different position in each. Pick the build you will look your rsID up in.")

        st.markdown("5. Download the optional corresponding MD5 file if you want to verify the download's integrity.", help = "A checksum file (.md5 or .txt) that PGS Grep can use to confirm the download arrived intact.")
        st.markdown(f"6. No need to unzip. Upload the `.txt.gz` file directly below and the MD5 file if you have it. {max_mb} megabytes is the maximum size for file uploads.", help=f"PGS Grep reads .gz files as-is. The {max_mb} MB cap is temporary while this site is hosted for free.")

    col_upload, col_md5 = st.columns([1,1])
    with col_upload:
        st.markdown("**Harmonized PGS File** *(`.txt.gz`)*")
        uploaded_file = st.file_uploader(
            "Upload file", type=["gz"], label_visibility="collapsed"
        )
    with col_md5:
        st.markdown("**MD5 Checksum** *(optional)*", help = "A checksum file (.md5 or .txt) that PGS Grep can use to confirm the download arrived intact.")
        uploaded_md5 = st.file_uploader(
            "Upload MD5", type=["md5"], key="md5_uploader",
            label_visibility="collapsed",
        )
    st.caption(
        f"{max_mb} MB maximum file size. Must be in harmonized format."
    )
    st.caption(
        f":grey[The {max_mb} MB limit is temporary while PGS Grep is in testing and "
        "hosted as a free site. Running the app locally lets you raise it.]"
    )

    if uploaded_file is not None:
        if uploaded_file.size > max_mb * 1024 * 1024:
            st.error(
                f"**'{uploaded_file.name}' is "
                f"{uploaded_file.size / (1024 * 1024):.1f} MB, which exceeds the "
                f"{max_mb} MB limit** for this hosted demo. This limit is temporary "
                "while PGS Grep is in testing and hosted as a free site. Please use a "
                "smaller score, or download and run PGS Grep locally (see the note on "
                "the Welcome page), where you can raise the limit."
            )
        else:
            file_bytes = uploaded_file.read()
            file_name = uploaded_file.name

            signature = f"{file_name}:{compute_md5(file_bytes)}"
            if signature != st.session_state["pgs_file_signature"]:
                reset_downstream_state()
                # Drop the previous file from disk before replacing it, so a
                # session that uploads repeatedly leaves only one file behind.
                discard_stored_upload()
                st.session_state["pgs_file_signature"] = signature
                st.session_state["pgs_file_path"] = store_uploaded_file(file_bytes)
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

    if stored_upload_path() is not None:
        st.success(f"📄 File ready: `{st.session_state['pgs_file_name']}`")
        if st.session_state["preview_metadata"] is None:
            try:
                # extract_pgs_metadata takes a path as readily as a file object,
                # so the upload never has to be read back into memory here.
                preview_meta = appV3.extract_pgs_metadata(stored_upload_path())
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
                "**This does not look like a harmonized PGS file.** No harmonized-build "
                "metadata was found, and the filename does not include `_hmPOS_`. "
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
            stored_upload_path() is not None
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
    target rsID, chromosome, center position, and genome build directly
    on this step. Genome build is grouped with the coordinates here (and
    is only editable here) because a chromosome/position pair is only
    meaningful together with the assembly it was looked up in.

    These four fields are bound to durable session_state keys
    ("target_rsid", "chromosome", "target_pos"/"target_pos_text",
    "genome_build") via the widget/durable key split described in the
    module docstring, so their values survive navigating to other steps
    or triggering a script rerun, and are later restated read-only in
    Step 5.
    """
    st.header("Locate Your Target RSID")
    st.markdown(
        "PGS files are indexed by **chromosomal position**, not by rsID, so your "
        "target variant (e.g. `rs1260326`) has to be translated into coordinates "
        "before the scan can run."
    )

    detected_build = get_display_build(st.session_state["preview_metadata"])

    st.subheader("How to find the coordinates")
    st.markdown(
        f"""
1. Open the [NCBI Genome Data Viewer](https://www.ncbi.nlm.nih.gov/gdv), a browser for
   genomic data run by the National Library of Medicine.
2. Search for your rsID (e.g. `rs1260326`).
3. **Set the reference assembly to `{detected_build}` to match your uploaded PGS file.**
   The same variant sits at different coordinates in different assemblies, so a
   mismatch here sends the scan to the wrong place.
4. Copy down the **chromosome number** and **base-pair position** listed for that
   assembly, or leave the tab open, and fill them in below.
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
    st.caption("Enter the rsID, chromosome, and position you looked up above. Default values are shown for rs1260326 as an example.")

    st.text_input(
        "Target rsID (must match center position)",
        value=st.session_state["target_rsid"],
        key="target_rsid_widget",
        help="Format: \"rs\" followed by digits, e.g. rs1260326. Labels your results "
             "and, with LD proxies on, is the variant sent to LDlink.",
    )
    st.session_state["target_rsid"] = st.session_state["target_rsid_widget"]
    rsid_valid = is_valid_rsid(st.session_state["target_rsid"])
    if not rsid_valid:
        st.caption("⚠️ Enter a valid rsID: \"rs\" followed by digits only (e.g. rs1260326)")

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
        help="Base-pair position from the Genome Data Viewer. Commas are fine "
             "(27,508,073 and 27508073 both work). Must be a positive whole number.",
    )
    st.session_state["target_pos_text"] = st.session_state["target_pos_text_widget"]
    cleaned_target_pos = st.session_state["target_pos_text"].replace(",", "").replace(" ", "").strip()
    position_valid = cleaned_target_pos.isdigit() and int(cleaned_target_pos) > 0
    if not position_valid:
        st.caption("⚠️ Enter a positive whole number, digits only (commas are fine)")

    build_options = ["GRCh38", "GRCh37"]
    st.selectbox(
        "Genome Build",
        build_options,
        index=build_options.index(st.session_state["genome_build"]),
        key="genome_build_widget",
        help="The reference assembly your PGS file is harmonized to, auto-detected from its "
             "metadata (defaults to GRCh38 if none found). A mismatch means no proxies will be found. Only change it if the "
             "detected value is wrong.",
    )
    st.session_state["genome_build"] = st.session_state["genome_build_widget"]

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = rsid_valid and position_valid
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
An **LD proxy** is a neighboring SNP that is usually inherited together with your
target, so it carries much of the same information. Turning this on creates a search
window around your target position and reports any variants linked closely to it,
which is useful when the target is missing from the PGS file. Leaving it off
searches only your target position.
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
        "Proxy lookups run through the **LDlink API**, a free service from the "
        "National Cancer Institute. Requests must carry a personal token, which "
        "you can register for in a couple of minutes."
    )
    st.caption(
        "Your token is used only for this session's LDlink queries. It is never stored "
        "or shared, and the "
        "[source code](https://github.com/CharlieBuhanan/PGS-Grep-Public) shows how it "
        "is handled."
    )
    st.markdown("[Get a free token here (ldlink.nih.gov)](https://ldlink.nih.gov/apiaccess)")

    st.text_input(
        "LDLink API Token",
        value=st.session_state["ldlink_token"],
        type="password",
        key="ldlink_token_widget",
        help="LD results are cached on the server, so repeating a query skips the API call.",
    )
    st.session_state["ldlink_token"] = st.session_state["ldlink_token_widget"]
    token_value = st.session_state["ldlink_token"].strip()
    token_has_whitespace = bool(re.search(r"\s", token_value))
    token_valid = bool(token_value) and not token_has_whitespace

    if not token_value:
        st.warning("A token is required to fetch LD proxies. You can leave this blank and go Back to skip LD proxies instead.")
    elif token_has_whitespace:
        st.warning("A token cannot contain spaces. Double-check that you copied it correctly.")

    st.write("")
    col_back, col_next = st.columns(2)
    with col_back:
        if st.button("Back", width='stretch'):
            go_to_prev_step()
            st.rerun()
    with col_next:
        can_advance = token_valid
        if st.button("Next", type="primary", width='stretch', disabled=not can_advance):
            go_to_next_step()
            st.rerun()


def render_step5_config() -> None:
    """
    Restate the target rsID, chromosome, position, and genome build
    collected in Step 3 as read-only fields (genome build is only
    user-editable there, right alongside the coordinates it applies to,
    since a chromosome/position pair is meaningless without knowing which
    assembly it was looked up in), then collect LD proxy settings and
    output filtering preferences. The flanking genomic search window is
    only user-configurable when LD proxies are enabled (Step 4); otherwise
    the scan searches for the single exact target position, and
    'window_size' is forced to 0.
    """
    st.header("Search Configuration")

    st.subheader("Target Variant")
    st.caption("Set in a previous step. Go Back to change these values.")

    st.text_input(
        "Target rsID",
        value=st.session_state["target_rsid"],
        disabled=True,
    )

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

    st.text_input(
        "Genome Build",
        value=st.session_state["genome_build"],
        disabled=True,
    )

    population_selected = True
    if wants_ld_proxies():
        st.subheader("LD Proxy Settings")
        st.number_input(
            "Flanking Size (+/- base pairs)",
            min_value=1,
            step=1_000,
            value=st.session_state["ld_window_size"],
            key="ld_window_size_widget",
            help="How far to either side of the target position the scan looks for "
                 "linked variants. Wider windows catch more proxies but take longer. "
                 "Must be a positive whole number.",
        )
        st.session_state["ld_window_size"] = st.session_state["ld_window_size_widget"]
        st.session_state["window_size"] = st.session_state["ld_window_size"]

        population_options = ["Choose...", "EUR", "AMR", "AFR", "EAS", "SAS"]
        st.selectbox(
            "LD Population (1000 Genomes)",
            population_options,
            index=population_options.index(st.session_state["population"]),
            key="population_widget",
            help="Reference population for the LD calculation. Linkage varies by "
                 "ancestry, so a strong proxy in one group may be unlinked in another. "
                 "Match your study cohort where possible. EUR = European, "
                 "AMR = Admixed American, AFR = African, EAS = East Asian, "
                 "SAS = South Asian.",
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
            help="Which statistic a proxy must clear to be kept. R² reflects how well "
                 "one variant predicts the other; D′ reflects how complete the linkage "
                 "is, ignoring allele frequency.",
        )
        st.session_state["ld_metric"] = st.session_state["ld_metric_widget"]
        use_r2 = st.session_state["ld_metric"] in ("R² only", "R² and D′ (both must pass)")
        use_dprime = st.session_state["ld_metric"] in ("D′ only", "R² and D′ (both must pass)")
        if use_r2:
            st.slider(
                "R² Threshold", 0.0, 1.0, step=0.05,
                value=st.session_state["r2_filter"], key="r2_filter_widget",
                help="Lowest R² a proxy may have with the target. Raise it for fewer, "
                     "tighter-linked proxies.",
            )
            st.session_state["r2_filter"] = st.session_state["r2_filter_widget"]
        if use_dprime:
            st.slider(
                "D′ Threshold", 0.0, 1.0, step=0.05,
                value=st.session_state["dprime_filter"], key="dprime_filter_widget",
                help="Lowest D′ a proxy may have with the target. Raise it for fewer, "
                     "tighter-linked proxies.",
            )
            st.session_state["dprime_filter"] = st.session_state["dprime_filter_widget"]

        start_window = st.session_state["target_pos"] - st.session_state["window_size"]
        end_window = st.session_state["target_pos"] + st.session_state["window_size"]
        st.caption(
            f"Search window: **Chr{st.session_state['chromosome']}:"
            f"{start_window:,} - {end_window:,}**"
        )

        st.checkbox(
            "Hide unlinked variants",
            value=st.session_state["proxies_only"],
            key="proxies_only_widget",
            help="Drop variants that are neither the target nor a qualifying proxy. "
                 "Affects the table and the exported CSV alike.",
        )
        st.session_state["proxies_only"] = st.session_state["proxies_only_widget"]
    else:
        st.session_state["window_size"] = 0
        st.session_state["proxies_only"] = True
        st.caption(
            "LD Proxies are turned off, so PGS Grep searches for the exact "
            "variant only."
        )

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

    search_window_line = (
        f"Chr{s['chromosome']}:{start_window:,} - {end_window:,} (±{s['window_size']:,} bp)"
        if s["window_size"] > 0
        else f"Chr{s['chromosome']}:{s['target_pos']:,} (exact position only, no flanking window)"
    )

    st.subheader("Summary")
    lines = [
        f"- **File:** `{s['pgs_file_name']}`",
        f"- **Target rsID:** `{s['target_rsid']}` at Chr{s['chromosome']}:{s['target_pos']:,}",
        f"- **Search window:** {search_window_line}",
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
        st.caption(f"LD cache found for **{s['target_rsid']}**: Repeated LDlink API calls can be skipped.")
    else:
        st.caption("No LD cache found yet: Will query the LDlink API on scan.")


def run_scan() -> None:
    """
    Execute the genomic scan using the cached PGS file and configuration,
    then store the results in session_state so they survive accidental
    reruns. The progress bar reserves its last stretch for real
    post-scan work rather than claiming completion early: reading the
    file is capped at PROGRESS_CAP (90%), and the remaining 10% covers
    building the filtered export DataFrame and encoding the downloadable
    CSV (both cached in session_state here so render_results() never has
    to redo them on a later rerun). The bar only reaches 1.0 once all of
    that has actually finished.
    """
    s = st.session_state
    start_window = s["target_pos"] - s["window_size"]
    end_window = s["target_pos"] + s["window_size"]

    scan_path = stored_upload_path()
    if scan_path is None:
        # The stored upload expired or was swept while this session sat idle.
        st.error(
            "**Your uploaded file is no longer available.** Uploads are held "
            "only temporarily on the server. Please go back to Step 2 and "
            "upload the score file again."
        )
        return

    engine = appV3.PGSScanEngine(token=s["ldlink_token"])
    status_placeholder = st.empty()

    with status_placeholder.container():
        if wants_ld_proxies():
            r2_threshold, dprime_threshold = active_ld_thresholds()
            with st.spinner("Fetching / loading LD proxy map…"):
                engine.fetch_ld_proxies(
                    target_rsid=s["target_rsid"],
                    genome_build=s["genome_build"],
                    r2_threshold=r2_threshold,
                    dprime_threshold=dprime_threshold,
                    population=s["population"],
                )
            if not engine.ld_map:
                st.warning("No LD proxies returned above the threshold(s). Only the exact target position will be searched.")
            else:
                st.success(f"{len(engine.ld_map)} proxy position(s) loaded")

        st.markdown("**Scanning PGS file…**")
        progress_bar = st.progress(0, text="Starting scan…")
        progress_status = st.empty()

        PROGRESS_CAP = 0.90

        def on_scan_progress(current_pos: int, percent_complete: float, variants_found: int) -> None:
            """
            Callback passed into engine.execute_scan(), invoked periodically during the scan.

            Args:
                current_pos:       Base-pair position currently being read (unused; scans are
                                    fast enough that a bare percentage is all that's worth showing).
                percent_complete:  Fraction (0.0-1.0) of the way through the file.
                variants_found:    Number of matching variants found so far.
            """
            clamped = max(0.0, min(PROGRESS_CAP, percent_complete))
            pct = int(clamped * 100)
            progress_bar.progress(clamped, text=f"Scanning… {pct}%")
            progress_status.caption(f"Variants in window so far: **{variants_found}**")

        exact_match, proxy_matches, results_df = engine.execute_scan(
            file_object=scan_path,
            chr_number=s["chromosome"],
            target_pos=s["target_pos"],
            start_window=start_window,
            end_window=end_window,
            target_rsid=s["target_rsid"],
            progress_callback=on_scan_progress,
        )

        progress_bar.progress(PROGRESS_CAP, text="Formatting results…")

        if s["proxies_only"]:
            export_df = results_df[results_df["Match_Status"].str.contains("PROXY|EXACT", na=False)]
        else:
            export_df = results_df

        r2_threshold, dprime_threshold = active_ld_thresholds()
        csv_credit_line = "# Charlie Buhanan, PGS Grep V1.1, 2026\n"
        header_comments = appV3.build_output_header(
            metadata=engine.last_metadata,
            target_rsid=s["target_rsid"],
            population=s.get("population") if wants_ld_proxies() else None,
            genome_build=s["genome_build"],
            r2_threshold=r2_threshold,
            dprime_threshold=dprime_threshold,
        )
        csv_data = export_df.to_csv(index=False)
        csv_bytes = (csv_credit_line + header_comments + csv_data).encode("utf-8")

        progress_bar.progress(1.0, text="Scan complete!")

    status_placeholder.empty()

    s["scan_results"] = {
        "exact_match": exact_match,
        "proxy_matches": proxy_matches,
        "results_df": results_df,
        "export_df": export_df,
        "csv_bytes": csv_bytes,
        "last_metadata": engine.last_metadata,
    }


def render_results() -> None:
    """Render the cached scan results: match banner, metrics, data table, and CSV download."""
    s = st.session_state
    results = s["scan_results"]
    exact_match = results["exact_match"]
    proxy_matches = results["proxy_matches"]
    results_df = results["results_df"]
    export_df = results["export_df"]
    csv_bytes = results["csv_bytes"]
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
              PROXY MATCH: Target absent at Chr{s['chromosome']}:{s['target_pos']:,}, but
              {len(proxy_matches)} LD-linked proxy variant(s) found in the window.
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.error(
            f"No match: Neither the target position (Chr{s['chromosome']}:{s['target_pos']:,}) "
            "nor any qualifying proxies were found in this window."
        )

    m1, m2, _ = st.columns([1, 1, 3])
    m1.metric("Exact Position Match", "Yes" if exact_match else "No")
    m2.metric("LD Proxies Captured", len(proxy_matches))

    if results_df.empty:
        st.warning("No variants fell within the specified genomic window.")
        return

    st.subheader("📋 Data Viewer")
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
        """Row-color a Match_Status cell based on whether it is an exact, proxy, or unlinked call."""
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
    st.caption("Most scans finish in about 10 seconds.")

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
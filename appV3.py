# TODO: GUI GUI GUI GUI GUI UI UX
# TODO: Some Test-Driven Development!

# Make a list of 12 one-line edge case test scenarios.
# Change some button colors. Image somewhere! Some more light blues. D' can be negative. Tooltip for LD Cache.
# TODO: Make sure README is accurate. Help message for setup LD API token, user side / release public
# TODO: Extend test cases. Testing suite: make sure this works for all harmonized formats including edge cases (missing columns, missing rsIDs, etc.)
# TODO: Read through all important code, check test suite. 
# TODO: Update & Check test suite. Test README file instructions, edit & revise README. Update requirements.txt if necessary
# TODO: Make this more efficient. Better clear cache button
# Future: option to see all LD proxies (view cache file basically)

import gzip
import io
import os
from datetime import date
from typing import Final
import streamlit as st # type: ignore
import pandas as pd # type: ignore
import requests # type: ignore
import constants

# LDlink uses numeric codes 23/24/25 for X/Y/MT; harmonized PGS files use letters.
LDLINK_SEX_MT_CODES: Final[dict[str, str]] = {"23": "X", "24": "Y", "25": "MT"}
PGS_SEX_MT_CODES: Final[dict[str, str]] = {"X": "23", "Y": "24", "MT": "25"}


def normalize_chr(chr_value: str) -> str:
    """
    Normalize a chromosome value to PGS harmonized notation (1-22, X, Y, MT).

    Accepts either LDlink-style numeric codes (23/24/25) or PGS-style letters
    (X/Y/MT, case-insensitive), as well as plain autosomes. Strips a leading
    'chr' prefix if present.
    """
    c = str(chr_value).strip().replace("chr", "").replace("Chr", "")
    upper = c.upper()
    if upper in ("X", "Y", "MT"):
        return upper
    if c in LDLINK_SEX_MT_CODES:
        return LDLINK_SEX_MT_CODES[c]
    return c


# Canonical output columns (always emitted in this order)
OUTPUT_COLS: Final[list[str]] = [
    "chr_name",
    "chr_position",
    "RS_ID",
    "effect_allele",
    "other_allele",
    "effect_weight",
    "Match_Status",
]

def _find_col(headers: list[str], candidates: list[str]) -> int | None:
    """
    Return the index of the first header that matches any candidate name
    (case-insensitive).  Returns None if not found.
    """
    lower_headers = [h.lower().strip() for h in headers]
    for cand in candidates:
        try:
            return lower_headers.index(cand.lower())
        except ValueError:
            continue
    return None


def resolve_column_map(headers: list[str]) -> dict:
    """
    Given a list of header strings from a PGS file, return a dict mapping
    logical field names → column indices.  Raises ValueError for required
    fields that cannot be found.

    Supports harmonized PGS formats as documented at https://github.com/PGScatalog/pgs-harmonizer
    """
    col = {}

    # ── chromosome ────────────────────────────────────────────────────────
    col["chr"] = _find_col(headers, ["hm_chr", "chr_name"])
    if col["chr"] is None:
        raise ValueError("Cannot find chromosome column (expected 'hm_chr' or 'chr_name').")

    # ── position ──────────────────────────────────────────────────────────
    col["pos"] = _find_col(headers, ["hm_pos", "chr_position"])
    if col["pos"] is None:
        raise ValueError("Cannot find position column (expected 'hm_pos' or 'chr_position').")

    # ── rsID (optional – may be absent or '.' for some rows) ──────────────
    col["rsid"] = _find_col(headers, ["hm_rsid", "rsid"])  # None is OK

    # ── alleles ───────────────────────────────────────────────────────────
    col["eff_al"] = _find_col(headers, ["effect_allele"])
    col["oth_al"] = _find_col(headers, ["other_allele", "reference_allele",
                                         "hm_inferotherallele"])

    # ── weight ────────────────────────────────────────────────────────────
    col["weight"] = _find_col(headers, ["effect_weight", "weight"])

    return col


def extract_pgs_metadata(file_object) -> dict:
    """
    Extract metadata from comment lines at the top of a harmonized PGS file.

    PGS Catalog harmonized files begin with '#'-prefixed comment lines in the form:
        #key=value
    This function parses those lines and returns a dict of key-value pairs,
    plus a 'date_accessed' key set to today's date.

    Args:
        file_object: File path (str) or file-like object (gzipped PGS file).

    Returns:
        Dict with metadata fields found in the file header, plus 'date_accessed'.
    """
    metadata: dict = {}

    if isinstance(file_object, str):
        opener = gzip.open(file_object, "rt", encoding="utf-8")
    else:
        file_object.seek(0)
        opener = gzip.open(io.BytesIO(file_object.read()), "rt", encoding="utf-8")

    with opener as f:
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if not line.startswith("#"):
                break  # comment block is over
            # Strip leading '#' and optional whitespace
            content = line.lstrip("#").strip()
            if "=" in content:
                key, _, value = content.partition("=")
                metadata[key.strip()] = value.strip()

    metadata["date_accessed"] = date.today().isoformat()
    return metadata


def build_output_header(
    metadata: dict,
    target_rsid: str,
    population: str,
    genome_build: str,
    r2_threshold: float | None = None,
    dprime_threshold: float | None = None,
) -> str:
    """
    Build a comment block for the top of the CSV output file.

    Args:
        metadata:          Dict returned by extract_pgs_metadata().
        target_rsid:       The queried rsID.
        population:        LD population code.
        genome_build:      Genome assembly string.
        r2_threshold:      R² linkage filter applied, or None if not used.
        dprime_threshold:  D′ linkage filter applied, or None if not used.

    Returns:
        A multi-line string of '#'-prefixed comment lines.
    """
    # Build human-readable LD filter description
    if r2_threshold is not None and dprime_threshold is not None:
        ld_filter_str = f"r^2 ≥ {r2_threshold} AND D' ≥ {dprime_threshold}"
    elif r2_threshold is not None:
        ld_filter_str = f"r^2 ≥ {r2_threshold}"
    elif dprime_threshold is not None:
        ld_filter_str = f"D' ≥ {dprime_threshold}"
    else:
        ld_filter_str = "None"

    lines = [
        "# PGS: Scan Output",
        "#",
        "# === Source File Metadata ===",
        f"# PGS ID: {metadata.get('pgs_id', metadata.get('PGS ID', 'Not found'))}",
        f"# PGS Name: {metadata.get('pgs_name', metadata.get('trait_mapped', 'Not found'))}",
        f"# Trait (EFO): {metadata.get('trait_efo', metadata.get('trait_efo_id', 'Not found'))}",
        f"# Genome Build: {metadata.get('genome_build', metadata.get('HmPOS_build', 'Not found'))}",
        f"# Original Author: {metadata.get('citation', metadata.get('Citation', 'Not found'))}",
        f"# Publication: {metadata.get('pgp_id', metadata.get('PGP ID', 'Not found'))}",
        f"# License: {metadata.get('license', 'Not found')}",
        f"# Date Accessed: {metadata.get('date_accessed', 'N/A')}",
        "#",
        "# === Query Parameters ===",
        f"# Target rsID: {target_rsid}",
        f"# LD Population: {population}",
        f"# Genome Assembly: {genome_build}",
        f"# LD Filter: {ld_filter_str}",
        "#",
    ]
    return "\n".join(lines) + "\n"


class PGSScanEngine:
    """ Engine for scanning PGS (Polygenic Score) files with LDlink LD proxy integration.

    Manages LD proxy data fetching/caching and scans harmonized PGS score files
    for target and proxy variant matches within specified genomic windows.
    """
    
    def __init__(self, token: str):
        """Initialize the PGS scan engine.
        Args:
            token: LDlink API authentication token for LD proxy queries.
        """
        self.token = token
        self.ld_map: dict[int, dict] = {}
        self.last_metadata: dict = {}
        os.makedirs(constants.LD_CACHE_DIR, exist_ok=True)

    def _cache_path(self, rsid: str, population: str, genome_build: str) -> str:
        """
        Generate a local cache file path for LD proxy data, as to not repeat API calls.
        """
        safe_build = genome_build.replace(" ", "_")
        filename = f"LD{rsid}_{population}_{safe_build}.txt"
        return os.path.join(constants.LD_CACHE_DIR, filename)

    def _parse_ld_text(
        self,
        text: str,
        r2_threshold: float | None,
        dprime_threshold: float | None,
    ) -> dict[int, dict]:
        """
        Parse LDlink API response text into a map of proxies.

        Filters rows by whichever combination of R² / D′ thresholds are active.
        Both thresholds must pass when both are provided.
        Raw D′ values from LDlink can be negative (indicating repulsion phase LD);
        we compare against the absolute value, which is standard practice.
        """
        ld_map: dict[int, dict] = {}
        if "RS_Number" not in text:
            return ld_map

        lines = text.strip().split("\n")
        headers = lines[0].split()

        try:
            rs_idx  = headers.index("RS_Number")
            pos_idx = headers.index("Coord")
            r2_idx  = headers.index("R2")
        except ValueError:
            return ld_map

        # D′ column is optional, gracefully absent in some LDlink responses
        try:
            dprime_idx: int | None = headers.index("Dprime")
        except ValueError:
            dprime_idx = None

        for line in lines[1:]:
            cols = line.split()
            if len(cols) <= max(rs_idx, pos_idx, r2_idx):
                continue
            try:
                r2_val = float(cols[r2_idx])

                dprime_val: float | None = None
                if dprime_idx is not None and dprime_idx < len(cols):
                    try:
                        dprime_val = abs(float(cols[dprime_idx]))
                    except ValueError:
                        dprime_val = None

                # Apply R² filter
                if r2_threshold is not None and r2_val < r2_threshold:
                    continue

                # Apply D′ filter (skip row if D′ is required but unavailable)
                if dprime_threshold is not None:
                    if dprime_val is None or dprime_val < dprime_threshold:
                        continue

                pos = int(cols[pos_idx].split(":")[-1])
                ld_map[pos] = {
                    "rsid":   cols[rs_idx],
                    "r2":     r2_val,
                    "dprime": dprime_val,  # None when column absent
                }
            except (ValueError, IndexError):
                continue
        return ld_map

    def fetch_ld_proxies(
        self,
        target_rsid:      str,
        genome_build:     str,
        population:       str = "EUR",
        r2_threshold:     float | None = None,
        dprime_threshold: float | None = None,
    ) -> dict[int, dict]:
        """
        Fetch LD proxy variants for a target rsID using LDlink API (or cache).

        Thresholds are applied after fetching. The raw API response is always
        cached unfiltered so threshold changes don't require a new API call.

        Args:
            target_rsid:      rsID to query (e.g. "rs10305420").
            genome_build:     "GRCh38" or "GRCh37".
            population:       1000 Genomes population code (e.g. "EUR").
            r2_threshold:     Minimum R² to include a proxy, or None to skip R² filter.
            dprime_threshold: Minimum |D′| to include a proxy, or None to skip D′ filter.
        """
        self.ld_map = {}
        cache_file  = self._cache_path(target_rsid, population, genome_build)

        if os.path.exists(cache_file):
            st.info(f"📂 Loading LD data from cache: `{cache_file}`")
            with open(cache_file, "r", encoding="utf-8") as f:
                cached_text = f.read()
            self.ld_map = self._parse_ld_text(cached_text, r2_threshold, dprime_threshold)
            return self.ld_map

        if not self.token:
            st.error("❌ No LDlink API token provided and no cached result found.")
            return self.ld_map

        build_param = genome_build.lower().strip()
        url = (
            f"https://ldlink.nih.gov/LDlinkRest/ldproxy?"
            f"var={target_rsid}&pop={population}&r2_d=r2&"
            f"genome_build={build_param}&token={self.token}"
        )

        try:
            response = requests.get(url, timeout=30)
        except Exception as e:
            st.error(f"❌ Network error contacting LDlink: {e}")
            return self.ld_map

        if response.status_code != 200:
            st.error(f"❌ LDlink server error: HTTP {response.status_code}")
            return self.ld_map

        text = response.text

        if text.strip().startswith("<!DOCTYPE") or "<html" in text.lower():
            st.error(
                "❌ LDlink returned an HTML error page — check your token, "
                "network proxy, or target rsID."
            )
            with open("error_debug_log.html", "w", encoding="utf-8") as f:
                f.write(text)
            return self.ld_map

        if "error" in text.lower() or "API token invalid" in text:
            st.error(f"❌ LDlink API error: {text.strip()}")
            return self.ld_map

        # Cache the raw response unfiltered — threshold changes reuse this file
        with open(cache_file, "w", encoding="utf-8") as f:
            f.write(text)
        st.success(f"✅ LD response saved to cache: `{cache_file}`")

        self.ld_map = self._parse_ld_text(text, r2_threshold, dprime_threshold)
        return self.ld_map

    def execute_scan(
        self,
        file_object,
        chr_number:        int,
        target_pos:        int,
        start_window:      int,
        end_window:        int,
        target_rsid:       str,
        progress_callback  = None,
    ) -> tuple[bool, list, pd.DataFrame]:
        """
        Scan a gzipped PGS file for target and LD proxy variants in a genomic window.
        Also populates self.last_metadata with the file's header metadata.

        Args:
            file_object:        File path (str) or Streamlit UploadedFile (.txt.gz).
            chr_number:         Chromosome to restrict scan to.
            target_pos:         Exact base-pair position of the target variant.
            start_window:       Start of the genomic search window.
            end_window:         End of the genomic search window.
            target_rsid:        rsID of the target variant (used as fallback label).
            progress_callback:  Optional callable(current_pos, percent_complete, variants_found).
                                Called every PROGRESS_INTERVAL data lines on the
                                target chromosome. percent_complete is the SNP's
                                fractional position within [start_window, end_window],
                                clamped to [0.0, 1.0], so the bar reflects how far
                                through the search window the scan actually is
                                rather than how many lines have been read.
        """
        PROGRESS_INTERVAL = 500

        exact_match:    bool = False
        proxy_matches:  list = []
        rows_processed: list = []
        lines_read:     int  = 0
        target_chr_str = normalize_chr(chr_number)

        # Extract metadata before scanning (resets seek position internally)
        self.last_metadata = extract_pgs_metadata(file_object)

        ld_rsid_by_pos: dict[int, str] = {
            pos: info["rsid"] for pos, info in self.ld_map.items()
        }

        if isinstance(file_object, str):
            opener = gzip.open(file_object, "rt", encoding="utf-8")
        else:
            file_object.seek(0)
            opener = gzip.open(io.BytesIO(file_object.read()), "rt", encoding="utf-8")

        col_map: dict | None = None

        with opener as f:
            for raw_line in f:
                line = raw_line.rstrip("\n")

                if line.startswith("#") or not line.strip():
                    continue

                columns = line.split("\t") if "\t" in line else line.split()

                if col_map is None:
                    first = columns[0].lower().strip()
                    if first in ("rsid", "chr_name", "hm_chr", "#chr_name"):
                        try:
                            col_map = resolve_column_map(columns)
                        except ValueError as e:
                            st.error(f"❌ Header parse error: {e}")
                            return False, [], pd.DataFrame()
                        continue
                    else:
                        continue

                lines_read += 1

                def safe_get(idx, default="N/A"):
                    if idx is None or idx >= len(columns):
                        return default
                    v = columns[idx].strip()
                    return v if v not in ("", ".", "NA", "nan", "NaN", "Null") else default

                c_name = normalize_chr(safe_get(col_map["chr"], ""))
                c_pos  = safe_get(col_map["pos"], "")

                try:
                    current_pos = int(c_pos)
                except ValueError:
                    continue

                if c_name != target_chr_str:
                    continue

                if progress_callback and lines_read % PROGRESS_INTERVAL == 0:
                    window_span = end_window - start_window
                    if window_span > 0:
                        pct = (current_pos - start_window) / window_span
                    else:
                        pct = 1.0
                    pct = max(0.0, min(1.0, pct))
                    progress_callback(current_pos, pct, len(rows_processed))

                if not (start_window <= current_pos <= end_window):
                    continue

                rs_id = safe_get(col_map.get("rsid"))
                if rs_id == "N/A" and current_pos in ld_rsid_by_pos:
                    rs_id = ld_rsid_by_pos[current_pos]

                if rs_id == "N/A" and current_pos == target_pos:
                    rs_id = target_rsid

                eff_al = safe_get(col_map.get("eff_al"))
                oth_al = safe_get(col_map.get("oth_al"))
                weight = safe_get(col_map.get("weight"), "0.0")

                if current_pos == target_pos:
                    status_flag = "🎯 EXACT TARGET MATCH"
                    exact_match = True
                elif current_pos in self.ld_map:
                    p = self.ld_map[current_pos]
                    # Include D′ in the label when it was fetched
                    if p.get("dprime") is not None:
                        status_flag = (
                            f"🔗 LD PROXY ({p['rsid']}, r²={p['r2']:.3f}, D′={p['dprime']:.3f})"
                        )
                    else:
                        status_flag = f"🔗 LD PROXY ({p['rsid']}, r²={p['r2']:.3f})"
                    proxy_matches.append((current_pos, p["rsid"], p["r2"], weight))
                else:
                    status_flag = "Unlinked Region Variant"

                rows_processed.append({
                    "chr_name":      c_name,
                    "chr_position":  current_pos,
                    "RS_ID":         rs_id,
                    "effect_allele": eff_al,
                    "other_allele":  oth_al,
                    "effect_weight": weight,
                    "Match_Status":  status_flag,
                })

        # Final callback so the GUI always reaches 100 %
        if progress_callback:
            progress_callback(end_window, 1.0, len(rows_processed))

        df = pd.DataFrame(rows_processed, columns=OUTPUT_COLS) if rows_processed else pd.DataFrame(columns=OUTPUT_COLS)
        return exact_match, proxy_matches, df
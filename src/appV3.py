
# TODO: Extend test cases
# TODO: Read through, check test suite. Testing suite: make sure this works for all harmonized formats including edge cases (missing columns, missing rsIDs, etc.)
# TODO: UI/UX as in notes, test interactions to make sure it won't crash with spam. Ease of use, little question marks everywhere
# TODO: Figure out chromosome 23 handling (X/Y/MT) in LDlink and harmonized PGS files.  LDlink uses 23 for X, 24 for Y, 25 for MT.  Harmonized PGS files use chr_name = "X", "Y", "MT".  Need to map these correctly in the scan.
# TODO: Update & Check test suite. Test README file instructions, edit & revise README.
# TODO: Make this more efficient. Better clear cache button
# TODO: Clone into public repo with deleted branch history to remove API token from history

import gzip
import io
import os
from datetime import date
from typing import Final
import streamlit as st # type: ignore
import pandas as pd # type: ignore
import requests # type: ignore
import constants

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


def build_output_header(metadata: dict, target_rsid: str, population: str,
                         genome_build: str, r2_threshold: float) -> str:
    """
    Build a comment block for the top of the CSV output file.

    Args:
        metadata:      Dict returned by extract_pgs_metadata().
        target_rsid:   The queried rsID.
        population:    LD population code.
        genome_build:  Genome assembly string.
        r2_threshold:  R² linkage filter applied.

    Returns:
        A multi-line string of '#'-prefixed comment lines.
    """
    lines = [
        "# LD-Aware PGS Grepper — Scan Output",
        "#",
        "# === Source File Metadata ===",
        f"# PGS ID           : {metadata.get('pgs_id', metadata.get('PGS ID', 'N/A'))}",
        f"# PGS Name         : {metadata.get('pgs_name', metadata.get('trait_mapped', 'N/A'))}",
        f"# Trait (EFO)      : {metadata.get('trait_efo', metadata.get('trait_efo_id', 'N/A'))}",
        f"# Genome Build     : {metadata.get('genome_build', metadata.get('HmPOS_build', genome_build))}",
        f"# Original Author  : {metadata.get('citation', metadata.get('Citation', 'N/A'))}",
        f"# Publication      : {metadata.get('pgp_id', metadata.get('PGP ID', 'N/A'))}",
        f"# License           : {metadata.get('license', 'N/A')}",
        f"# Date Accessed    : {metadata.get('date_accessed', 'N/A')}",
        "#",
        "# === Query Parameters ===",
        f"# Target rsID      : {target_rsid}",
        f"# LD Population    : {population}",
        f"# Genome Assembly  : {genome_build}",
        f"# R² Threshold     : {r2_threshold}",
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

    def _parse_ld_text(self, text: str, r2_threshold: float) -> dict[int, dict]:
        """
        Parse LDlink API response text into a map of proxies.
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

        for line in lines[1:]:
            cols = line.split()
            if len(cols) <= max(rs_idx, pos_idx, r2_idx):
                continue
            try:
                r2_val = float(cols[r2_idx])
                if r2_val >= r2_threshold:
                    pos = int(cols[pos_idx].split(":")[-1])
                    ld_map[pos] = {
                        "rsid": cols[rs_idx],
                        "r2":   r2_val,
                    }
            except (ValueError, IndexError):
                continue
        return ld_map

    def fetch_ld_proxies(
        self,
        target_rsid:   str,
        genome_build:  str,
        r2_threshold:  float,
        population:    str = "EUR",
    ) -> dict[int, dict]:
        """
        Fetch LD proxy variants for a target rsID using LDlink API (or cache).
        """
        self.ld_map = {}
        cache_file  = self._cache_path(target_rsid, population, genome_build)

        if os.path.exists(cache_file):
            st.info(f"📂 Loading LD data from cache: `{cache_file}`")
            with open(cache_file, "r", encoding="utf-8") as f:
                cached_text = f.read()
            self.ld_map = self._parse_ld_text(cached_text, r2_threshold)
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

        with open(cache_file, "w", encoding="utf-8") as f:
            f.write(text)
        st.success(f"✅ LD response saved to cache: `{cache_file}`")

        self.ld_map = self._parse_ld_text(text, r2_threshold)
        return self.ld_map

    def execute_scan(
        self,
        file_object,
        chr_number:   int,
        target_pos:   int,
        start_window: int,
        end_window:   int,
        target_rsid: str 
    ) -> tuple[bool, list, pd.DataFrame]:
        """
        Scan a gzipped PGS file for target and LD proxy variants in a genomic window.
        Also populates self.last_metadata with the file's header metadata.
        """

        exact_match:   bool  = False
        proxy_matches: list  = []
        rows_processed: list = []
        target_chr_str = str(chr_number).replace("chr", "")

        # Extract metadata before scanning
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

                def safe_get(idx, default="N/A"):
                    if idx is None or idx >= len(columns):
                        return default
                    v = columns[idx].strip()
                    return v if v not in ("", ".", "NA", "nan", "NaN", "Null") else default

                c_name = safe_get(col_map["chr"], "").replace("chr", "")
                c_pos  = safe_get(col_map["pos"], "")

                try:
                    current_pos = int(c_pos)
                except ValueError:
                    continue

                if c_name != target_chr_str:
                    continue
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

        df = pd.DataFrame(rows_processed, columns=OUTPUT_COLS) if rows_processed else pd.DataFrame(columns=OUTPUT_COLS)
        return exact_match, proxy_matches, df
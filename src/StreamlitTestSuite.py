"""
StreamlitTestSuite.py
====================
pytest suite for appV3.PGSScanEngine and resolve_column_map.
All file I/O uses in-memory gzip streams; all HTTP calls are mocked.

Install Dependencies:
    pip install pytest pytest-mock requests

Run Tests:
    python -m pytest test_pgs_explorer.py -v
"""

import gzip
import io
import os
import textwrap

import pandas as pd
import pytest
import requests

# appV3 calls streamlit at import time; stub it before importing.
import sys
from unittest.mock import MagicMock, patch

sys.modules.setdefault("streamlit", MagicMock())

import appV3  # noqa: E402 — must come after the st stub

# Shared fixtures
LDLINK_RESPONSE = textwrap.dedent("""\
    RS_Number\tCoord\tAlleles\tMA\tMAF\tR2\tDprime
    rs7903146\tchr10:112998590\tC/T\tT\t0.23\t1.0\t1.0
    rs4506565\tchr10:113005746\tA/T\tT\t0.23\t0.97\t0.99
    rs12255372\tchr10:113037807\tG/T\tT\t0.23\t0.95\t0.99
    rs7901695\tchr10:112993858\tC/T\tT\t0.23\t0.75\t0.97
    rs11196205\tchr10:113012271\tC/T\tT\t0.23\t0.62\t0.94
""")

PGS_TSV = textwrap.dedent("""\
    # PGS Catalog Score File
    chr_name\tchr_position\trsID\teffect_allele\tother_allele\teffect_weight
    10\t112998590\trs7903146\tT\tC\t0.42
    10\t113005746\trs4506565\tT\tA\t0.31
    10\t113037807\trs12255372\tT\tG\t0.18
    10\t999999999\trs0000000\tA\tG\t0.05
""")


def make_gz(content: str) -> io.BytesIO:
    """Return an in-memory gzipped BytesIO of *content*."""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content.encode())
    buf.seek(0)
    return buf


@pytest.fixture()
def engine(tmp_path, monkeypatch):
    """PGSScanEngine with its cache dir redirected to a temp folder."""
    monkeypatch.setattr(appV3.constants, "LD_CACHE_DIR", str(tmp_path))
    return appV3.PGSScanEngine(token="valid-token-123")


@pytest.fixture()
def engine_no_token(tmp_path, monkeypatch):
    monkeypatch.setattr(appV3.constants, "LD_CACHE_DIR", str(tmp_path))
    return appV3.PGSScanEngine(token="")


@pytest.fixture()
def mock_200(monkeypatch):
    """Patch requests.get to return a valid 200 LDlink response."""
    resp = MagicMock()
    resp.status_code = 200
    resp.text = LDLINK_RESPONSE
    monkeypatch.setattr(requests, "get", MagicMock(return_value=resp))
    return requests.get


@pytest.fixture()
def mock_401(monkeypatch):
    resp = MagicMock()
    resp.status_code = 401
    resp.text = "Unauthorized"
    monkeypatch.setattr(requests, "get", MagicMock(return_value=resp))
    return requests.get


@pytest.fixture()
def mock_429(monkeypatch):
    resp = MagicMock()
    resp.status_code = 429
    resp.text = "Too Many Requests"
    monkeypatch.setattr(requests, "get", MagicMock(return_value=resp))
    return requests.get


@pytest.fixture()
def pgs_file():
    return make_gz(PGS_TSV)

# TC-01: Empty API token → no HTTP call, ld_map stays empty
def test_empty_token_skips_api_call(engine_no_token, mock_200):
    """fetch_ld_proxies returns an empty map and never calls requests.get."""
    result = engine_no_token.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.7,
        population="EUR",
    )

    assert result == {}
    mock_200.assert_not_called()


# ─────────────────────────────────────────────────────────────────────────────
# TC-02: HTTP 401 → ld_map stays empty, no exception raised to caller
# ─────────────────────────────────────────────────────────────────────────────
def test_http_401_returns_empty_map(engine, mock_401):
    """A 401 response causes fetch_ld_proxies to return {} without crashing."""
    result = engine.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.7,
        population="EUR",
    )

    assert result == {}
    assert engine.ld_map == {}


# ─────────────────────────────────────────────────────────────────────────────
# TC-03: Valid 200 response → _parse_ld_text produces correct ld_map
# ─────────────────────────────────────────────────────────────────────────────

def test_valid_response_populates_ld_map(engine, mock_200):
    """200 response text is parsed into a position-keyed ld_map."""
    result = engine.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.0,
        population="EUR",
    )

    # All five rows from LDLINK_RESPONSE should be present
    assert len(result) == 5
    # Spot-check a known entry
    assert 112998590 in result
    assert result[112998590]["rsid"] == "rs7903146"
    assert result[112998590]["r2"] == 1.0


def test_parse_ld_text_direct(engine):
    """Unit-test _parse_ld_text independently of HTTP."""
    ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0)

    assert isinstance(ld_map, dict)
    assert 113005746 in ld_map
    assert ld_map[113005746]["rsid"] == "rs4506565"


# ─────────────────────────────────────────────────────────────────────────────
# TC-04: Malformed RSID → LDlink returns error text → ld_map empty
# (appV3 guards on the API's own error text; no client-side regex)
# ─────────────────────────────────────────────────────────────────────────────

def test_malformed_rsid_error_text_returns_empty_map(engine, monkeypatch):
    """When LDlink returns an error message body, ld_map stays empty."""
    resp = MagicMock()
    resp.status_code = 200
    resp.text = "error: variant not found in 1000 Genomes"
    monkeypatch.setattr(requests, "get", MagicMock(return_value=resp))

    result = engine.fetch_ld_proxies(
        target_rsid="MYSNP123",
        genome_build="GRCh38",
        r2_threshold=0.7,
        population="EUR",
    )

    assert result == {}


def test_parse_ld_text_missing_header_returns_empty(engine):
    """_parse_ld_text gracefully returns {} when 'RS_Number' is absent."""
    bad_text = "Col1\tCol2\nfoo\tbar\n"
    result = engine._parse_ld_text(bad_text, r2_threshold=0.0)
    assert result == {}


# ─────────────────────────────────────────────────────────────────────────────
# TC-05: r² threshold filters rows correctly inside _parse_ld_text
# ─────────────────────────────────────────────────────────────────────────────
def test_r2_filter_08_keeps_three_rows(engine):
    # R2 values: 1.0, 0.97, 0.95, 0.75, 0.62 → at ≥0.8: first three only
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.8)
    assert len(result) == 3
    assert all(v["r2"] >= 0.8 for v in result.values())


def test_r2_filter_10_keeps_one_row(engine):
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=1.0)
    assert len(result) == 1
    assert list(result.values())[0]["r2"] == 1.0


def test_r2_filter_00_keeps_all_rows(engine):
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0)
    assert len(result) == 5


# ─────────────────────────────────────────────────────────────────────────────
# TC-06: File-based caching → second call reads from disk, not the API
# ─────────────────────────────────────────────────────────────────────────────

def test_cache_hit_skips_network(engine, mock_200, tmp_path, monkeypatch):
    """Pre-seed the cache file; requests.get must never be called."""
    monkeypatch.setattr(appV3.constants, "LD_CACHE_DIR", str(tmp_path))
    cache_file = engine._cache_path("rs7903146", "EUR", "GRCh38")
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    with open(cache_file, "w") as f:
        f.write(LDLINK_RESPONSE)

    result = engine.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.0,
        population="EUR",
    )

    mock_200.assert_not_called()
    assert len(result) == 5


def test_two_live_calls_write_then_read_cache(engine, mock_200, tmp_path, monkeypatch):
    """First call hits the API and writes cache; second call reads it."""
    monkeypatch.setattr(appV3.constants, "LD_CACHE_DIR", str(tmp_path))

    engine.fetch_ld_proxies(
        target_rsid="rs7903146", genome_build="GRCh38",
        r2_threshold=0.0, population="EUR",
    )
    engine.fetch_ld_proxies(
        target_rsid="rs7903146", genome_build="GRCh38",
        r2_threshold=0.0, population="EUR",
    )

    # requests.get called exactly once across both invocations
    assert mock_200.call_count == 1


# ─────────────────────────────────────────────────────────────────────────────
# TC-07: Population parameter is embedded correctly in the API URL
# ─────────────────────────────────────────────────────────────────────────────

def test_population_eur_in_url(engine, mock_200):
    engine.fetch_ld_proxies(
        target_rsid="rs7903146", genome_build="GRCh38",
        r2_threshold=0.7, population="EUR",
    )
    called_url = mock_200.call_args[0][0]
    assert "pop=EUR" in called_url


def test_population_afr_in_url(engine, mock_200):
    engine.fetch_ld_proxies(
        target_rsid="rs7903146", genome_build="GRCh38",
        r2_threshold=0.7, population="AFR",
    )
    called_url = mock_200.call_args[0][0]
    assert "pop=AFR" in called_url


# ─────────────────────────────────────────────────────────────────────────────
# TC-08: execute_scan returns a valid CSV-serialisable DataFrame
# ─────────────────────────────────────────────────────────────────────────────

def test_execute_scan_returns_dataframe(engine, pgs_file):
    # Pre-load ld_map so we don't need a network call
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.7)

    _, _, df = engine.execute_scan(
        file_object=pgs_file,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert isinstance(df, pd.DataFrame)
    assert list(df.columns) == appV3.OUTPUT_COLS


def test_results_df_serialises_to_csv(engine, pgs_file):
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.7)

    _, _, df = engine.execute_scan(
        file_object=pgs_file,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    csv_bytes = df.to_csv(index=False).encode("utf-8")
    assert isinstance(csv_bytes, bytes)
    roundtrip = pd.read_csv(io.StringIO(csv_bytes.decode()))
    assert list(roundtrip.columns) == appV3.OUTPUT_COLS


# ─────────────────────────────────────────────────────────────────────────────
# TC-09: Unrecognised file (non-gzip) → execute_scan raises an exception
# ─────────────────────────────────────────────────────────────────────────────

def test_pdf_bytes_raise_on_scan(engine):
    """A plain non-gzip stream causes gzip.open to raise inside execute_scan."""
    fake_pdf = io.BytesIO(b"%PDF-1.4 this is not gzip data at all")

    with pytest.raises(Exception):
        engine.execute_scan(
            file_object=fake_pdf,
            chr_number=10,
            target_pos=112998590,
            start_window=112990000,
            end_window=113050000,
            target_rsid="rs7903146",
        )


def test_plain_text_gz_with_no_tsv_structure_yields_empty_df(engine):
    """Gzipped content with no recognisable header → empty DataFrame."""
    junk_gz = make_gz("this,is,a,csv,not,a,pgs,file\n1,2,3,4,5,6,7\n")

    _, _, df = engine.execute_scan(
        file_object=junk_gz,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert df.empty


# ─────────────────────────────────────────────────────────────────────────────
# TC-10: HTTP 429 → ld_map stays empty, no exception raised to caller
# ─────────────────────────────────────────────────────────────────────────────

def test_http_429_returns_empty_map(engine, mock_429):
    """A 429 Too Many Requests is handled without crashing."""
    result = engine.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.7,
        population="EUR",
    )

    assert result == {}
    assert engine.ld_map == {}


def test_http_429_does_not_write_cache(engine, mock_429, tmp_path, monkeypatch):
    """No cache file should be written on a 429."""
    monkeypatch.setattr(appV3.constants, "LD_CACHE_DIR", str(tmp_path))

    engine.fetch_ld_proxies(
        target_rsid="rs7903146",
        genome_build="GRCh38",
        r2_threshold=0.7,
        population="EUR",
    )

    cache_file = engine._cache_path("rs7903146", "EUR", "GRCh38")
    assert not os.path.exists(cache_file)


# ─────────────────────────────────────────────────────────────────────────────
# resolve_column_map unit tests
# ─────────────────────────────────────────────────────────────────────────────

def test_resolve_column_map_standard_headers():
    headers = ["chr_name", "chr_position", "rsID", "effect_allele",
               "other_allele", "effect_weight"]
    col = appV3.resolve_column_map(headers)
    assert col["chr"] == 0
    assert col["pos"] == 1
    assert col["rsid"] == 2


def test_resolve_column_map_harmonized_headers():
    headers = ["hm_chr", "hm_pos", "hm_rsID", "effect_allele", "other_allele", "effect_weight"]
    col = appV3.resolve_column_map(headers)
    assert col["chr"] == 0
    assert col["pos"] == 1
    assert col["rsid"] == 2


def test_resolve_column_map_missing_chr_raises():
    with pytest.raises(ValueError, match="chromosome"):
        appV3.resolve_column_map(["chr_position", "effect_allele"])


def test_resolve_column_map_missing_pos_raises():
    with pytest.raises(ValueError, match="position"):
        appV3.resolve_column_map(["chr_name", "effect_allele"])


# ─────────────────────────────────────────────────────────────────────────────
# execute_scan match-classification tests
# ─────────────────────────────────────────────────────────────────────────────

def test_exact_match_detected(engine, pgs_file):
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0)

    exact, proxies, df = engine.execute_scan(
        file_object=pgs_file,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert exact is True
    assert any("EXACT" in str(s) for s in df["Match_Status"])


def test_proxy_matches_detected(engine, pgs_file):
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0)

    exact, proxies, df = engine.execute_scan(
        file_object=pgs_file,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert len(proxies) >= 1
    assert any("PROXY" in str(s) for s in df["Match_Status"])


def test_out_of_window_variants_excluded(engine, pgs_file):
    """rs0000000 at position 999999999 must not appear in results."""
    engine.ld_map = {}

    _, _, df = engine.execute_scan(
        file_object=pgs_file,
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert 999999999 not in df["chr_position"].values


def test_wrong_chromosome_excluded(engine):
    content = textwrap.dedent("""\
        chr_name\tchr_position\trsID\teffect_allele\tother_allele\teffect_weight
        9\t112998590\trs9999999\tT\tC\t0.10
    """)
    engine.ld_map = {}

    _, _, df = engine.execute_scan(
        file_object=make_gz(content),
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    assert df.empty
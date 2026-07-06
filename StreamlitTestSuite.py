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
    ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0, dprime_threshold=None)

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
    result = engine._parse_ld_text(bad_text, r2_threshold=0.0, dprime_threshold=None)
    assert result == {}


# ─────────────────────────────────────────────────────────────────────────────
# TC-05: r² threshold filters rows correctly inside _parse_ld_text
# ─────────────────────────────────────────────────────────────────────────────
def test_r2_filter_08_keeps_three_rows(engine):
    # R2 values: 1.0, 0.97, 0.95, 0.75, 0.62 → at ≥0.8: first three only
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.8, dprime_threshold=None)
    assert len(result) == 3
    assert all(v["r2"] >= 0.8 for v in result.values())


def test_r2_filter_10_keeps_one_row(engine):
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=1.0, dprime_threshold=None)
    assert len(result) == 1
    assert list(result.values())[0]["r2"] == 1.0


def test_r2_filter_00_keeps_all_rows(engine):
    result = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0, dprime_threshold=None)
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
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.7, dprime_threshold=None)

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
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.7, dprime_threshold=None)

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
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0, dprime_threshold=None)

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
    engine.ld_map = engine._parse_ld_text(LDLINK_RESPONSE, r2_threshold=0.0, dprime_threshold=None)

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
# ═════════════════════════════════════════════════════════════════════════
# StreamlitGUI.py (wizard UI) tests
# ═════════════════════════════════════════════════════════════════════════
#
# appV3.py is the scan engine (already tested above); StreamlitGUI.py is the
# actual multi-step wizard UI built on top of it. Its steps aren't
# st.tabs() — they're a hand-rolled wizard driven by st.session_state
# ("wizard_step"), advanced/reversed via go_to_next_step() / go_to_prev_step().
# That's the direct analogue of "switching tabs and back" in this app: each
# render_stepN_* function is what gets drawn on a given rerun, and only
# values written into session_state's durable keys are expected to survive
# moving away from a step and back.
#
# streamlit.testing.v1.AppTest (Streamlit's official headless UI test tool)
# is deliberately NOT used here: as of this Streamlit release it has no
# support for simulating st.file_uploader, which is central to two of the
# requested scenarios (non-harmonized file, oversized file). Instead these
# tests call the render_stepN_* functions directly against a stubbed
# `streamlit` module (same technique the existing suite already uses for
# appV3), with session_state replaced by a real dict and st.button/
# st.columns/st.file_uploader given controllable fakes.

import streamlit as _st_stub  # same MagicMock stub already in sys.modules

# StreamlitGUI executes wizard-rendering code at import time (its bottom line
# unconditionally calls STEP_RENDERERS[st.session_state["wizard_step"]]()), so
# session_state must be a real dict *before* the module is imported, or that
# first render pass KeyErrors against a MagicMock instead of a real string.
_st_stub.session_state = {}

import StreamlitGUI as gui


class FakeUploadedFile:
    """Minimal stand-in for Streamlit's UploadedFile: bytes + name + size."""

    def __init__(self, content: bytes, name: str, size: int | None = None):
        self._buf = io.BytesIO(content)
        self.name = name
        self.size = size if size is not None else len(content)

    def read(self, *a, **kw):
        return self._buf.read(*a, **kw)

    def seek(self, *a, **kw):
        return self._buf.seek(*a, **kw)


def _fake_columns(spec, **kwargs):
    """st.columns(2) and st.columns([1,1]) both need to be unpackable."""
    n = spec if isinstance(spec, int) else len(spec)
    return tuple(MagicMock() for _ in range(n))


@pytest.fixture()
def gui_state():
    """
    Fresh, isolated session_state + a controllable `st` mock for each
    StreamlitGUI test. Resets call history so assertions on st.error /
    st.warning / st.button etc. only see calls made by the test itself,
    and gives every test a clean wizard state via init_session_state().
    """
    gui.st.reset_mock()
    gui.st.session_state = {}
    gui.st.columns = MagicMock(side_effect=_fake_columns)
    gui.st.button = MagicMock(return_value=False)
    gui.st.file_uploader = MagicMock(return_value=None)
    gui.init_session_state()
    return gui.st.session_state


def _button_call(label):
    """Find the recorded st.button(...) call for a given label."""
    return next(c for c in gui.st.button.call_args_list if c.args[0] == label)


def make_harmonized_gz() -> bytes:
    content = textwrap.dedent("""\
        ###PGS CATALOG SCORING FILE
        #pgs_id=PGS000014
        #hmpos_build=GRCh38
        hm_chr\thm_pos\thm_rsID\teffect_allele\tother_allele\teffect_weight
        10\t112998590\trs7903146\tT\tC\t0.42
        """)
    return _gz_bytes(content)


def make_non_harmonized_gz() -> bytes:
    content = textwrap.dedent("""\
        # PGS Catalog Score File
        chr_name\tchr_position\trsID\teffect_allele\tother_allele\teffect_weight
        10\t112998590\trs7903146\tT\tC\t0.42
        """)
    return _gz_bytes(content)


def _gz_bytes(content: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write(content.encode())
    return buf.getvalue()


# ── Non-harmonized upload: must error and block progress ───────────────────

def test_non_harmonized_upload_shows_error_and_blocks_next(gui_state):
    """Uploading a non-harmonized file (no hm_pos metadata, filename lacks
    '_hmPOS_') must surface a UI error and disable the Next button."""
    raw = make_non_harmonized_gz()
    gui.st.file_uploader = MagicMock(
        return_value=FakeUploadedFile(raw, "PGS000014.txt.gz")
    )

    gui.render_step2_upload()

    assert gui_state["pgs_file_is_harmonized"] is False
    error_text = " ".join(str(c.args[0]) for c in gui.st.error.call_args_list)
    assert "doesn't look like a harmonized" in error_text
    assert _button_call("Next").kwargs["disabled"] is True


def test_harmonized_upload_has_no_error_and_enables_next(gui_state):
    """A properly harmonized file must not trigger the error and must leave
    Next enabled, as a contrasting control case for the test above."""
    raw = make_harmonized_gz()
    gui.st.file_uploader = MagicMock(
        return_value=FakeUploadedFile(raw, "PGS000014_hmPOS_GRCh38.txt.gz")
    )

    gui.render_step2_upload()

    assert gui_state["pgs_file_is_harmonized"] is True
    assert gui.st.error.call_args_list == []
    assert _button_call("Next").kwargs["disabled"] is False


# ── Oversized upload: must error and block progress ─────────────────────────

def test_oversized_upload_shows_error_and_blocks_next(gui_state):
    """A file over the 200 MB cap must be rejected with a clear error,
    never stored in session_state, and must leave Next disabled."""
    oversized = FakeUploadedFile(
        b"tiny content",
        "PGS000014_hmPOS_GRCh38.txt.gz",
        size=gui.MAX_PGS_FILE_SIZE_BYTES + 1,
    )
    gui.st.file_uploader = MagicMock(return_value=oversized)

    gui.render_step2_upload()

    assert gui_state["pgs_file_bytes"] is None
    error_text = " ".join(str(c.args[0]) for c in gui.st.error.call_args_list)
    assert "exceeds the 200 MB limit" in error_text
    assert _button_call("Next").kwargs["disabled"] is True


def test_file_exactly_at_size_limit_is_accepted(gui_state):
    """Boundary check: the cap is a strict '>' comparison, so a file of
    exactly MAX_PGS_FILE_SIZE_BYTES must be accepted, not rejected."""
    raw = make_harmonized_gz()
    exact = FakeUploadedFile(
        raw, "PGS000014_hmPOS_GRCh38.txt.gz", size=gui.MAX_PGS_FILE_SIZE_BYTES
    )
    gui.st.file_uploader = MagicMock(return_value=exact)

    gui.render_step2_upload()

    assert gui_state["pgs_file_bytes"] == raw
    assert gui.st.error.call_args_list == []


# ── Session state survives moving back and forth between wizard steps ──────

def test_session_state_persists_across_step_navigation(gui_state):
    """Values collected on one step must still be there after navigating
    forward and back through other steps — the wizard's equivalent of
    switching tabs and coming back."""
    gui_state["wizard_step"] = "rsid_guide"
    gui_state["target_rsid"] = "rs123456"
    gui_state["chromosome"] = 7
    gui_state["target_pos"] = 555555

    gui.go_to_next_step()  # -> ld_choice
    gui.go_to_next_step()  # -> ld_auth or config, depending on want_ld_proxies
    gui.go_to_prev_step()
    gui.go_to_prev_step()  # back to rsid_guide

    assert gui_state["wizard_step"] == "rsid_guide"
    assert gui_state["target_rsid"] == "rs123456"
    assert gui_state["chromosome"] == 7
    assert gui_state["target_pos"] == 555555


def test_reset_downstream_state_only_clears_scan_fields(gui_state):
    """Detecting a newly-uploaded file clears stale scan results/metadata,
    but must not wipe unrelated session state (target rsID, thresholds,
    etc.) that the user already entered on other steps."""
    gui_state["scan_results"] = {"exact_match": True}
    gui_state["preview_metadata"] = {"pgs_id": "PGS1"}
    gui_state["target_rsid"] = "rs999"
    gui_state["r2_filter"] = 0.9

    gui.reset_downstream_state()

    assert gui_state["scan_results"] is None
    assert gui_state["preview_metadata"] is None
    assert gui_state["target_rsid"] == "rs999"
    assert gui_state["r2_filter"] == 0.9


# ── "Next step" behavior around the LDLink API key choice ──────────────────

def test_next_step_skips_token_page_when_ld_proxies_declined(gui_state):
    """Choosing 'No' to LD proxies must skip straight from ld_choice to
    config — the token page should never be shown."""
    gui_state["wizard_step"] = "ld_choice"
    gui_state["want_ld_proxies"] = "No, scan target position only"

    gui.go_to_next_step()

    assert gui_state["wizard_step"] == "config"


def test_next_step_visits_token_page_when_ld_proxies_requested(gui_state):
    """Choosing 'Yes' to LD proxies must route through ld_auth (the token
    entry page) before reaching config."""
    gui_state["wizard_step"] = "ld_choice"
    gui_state["want_ld_proxies"] = "Yes, also search LD proxies (requires a free token)"

    gui.go_to_next_step()

    assert gui_state["wizard_step"] == "ld_auth"


def test_prev_step_skips_token_page_when_ld_proxies_declined(gui_state):
    """The same skip logic must apply symmetrically going backward."""
    gui_state["wizard_step"] = "config"
    gui_state["want_ld_proxies"] = "No, scan target position only"

    gui.go_to_prev_step()

    assert gui_state["wizard_step"] == "ld_choice"


def test_token_page_next_disabled_until_token_entered(gui_state):
    """The API-token step must not let the user advance with a blank
    token, regardless of which LD-proxies choice got them there."""
    gui_state["wizard_step"] = "ld_auth"
    gui_state["ldlink_token"] = ""
    gui_state["ldlink_token_widget"] = ""

    gui.render_step4_5_ld_auth()

    assert _button_call("Next").kwargs["disabled"] is True
    warning_text = " ".join(str(c.args[0]) for c in gui.st.warning.call_args_list)
    assert "token is required" in warning_text


def test_token_page_next_enabled_once_token_entered(gui_state):
    """Once a non-blank token is present, Next must become enabled."""
    gui_state["wizard_step"] = "ld_auth"
    gui_state["ldlink_token"] = "my-real-token"
    gui_state["ldlink_token_widget"] = "my-real-token"

    gui.render_step4_5_ld_auth()

    assert _button_call("Next").kwargs["disabled"] is False
    assert gui.st.warning.call_args_list == []


# ═════════════════════════════════════════════════════════════════════════
# Additional appV3.py edge cases
# ═════════════════════════════════════════════════════════════════════════

def test_dprime_only_filter(engine):
    """dprime_threshold alone (r2_threshold=None) must filter purely on D',
    ignoring R² entirely."""
    # Dprime column values in LDLINK_RESPONSE: 1.0, 0.99, 0.99, 0.97, 0.94
    result = engine._parse_ld_text(
        LDLINK_RESPONSE, r2_threshold=None, dprime_threshold=0.98
    )
    assert len(result) == 3
    assert all(v["dprime"] >= 0.98 for v in result.values())


def test_combined_r2_and_dprime_filter_is_stricter_than_either_alone(engine):
    """When both thresholds are supplied, a row must pass BOTH to be kept."""
    r2_only = engine._parse_ld_text(
        LDLINK_RESPONSE, r2_threshold=0.9, dprime_threshold=None
    )
    combined = engine._parse_ld_text(
        LDLINK_RESPONSE, r2_threshold=0.9, dprime_threshold=0.995
    )
    assert len(r2_only) == 3
    assert len(combined) == 1
    assert list(combined.values())[0]["dprime"] == 1.0


def test_negative_dprime_is_compared_by_absolute_value(engine):
    """Raw LDlink D' can be negative (repulsion phase); filtering must use
    |D'|, not the signed value."""
    text = textwrap.dedent("""\
        RS_Number\tCoord\tAlleles\tMA\tMAF\tR2\tDprime
        rs0000001\tchr1:1000\tC/T\tT\t0.2\t0.9\t-0.85
    """)
    kept = engine._parse_ld_text(text, r2_threshold=None, dprime_threshold=0.8)
    excluded = engine._parse_ld_text(text, r2_threshold=None, dprime_threshold=0.9)
    assert len(kept) == 1
    assert kept[1000]["dprime"] == 0.85
    assert len(excluded) == 0


def test_html_error_page_response_is_handled_gracefully(engine, monkeypatch, tmp_path):
    """An HTML error page (e.g. from a proxy/WAF) instead of TSV must not
    crash fetch_ld_proxies, must not populate ld_map, and must not write
    a cache file."""
    monkeypatch.chdir(tmp_path)
    resp = MagicMock()
    resp.status_code = 200
    resp.text = "<!DOCTYPE html><html><body>Access Denied</body></html>"
    monkeypatch.setattr(requests, "get", MagicMock(return_value=resp))

    result = engine.fetch_ld_proxies(
        target_rsid="rs7903146", genome_build="GRCh38",
        r2_threshold=0.7, population="EUR",
    )

    assert result == {}
    cache_file = engine._cache_path("rs7903146", "EUR", "GRCh38")
    assert not os.path.exists(cache_file)


def test_normalize_chr_ldlink_numeric_sex_codes(engine):
    """LDlink's numeric 23/24/25 codes must map to PGS harmonized X/Y/MT."""
    assert appV3.normalize_chr("23") == "X"
    assert appV3.normalize_chr("24") == "Y"
    assert appV3.normalize_chr("25") == "MT"


def test_normalize_chr_strips_prefix_and_normalizes_case(engine):
    """A leading 'chr'/'Chr' prefix is stripped, and sex/MT letters are
    normalized to uppercase regardless of input case."""
    assert appV3.normalize_chr("chr7") == "7"
    assert appV3.normalize_chr("Chr10") == "10"
    assert appV3.normalize_chr("x") == "X"
    assert appV3.normalize_chr("mt") == "MT"


def test_extract_pgs_metadata_normalizes_keys(engine):
    """Comment-line keys are lowercased and spaces become underscores, so
    downstream lookups (e.g. metadata.get('trait_mapped')) work regardless
    of how the source file capitalized/spaced its header."""
    content = textwrap.dedent("""\
        #PGS ID=PGS000014
        #Trait Mapped=Type 2 Diabetes
        hm_chr\thm_pos\thm_rsID\teffect_allele\tother_allele\teffect_weight
        10\t112998590\trs7903146\tT\tC\t0.42
    """)
    meta = appV3.extract_pgs_metadata(make_gz(content))
    assert meta["pgs_id"] == "PGS000014"
    assert meta["trait_mapped"] == "Type 2 Diabetes"


def test_resolve_column_map_alternate_other_allele_headers(engine):
    """'other_allele' has two accepted alternates; both must resolve."""
    headers = ["hm_chr", "hm_pos", "hm_rsID", "effect_allele",
               "reference_allele", "effect_weight"]
    col = appV3.resolve_column_map(headers)
    assert col["oth_al"] == 4

    headers2 = ["hm_chr", "hm_pos", "hm_rsID", "effect_allele",
                "hm_inferotherallele", "effect_weight"]
    col2 = appV3.resolve_column_map(headers2)
    assert col2["oth_al"] == 4


def test_missing_rsid_column_falls_back_to_ld_map_and_target_rsid(engine):
    """When the PGS file has no rsID column at all, RS_ID must fall back to
    the LD-proxy rsID for proxy positions, and to target_rsid for the exact
    target position."""
    content = textwrap.dedent("""\
        hm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight
        10\t112998590\tT\tC\t0.42
        10\t113005746\tT\tA\t0.31
    """)
    engine.ld_map = engine._parse_ld_text(
        LDLINK_RESPONSE, r2_threshold=0.0, dprime_threshold=None
    )

    _, _, df = engine.execute_scan(
        file_object=make_gz(content),
        chr_number=10,
        target_pos=112998590,
        start_window=112990000,
        end_window=113050000,
        target_rsid="rs7903146",
    )

    row_target = df[df["chr_position"] == 112998590].iloc[0]
    row_proxy = df[df["chr_position"] == 113005746].iloc[0]
    assert row_target["RS_ID"] == "rs7903146"
    assert row_proxy["RS_ID"] == "rs4506565"
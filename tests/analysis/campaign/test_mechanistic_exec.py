"""Mechanistic execution layer: mdfit.xtc policy + provenance + preflight."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from analysis.campaign import mechanistic_exec as mx


def test_registry_shape():
    assert mx.PROFILE == "xanthone_mechanistic"
    assert {s.name for s in mx.SYSTEMS} == {
        "A1", "APO", "COA", "COMP", "A3-HMG-R", "HMG-R-200ns-A6",
        "system_A3_COA", "system_A6_COA"}
    trio = [s for s in mx.SYSTEMS if s.comparison_label in ("A1", "A3", "A6")]
    assert {s.name for s in trio} == set(mx.SHARED_BASIS_TRIO)
    # references are the ones that need reprocessing
    assert all(s.needs_reprocess for s in mx.SYSTEMS if s.role == "reference")
    assert not any(s.needs_reprocess for s in mx.SYSTEMS if s.role == "new")


def test_preflight_blocks_when_nothing_analysable(tmp_path):
    d = tmp_path / "sys"
    d.mkdir()
    sd = mx.SystemDef("X", str(d), "new", "HMG-CoA reductase", "LIG", "xanthone",
                      None, "ActiveSite_HMG", "Catalytic_HMG", None, False)
    c = mx.preflight(sd, gmx="/nonexistent-gmx")
    assert c.status == "REVIEW_REQUIRED"
    assert any("no mdfit.xtc / mdcenter.xtc / md.xtc" in r for r in c.reasons)


def test_preflight_flags_missing_site_and_calpha(tmp_path):
    d = tmp_path / "sys"
    d.mkdir()
    (d / "md.tpr").write_bytes(b"tpr")
    (d / "mdfit.xtc").write_bytes(b"fit")
    (d / "index.ndx").write_text("[ System ]\n1\n[ Protein ]\n1 2\n[ LIG ]\n3\n")
    sd = mx.SystemDef("X", str(d), "new", "HMG-CoA reductase", "LIG", "xanthone",
                      None, "ActiveSite_HMG", "Catalytic_HMG", None, False)
    c = mx.preflight(sd, gmx="/nonexistent-gmx")
    assert any("Catalytic_HMG" in r or "ActiveSite_HMG" in r for r in c.reasons)
    assert any("C-alpha" in r for r in c.reasons)
    assert c.status == "REVIEW_REQUIRED"          # no C-alpha -> hard block
    assert c.fit_source.endswith("mdfit.xtc")


def test_provenance_is_mdfit_and_records_raw(tmp_path):
    src = tmp_path / "src"
    (src).mkdir()
    fit = src / "mdfit.xtc"; fit.write_bytes(b"fit")
    raw = src / "md.xtc"; raw.write_bytes(b"raw")
    tpr = src / "md.tpr"; tpr.write_bytes(b"tpr")
    ndx = src / "index.ndx"; ndx.write_text("[ Protein ]\n1\n")
    out = tmp_path / "o"
    out.mkdir()
    o1 = out / "rmsd.xvg"
    o1.write_text("@\n0 0.1\n")
    prov = mx.write_provenance(
        out / "provenance.protein_rmsd.json", analysis="protein_rmsd", system="A1",
        argv=["gmx", "rms", "-f", str(fit)], stdin="Backbone\nBackbone\n",
        source_trajectory=str(fit), tpr=str(tpr), index=str(ndx),
        semantic_groups={"protein": "Backbone"}, gmx="gmx", returncode=0,
        outputs=[o1], raw_md_xtc=str(raw))
    assert prov["source_trajectory"].endswith("mdfit.xtc")
    assert prov["raw_production_trajectory"].endswith("md.xtc")
    assert prov["preprocessing_chain"][0] == "md.xtc"
    assert prov["source_fingerprints"]["analysis_trajectory"]["sha256"]
    assert prov["source_fingerprints"]["raw_production_trajectory"]["sha256"]
    on_disk = json.loads((out / "provenance.protein_rmsd.json").read_text())
    assert on_disk["profile"] == "xanthone_mechanistic"


def test_shell_preview_roundtrip():
    s = mx.shell_preview(["gmx", "rms", "-f", "a b.xtc"], "Backbone\n")
    assert "a b.xtc" in s and s.startswith("printf %s")


def test_real_sources_preflight_ok_when_present():
    """If the scientific source trees are on this machine, preflight must resolve
    a Cα group.  Uses gmx='' to skip the (slow) trajectory span probe."""
    sd = mx.SYSTEM_BY_NAME["A1"]
    if not Path(sd.directory).is_dir():
        pytest.skip("scientific source tree not on this machine")
    c = mx.preflight(sd, gmx="/nonexistent-gmx")   # no span probe
    assert c.c_alpha == 1614
    assert c.ligand_atoms == 35
    assert c.fit_source in ("reprocess",) or c.fit_source.endswith("mdfit.xtc")
    assert "md.xtc is not an acceptable" not in " ".join(c.reasons)

import json
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = ROOT / "src" / "gff_stats.py"


def load_module():
    spec = importlib.util.spec_from_file_location("gff_stats", str(SRC_PATH))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_compute_stats_from_sample():
    mod = load_module()
    sample = ROOT / "data" / "sample.gff"
    stats = mod.compute_stats_from_gff(str(sample))

    assert stats["total_features"] == 6
    assert stats["by_type"] == {"gene": 2, "CDS": 3, "mRNA": 1}
    # averages: gene lengths 900 and 900 -> 900.0; CDS lengths 101,51,101 -> 84.3
    assert stats["avg_length"]["gene"] == 900.0
    assert stats["avg_length"]["CDS"] == 84.3
    assert stats["avg_length"]["mRNA"] == 150.0
    # strand distribution: + = 4, - = 2 -> 67 and 33 (rounded)
    assert stats["strand_distribution"]["+"] in (66, 67, 67)
    assert stats["strand_distribution"]["-"] in (33, 34, 33)


def test_compute_stats_with_filter():
    mod = load_module()
    sample = ROOT / "data" / "sample.gff"
    
    # Filter by CDS type
    stats = mod.compute_stats_from_gff(str(sample), filter_type="CDS")
    assert stats["total_features"] == 3
    assert stats["by_type"] == {"CDS": 3}
    assert stats["avg_length"]["CDS"] == 84.3
    assert stats["filter_type"] == "CDS"
    
    # Filter by gene type
    stats = mod.compute_stats_from_gff(str(sample), filter_type="gene")
    assert stats["total_features"] == 2
    assert stats["by_type"] == {"gene": 2}
    assert stats["avg_length"]["gene"] == 900.0
    assert stats["filter_type"] == "gene"


def test_cli_writes_json(tmp_path):
    mod = load_module()
    sample = ROOT / "data" / "sample.gff"
    out = tmp_path / "out.json"
    # call cli with argv
    res = mod.cli([f"--gff={str(sample)}", f"--out={str(out)}"])
    assert res == 0
    assert out.exists()
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["total_features"] == 6
    assert "by_type" in data and "avg_length" in data and "strand_distribution" in data


def test_cli_with_filter(tmp_path):
    mod = load_module()
    sample = ROOT / "data" / "sample.gff"
    out = tmp_path / "out_filtered.json"
    # call cli with filter
    res = mod.cli([f"--gff={str(sample)}", f"--out={str(out)}", "--filter-type=CDS"])
    assert res == 0
    assert out.exists()
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["total_features"] == 3
    assert data["by_type"] == {"CDS": 3}
    assert data["filter_type"] == "CDS"

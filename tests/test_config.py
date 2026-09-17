from pathlib import Path

from chartseq_tools.config import load_config, write_example


def test_config_merge_and_example(tmp_path: Path):
    example = write_example(tmp_path / "chartseq.yaml")
    assert "bc_pattern" in example.read_text(encoding="utf-8")
    example.write_text("pipeline:\n  threads: 3\n", encoding="utf-8")
    config = load_config(example)
    assert config["pipeline"]["threads"] == 3
    assert config["library"]["bc_pattern"].startswith("C")


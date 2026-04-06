from pathlib import Path

from src.config import load_config


def test_load_config_resolves_relative_paths(tmp_path: Path) -> None:
    config_dir = tmp_path / "config"
    config_dir.mkdir()
    config_file = config_dir / "local.yaml"
    config_file.write_text(
        "\n".join(
            [
                "project_name: demo",
                "paths:",
                "  metadata_file: ../inputs/meta.csv",
                "  query_cif_dir: ../inputs/cifs",
                "  foldseek_bin: ../bin/foldseek",
                "  target_database: ../db/test_db",
                "  foldseek_output_tsv: ../outputs/results.tsv",
                "  tmp_dir: ../tmp/work",
                "run:",
                "  top_hits: 5",
            ]
        ),
        encoding="utf-8",
    )

    config = load_config(config_file)

    assert config.project_name == "demo"
    assert config.run.top_hits == 5
    assert config.paths.metadata_file == (tmp_path / "inputs" / "meta.csv").resolve()

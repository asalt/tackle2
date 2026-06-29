import tomllib
from pathlib import Path

from click.testing import CliRunner

from python.cli import APP_NAME, main


def test_get_config_creates_default_files():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(main, ["get-config"])
        assert result.exit_code == 0
        assert "Writing" in result.output
        assert "done" in result.output
        config_path = Path(f"{APP_NAME}.toml")
        assert config_path.exists()
        config = tomllib.loads(config_path.read_text())
        assert config["params"]["bubbleplot"]["advanced"]["height_scale"] == 0.8


def test_get_config_with_colormap():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(
            main,
            [
                "get-config",
                "--include-colormap",
                "--name",
                "custom.toml",
                "--colormap-name",
                "custom-map.json",
            ],
        )
        assert result.exit_code == 0
        assert "custom.toml" in result.output
        assert "custom-map.json" in result.output
        assert Path("custom.toml").exists()
        assert Path("custom-map.json").exists()

import subprocess
import sys
import tempfile
from pathlib import Path


def test_package_install_and_import():
    repo_root = Path(__file__).resolve().parent.parent  # One level up from tests/
    dist_dir = tempfile.mkdtemp()

    # Step 1: Build the wheel
    subprocess.check_call(
        [sys.executable, "-m", "build", "--wheel", "--outdir", dist_dir],
        cwd=repo_root,
    )

    wheel_files = list(Path(dist_dir).glob("*.whl"))
    assert wheel_files, "No wheel file was built."
    wheel_path = wheel_files[0]

    # Step 2: Create a virtual environment
    venv_dir = Path(dist_dir) / "venv"
    subprocess.check_call([sys.executable, "-m", "venv", str(venv_dir)])

    pip = venv_dir / "bin" / "pip"
    python = venv_dir / "bin" / "python"

    # Step 3: Install the wheel
    subprocess.check_call([str(pip), "install", str(wheel_path)])

    # Step 4: Test import of installed package
    subprocess.check_call([str(python), "-c", "import cellxgene_schema"])

    # Step 5 (optional): Test CLI entry point
    subprocess.check_call([str(python), "-m", "cellxgene_schema.cli", "--help"])

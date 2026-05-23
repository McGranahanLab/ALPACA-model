import subprocess
from pathlib import Path

import pytest

from golden_utils import GUROBI_VOLATILE_COLS, VOLATILE_COLS, assert_csv_matches_golden  # noqa: F401


def pytest_addoption(parser):
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Overwrite golden files with the current test outputs instead of comparing.",
    )


@pytest.fixture(scope="session")
def update_golden(request):
    return request.config.getoption("--update-golden")


@pytest.fixture(scope="session")
def repo_root():
    return Path(__file__).parent.parent


@pytest.fixture
def run_alpaca(repo_root):
    def _run(*args, expect_failure=False):
        result = subprocess.run(
            ["alpaca", *args],
            capture_output=True,
            text=True,
            cwd=repo_root,
        )
        if not expect_failure:
            assert result.returncode == 0, (
                f"alpaca exited with code {result.returncode}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )
        return result

    return _run

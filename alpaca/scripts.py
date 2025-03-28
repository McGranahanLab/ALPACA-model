import os
import subprocess
import sys
import pkg_resources


def input_conversion():
    """
    Wrapper function to execute input_conversion.sh from submodule
    """
    print("Running input_conversion - it may take a few minutes")
    try:
        # Locate the input_conversion.sh script
        script_path = pkg_resources.resource_filename(
            "alpaca", "scripts/submodules/alpaca_input_formatting/input_conversion.sh"
        )

        # Locate the submodules directory
        submodules_path = pkg_resources.resource_filename(
            "alpaca", "scripts/submodules"
        )

        # Ensure the script is executable
        os.chmod(script_path, 0o755)

        # Set environment variable to help script locate its dependencies
        env = os.environ.copy()
        env["SUBMODULES_PATH"] = submodules_path

        # Execute the shell script with all passed arguments
        result = subprocess.run(
            [script_path] + sys.argv[1:],
            check=True,
            text=True,
            capture_output=True,
            env=env,
        )
        print(result.stdout)
        return 0
    except subprocess.CalledProcessError as e:
        print(f"Error executing input_conversion: {e.stderr}", file=sys.stderr)
        return e.returncode
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(input_conversion())

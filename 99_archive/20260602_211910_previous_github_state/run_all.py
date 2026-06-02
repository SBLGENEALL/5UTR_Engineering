import subprocess
import sys


def main():
    cmd = [
        sys.executable,
        "run_from_step.py",
        "--start",
        "0",
        "--end",
        "16",
        "--config",
        "config/config.yaml",
    ]

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
SRC = HERE / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from organpath.cli import phyview_main


if __name__ == "__main__":
    raise SystemExit(phyview_main())

# Damit die Tests die Skripte als Paket "skripte" finden.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

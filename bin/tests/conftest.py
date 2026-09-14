import sys
from pathlib import Path

# vcfstats.py / combine_stats.py live in bin/, one directory up from here, and
# aren't part of an installed package -- make them importable as plain modules.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

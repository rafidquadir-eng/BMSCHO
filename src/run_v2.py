"""Phase 2, one command: python src/run_v2.py

Needs the phase-1 downloads (python cho_poc.py --download) and phase-1 results
(committed in results/). Stages, in order:
  1. P1  Yao 2021 Table S3 pathway-recovery test            (src/yao.py)
  2. P2  knockout rankings, then scoring against literature (src/knockout.py, src/knockout_eval.py)
  3. P3  other endpoints on GSE30321                         (src/endpoints.py; ~20 min on 4 cores)
  4. P4  cross-modality alignment                            (src/alignment.py)
  5. P5  culture-format contrast                             (src/scale_format.py)
  6. metrics_v2.csv and figures
Fetches what it needs into raw/ (git-ignored).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))


def main(quick=False):
    import alignment
    import assemble_v2
    import endpoints
    import figures_v2
    import knockout
    import knockout_eval
    import scale_format
    import yao

    for name, fn in [("P1 Yao pathway recovery", yao.main), ("P2 knockout freeze", knockout.freeze),
                     ("P2 knockout evaluation", knockout_eval.main), ("P3 endpoints", lambda: endpoints.main(quick=quick)),
                     ("P4 alignment", alignment.main), ("P5 culture format", scale_format.main),
                     ("metrics_v2", assemble_v2.main), ("figures", figures_v2.make_all)]:
        t = time.time()
        print(f"=== {name}")
        fn()
        print(f"=== {name} done in {time.time() - t:.0f}s")


if __name__ == "__main__":
    main(quick="--quick" in sys.argv)

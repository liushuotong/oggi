"""Isolated worker for the copied OrthoFinder single-family HOG pipeline."""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cluster_engine.hog_tree import infer_hogs


if __name__ == '__main__':
    try:
        request = json.loads(Path(sys.argv[1]).read_text(encoding='utf-8'))
        result = infer_hogs(**request)
        print('HOG inference complete: %d HOGs at %s, %d unresolved targets' %
              (result['selected_hog_count'], result['selected_level'], len(result['unsupported'])), flush=True)
    except Exception as exc:
        print('HOG inference failed: %s: %s' % (type(exc).__name__, exc), file=sys.stderr)
        sys.exit(2)

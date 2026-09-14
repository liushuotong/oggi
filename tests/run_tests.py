"""Strict CI entry point: failed, skipped or undiscovered tests cannot pass."""
import sys
import unittest
from pathlib import Path


def exit_code(result):
    return 0 if result.testsRun > 0 and result.wasSuccessful() and not result.skipped else 1


def main():
    directory = Path(__file__).resolve().parent
    suite = unittest.defaultTestLoader.discover(str(directory), pattern='test_*.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    code = exit_code(result)
    if result.skipped:
        print('STRICT TEST FAILURE: skipped tests are not accepted. Install all test dependencies.', file=sys.stderr)
    if result.testsRun == 0:
        print('STRICT TEST FAILURE: no tests discovered.', file=sys.stderr)
    return code


if __name__ == '__main__':
    sys.exit(main())

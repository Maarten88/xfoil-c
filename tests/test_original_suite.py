"""The original xfoil-python test suite (xfoil/test.py of the reference checkout with docs/xfoil-python.patch;
skipped for an unpatched upstream checkout, whose test.py expects 4 return values from a()/cl()), unchanged
except for its two package-relative imports, run against xfoil_c: `from .xfoil import XFoil` ->
`from xfoil_c.xfoil import XFoil` and `from .model import Airfoil` -> `from xfoil_c.model import Airfoil`."""
import types
from ctypes import byref, c_bool
from pathlib import Path

import pytest

from conftest import require_xfoil_python

ORIGINAL = Path(require_xfoil_python('xfoil', 'test.py', patched=True))


def _load():
    src = ORIGINAL.read_text()
    for old, new in (('from .xfoil import XFoil', 'from xfoil_c.xfoil import XFoil'),
                     ('from .model import Airfoil', 'from xfoil_c.model import Airfoil')):
        assert src.count(old) == 1, old
        src = src.replace(old, new)
    mod = types.ModuleType('xfoil_c_original_test')
    exec(compile(src, str(ORIGINAL), 'exec'), mod.__dict__)
    return mod


TestXFoil = _load().TestXFoil          # collected by pytest as a unittest.TestCase


@pytest.fixture(autouse=True)
def _printing_off_afterwards(clib):
    """The original tests leave printing on (library default); restore the fixture contract."""
    yield
    clib.set_print(byref(c_bool(False)))

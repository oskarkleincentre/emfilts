from __future__ import absolute_import
from .version import __VERSION__ as __version__
from .spatial import *

here = __file__
basedir = os.path.split(here)[0]
example_data = os.path.join(basedir, 'example_data')

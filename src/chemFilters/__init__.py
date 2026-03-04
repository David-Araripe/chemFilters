# -*- coding: utf-8 -*-

from .version_helper import get_version

from .filters.rdkit_filters import RdkitFilters  # noqa F401

try:
    from .filters.bloom_filters import MolbloomFilters  # noqa F401
except ImportError:
    pass

try:
    from .filters.pep_filters import PeptideFilters  # noqa F401
except ImportError:
    pass

try:
    from .filters.silly_filters import SillyMolSpotterFilter  # noqa F401
except ImportError:
    pass

__version__ = get_version()

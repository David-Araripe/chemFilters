# ABOUTME: Sphinx configuration for chemFilters documentation.
# ABOUTME: Configures extensions, theme (Furo), and autodoc settings.

import chemFilters

project = "chemFilters"
copyright = "2024, David Araripe"
author = "David Araripe"
version = chemFilters.__version__
release = version

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

# MyST settings
myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
]

# Autodoc settings
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

# Napoleon settings (Google-style docstrings)
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True

# Intersphinx
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "rdkit": ("https://www.rdkit.org/docs/", None),
}

# Theme
html_theme = "furo"
html_title = " "
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "light_logo": "logo-light.svg",
    "dark_logo": "logo-dark.svg",
}

# Source settings
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

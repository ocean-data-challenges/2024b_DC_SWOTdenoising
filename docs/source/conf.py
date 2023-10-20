# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../src'))
for x in os.walk('../../src'):
  sys.path.insert(0, x[0])

# this is to tell reathedocs not to try to document numpy which is external.
autodoc_mock_imports = ['numpy','xarray','matplotlib','pyinterp','netCDF4','scipy','numba','pandas','datetime','cartopy','hvplot','warnings','packaging','dask','src']


# -- Project information -----------------------------------------------------

project = 'SWOT denoising data challenge'
copyright = '2023, Datlas'
author = 'Sammy Metref'

# The full version, including alpha/beta/rc tags
release = ''

master_doc = 'index'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'myst_parser',
    "nbsphinx",
    "sphinx_gallery.load_style",
]


napoleon_google_docstring = False

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}
 
myst_enable_extensions = ["dollarmath", "amsmath"]

myst_footnote_transition = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['gallery']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
#html_theme = 'classic'
html_theme = 'sphinx_rtd_theme'
#html_theme = 'alabaster'
#html_theme = 'nature'
#html_theme = 'pyramid'

html_logo = "figures/dc_2024_WOC-ESA_map2.jpg"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


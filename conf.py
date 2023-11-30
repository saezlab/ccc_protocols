# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information -----------------------------------------------------

project = 'LIANA x Tensor-cell2cell'
copyright = '2023, Hratch Baghdassarian, Daniel Dimitrov, Erick Armingol, et al'
author = 'Hratch Baghdassarian, Daniel Dimitrov, Erick Armingol, et al'

# The full version, including alpha/beta/rc tags
release = '0.0.3'


# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosectionlabel',
    'numpydoc',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []




# -- Options for HTML output -------------------------------------------------

master_doc = 'index'

html_theme = 'sphinx_rtd_theme'

html_static_path = ["_static"]
html_theme_options = dict(
    logo_only=True,
    display_version=True,
)
html_context = dict(
    display_github=False,
    github_user='hmbaghdassarian/',
    github_repo='ccc_protocols',
    github_version='main',
    conf_py_path='.',
)

html_show_sphinx = True
html_logo = '_static/logo.png'
html_favicon = '_static/logo.png'
html_css_files = ['custom.css']

# -- Options for EPUB output
epub_show_urls = 'footnote'

nbsphinx_execute = 'never'

nbsphinx_allow_errors = True

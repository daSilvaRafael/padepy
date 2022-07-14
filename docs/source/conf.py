# Configuration file for the Sphinx documentation builder.
import sys
import os
sys.path.insert(0, os.path.abspath('../'))
# -- Project information

project = 'padepy'
copyright = '2022, Rafael'
author = 'Rafael Pereira da Silva'

release = '0.1'
version = '0.1'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'prev_next_buttons_location': None
}
# -- Options for EPUB output
epub_show_urls = 'footnote'

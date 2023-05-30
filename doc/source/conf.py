# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

from importlib.metadata import version as metadata_version

project = 'ENDF Python Interface'
copyright = '2023, Paul Romano'
author = 'Paul Romano'

release = metadata_version('endf')
version = '.'.join(release.split('.')[:2])

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_autodoc_typehints',
    'sphinx_design',
]

templates_path = ['_templates']
exclude_patterns = []

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']

html_theme_options = {
    "github_url": "https://github.com/paulromano/endf-python",
    #"navbar_end": ["navbar-icon-links"],
    "show_toc_level": 3,
    "logo": {"text": project}
}

napoleon_use_rtype = False

# TODO: Using ivar results in better looking documentation, but it breaks
# cross-references to attributes which is annoying. Figure out a way to get it
# to work with cross-references

#napoleon_use_ivar = True

# -- Options for LaTeX output ------------------------------------------------
latex_domain_indices = False

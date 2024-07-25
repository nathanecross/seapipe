# Configuration file for the Sphinx documentation builder.

# -- Project information

project = 'seapipe'
copyright = '2024, N. Cross & A. Perrault'
author = 'Nathan Cross'

release = '0.1'
version = '0.1.0'

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
html_static_path = ['_static']
html_logo = 'logo.png'
html_theme_options = {
    'style_nav_header_background': '#030D7B',
    'logo_only': True,
    'logo_only': True,
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': True,
    'titles_only': False,
    'style_external_links': True,
    'display_version': True,
}

# -- Options for EPUB output
epub_show_urls = 'footnote'

# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# Project information
project = 'Physica'
copyright = '2021-2024, WeiBo He'
author = 'WeiBo He'
release = '0.0.1'

# General
extensions = ['myst_parser']
myst_enable_extensions = ['dollarmath']

'''
Mathjax may be blocked by the GFW and math equations cannot render correctly. Use CDNs provided by [1].

Reference:
[1] https://blog.csdn.net/qq_29654777/article/details/122903558
'''
mathjax_path = 'https://cdn.staticfile.net/mathjax/3.2.2/es5/tex-mml-chtml.js'

templates_path = ['_templates']
language = 'zh_CN'
exclude_patterns = []

# HTML Options
html_theme = 'conestack'
html_static_path = []

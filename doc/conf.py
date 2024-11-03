"""
Copyright 2024 Weibo He.

This file is part of Physica.

Permission is granted to copy, distribute and/or modify this document
under the terms of the GNU Free Documentation License, Version 1.3
or any later version published by the Free Software Foundation;
with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You should have received a copy of the GNU Free Documentation License
along with Physica.  If not, see <https://www.gnu.org/licenses/>.
"""
# Project information
project = 'Physica'
copyright = '2021-2024, Weibo He'
author = 'Weibo He'
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

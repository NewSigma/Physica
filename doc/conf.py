# This file is part of Physica.
#
# SPDX-FileCopyrightText: 2024 Weibo He
# SPDX-License-Identifier: GFDL-1.3-no-invariants-or-later
# Project information
project = 'Physica'
release = '2026'
author = 'Weibo He'
copyright = f'2021-{release}, {author}'

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
language = 'en'
exclude_patterns = []

# HTML Options
html_theme = 'conestack'
html_static_path = []

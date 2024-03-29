# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Astropy documentation build configuration file.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this file.
#
# All configuration values have a default. Some values are defined in
# the global Astropy configuration which is loaded here before anything else.
# See astropy.sphinx.conf for which values are set there.

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
# sys.path.insert(0, os.path.abspath('..'))
# IMPORTANT: the above commented section was generated by sphinx-quickstart, but
# is *NOT* appropriate for astropy or Astropy affiliated packages. It is left
# commented out with this explanation to make it clear why this should not be
# done. If the sys.path entry above is added, when the astropy.sphinx.conf
# import occurs, it will import the *source* version of astropy instead of the
# version installed (if invoked as "make html" or directly with sphinx), or the
# version in the build directory (if "python setup.py build_sphinx" is used).
# Thus, any C-extensions that are needed to build the documentation will *not*
# be accessible, and the documentation will not build correctly.

# import os
import sys
import datetime
# from importlib import import_module

try:
    from sphinx_astropy.conf.v1 import *  # noqa
except ImportError:
    print('ERROR: the documentation requires the sphinx-astropy '
          'package to be installed')
    sys.exit(1)

# Get configuration information from setup.cfg
from configparser import ConfigParser
conf = ConfigParser()
conf.read([os.path.join(os.path.dirname(__file__), '..', 'setup.cfg')])
setup_cfg = dict(conf.items('metadata'))

# -- General configuration ----------------------------------------------------

# Add package root directory (one level above this file) to sys.path:
sys.path.insert(0, os.path.abspath('..'))

graphviz_dot = 'C:/Programs/graphviz/bin/dot.exe'

# By default, highlight as Python 3.
highlight_language = 'python3'

# add any custom intersphinx for this package
intersphinx_mapping['astropy'] = ('https://docs.astropy.org/en/latest/', None)  # noqa: F405
intersphinx_mapping['astropy-dev'] = ('https://docs.astropy.org/en/latest/', None)  # noqa: F405
# intersphinx_mapping['pyerfa'] = ('https://pyerfa.readthedocs.io/en/stable/', None)  # noqa: F405
intersphinx_mapping['pytest'] = ('https://docs.pytest.org/en/stable/', None)  # noqa: F405
# intersphinx_mapping['ipython'] = ('https://ipython.readthedocs.io/en/stable/', None)  # noqa: F405
intersphinx_mapping['pandas'] = ('https://pandas.pydata.org/pandas-docs/stable/', None)  # noqa: F405, E501
intersphinx_mapping['sphinx_automodapi'] = ('https://sphinx-automodapi.readthedocs.io/en/stable/', None)  # noqa: F405, E501
intersphinx_mapping['packagetemplate'] = ('https://docs.astropy.org/projects/package-template/en/latest/', None)  # noqa: F405, E501
# intersphinx_mapping['h5py'] = ('https://docs.h5py.org/en/stable/', None)  # noqa: F405
intersphinx_mapping['statsmodels'] = ('https://www.statsmodels.org/stable/', None)  # noqa: F405
intersphinx_mapping['photutils'] = ('https://photutils.readthedocs.io/en/stable/', None)  # noqa: F405
# intersphinx_mapping['skyfield'] = ('https://rhodesmill.org/skyfield/documentation/', None)  # noqa: F405

# If your documentation needs a minimal Sphinx version, state it here.
# needs_sphinx = '1.2'

# To perform a Sphinx version check that needs to be more specific than
# major.minor, call `check_sphinx_version("X.Y.Z")` here.
# check_sphinx_version("1.2.1")

# Whether to create cross-references for the parameter types in the
# Parameters, Other Parameters, Returns and Yields sections of the docstring.
numpydoc_xref_param_type = True

# Additional, package-specific cross-references to fully qualified paths
# (or correct ReST references) for the
# aliases/shortcuts used when specifying the types of parameters.
# (Numpy provides some defaults
# https://github.com/numpy/numpydoc/blob/b352cd7635f2ea7748722f410a31f937d92545cc/numpydoc/xref.py#L62-L94
# and a base set comes from sphinx-astropy.)
numpydoc_xref_aliases.update({
    # python & adjacent
    "Any": "`~typing.Any`",
    "file-like": ":term:`python:file-like object`",
    "file": ":term:`python:file object`",
    "path-like": ":term:`python:path-like object`",
    "module": ":term:`python:module`",
    "buffer-like": ":term:buffer-like",
    "hashable": ":term:`python:hashable`",
    # for matplotlib
    "color": ":term:`color`",
    # for numpy and python primitives
    # "ints": ":class:`python:int` values",
    # "floats": ":class:`python:float` values",
    # "tuples": ":class:`python:tuple` values",
    # for astropy
    "number": ":term:`number`",
    "writable": ":term:`writable file-like object`",
    "readable": ":term:`readable file-like object`",
    "BaseHDU": ":doc:`HDU </io/fits/api/hdus>`"
})

# Add from sphinx-astropy 1) glossary aliases 2) physical types.
numpydoc_xref_aliases.update(numpydoc_xref_astropy_aliases)


# -- Project information ------------------------------------------------------

# This does not *have* to match the package name, but typically does
project = setup_cfg['name']
author = setup_cfg['author']
copyright_date = datetime.datetime.utcnow().year
if copyright_date > 2022:
    copyright = f'2022–' + str(copyright_date) + '  ' + author
else:
    copyright = f'2022  ' + author
today_fmt = '%Y-%m-%d'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

#
# import_module(setup_cfg['name'])
# package = sys.modules[setup_cfg['name']]
#
# # The short X.Y version.
# version = package.__version__.split('-', 1)[0]
# # The full version, including alpha/beta/rc tags.
# release = package.__version__

version = '0.1' + 'beta'
release = version
root_doc = 'astropack/index'

# -- Options for HTML output --------------------------------------------------

# A NOTE ON HTML THEMES
# The global astropy configuration uses a custom theme, 'bootstrap-astropy',
# which is installed along with astropy. A different theme can be used or
# the options for this theme can be modified by overriding some of the
# variables set in the global configuration. The variables set in the
# global configuration are listed below, commented out.

# Add any paths that contain templates here, relative to this directory.
# if 'templates_path' not in locals():  # in case parent conf.py defines it
#     templates_path = []
# templates_path.append('_templates')

# Add any paths that contain custom themes here, relative to this directory.
# To use a different custom theme, add the directory containing the theme.
# html_theme_path = []

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes. To override the custom theme, set this to the
# name of a builtin theme or the name of a custom theme in html_theme_path.
html_theme = "bootstrap-astropy"
html_theme_options = {
    'logotext1': 'astropack',   # white, semi-bold
    'logotext2': '',  # orange, light
    'logotext3': ':docs'   # white, light
}
# html_static_path = ['_static']
html_style = 'astropack.css'

# html_short_title = 'astropack'
html_search_language = 'en'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns.append('_templates')
exclude_patterns.append('changes')  # noqa: F405
exclude_patterns.append('_pkgtemplate.rst')  # noqa: F405
# .inc.rst means *include* files, don't have sphinx process them
exclude_patterns.append('**/*.inc.rst')    # noqa: F405, E501

# This is added to the end of RST files - a good place to put substitutions to
# be used globally.
with open("common_links.txt", "r") as cl:
    rst_epilog += cl.read()

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# Add any paths that contain templates here, relative to this directory.
# templates_path = ['_templates']

html_permalinks = False  # EVD 2022-03-21, was True.
source_suffix = '.rst'

# Custom sidebar templates, maps document names to template names.
# html_sidebars = {}

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = ''

# The name of an image file (in conf.py's own directory) to use as favicon of the
# docs. This file should be a Windows icon file (.ico) of size 16x16
# or preferably 32x32 pixels.
html_favicon = 'backpack-star-favicon3.ico'

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = '%Y-%m-%d'

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = '{0} v{1}'.format(project, release)

# Output file base name for HTML help builder.
htmlhelp_basename = project + 'doc'

# Prefixes that are ignored for sorting the Python module index
modindex_common_prefix = ["astropack."]

# Required by sphinx-automodapi extension:
numpydoc_show_class_members = False

autosummary_generate = True
automodsumm_inherited_members = True

# autodoc_member_order = 'alphabetical'  # doesn't work with automodapi.

# autodoc_default_options = {'member-order': 'alphabetical'}  # not with automodapi.

# print()
# for e in extensions:
#     print(str(e))

# -- Options for LaTeX output -------------------------------------------------

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [('index', project + '.tex', project + u' Documentation',
                    author, 'manual')]


# -- Options for manual page output -------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', project.lower(), project + u' Documentation',
              [author], 1)]


# -- Options for the edit_on_github extension ---------------------------------

if setup_cfg.get('edit_on_github').lower() == 'true':

    extensions += ['sphinx_astropy.ext.edit_on_github']

    edit_on_github_project = setup_cfg['github_project']
    edit_on_github_branch = "main"

    edit_on_github_source_root = ""
    edit_on_github_doc_root = "docs"

# -- Resolving issue number to links in changelog -----------------------------
github_issues_url = 'https://github.com/{0}/issues/'.format(setup_cfg['github_project'])


# -- Options for linkcheck output -------------------------------------------
linkcheck_retry = 5
linkcheck_ignore = [
    r'https://github\.com/edose/astropack/(?:issues|pull)/\d+',
]
linkcheck_timeout = 180
linkcheck_anchors = False

# -- Turn on nitpicky mode for sphinx (to warn about references not found) ----
#
# nitpicky = True
# nitpick_ignore = []
#
# Some warnings are impossible to suppress, and you can list specific references
# that should be ignored in a nitpick-exceptions file which should be inside
# the docs/ directory. The format of the file should be:
#
# <type> <class>
#
# for example:
#
# py:class astropy.io.votable.tree.Element
# py:class astropy.io.votable.tree.SimpleElement
# py:class astropy.io.votable.tree.SimpleElementWithContent
#
# Uncomment the following lines to enable the exceptions:
#
# for line in open('nitpick-exceptions'):
#     if line.strip() == "" or line.startswith("#"):
#         continue
#     dtype, target = line.split(None, 1)
#     target = target.strip()
#     nitpick_ignore.append((dtype, six.u(target)))

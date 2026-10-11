"""Sphinx configuration for the heatwave3 Python documentation."""

from importlib.metadata import version as _version

project = "heatwave3"
author = "Robert W. Schlegel and Albertus J. Smit"
copyright = "2023-2026, heatwave3 authors"
release = _version("heatwave3")
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "numpydoc",
    "myst_nb",
]
source_suffix = {".md": "myst-nb", ".rst": "restructuredtext", ".ipynb": "myst-nb"}
exclude_patterns = ["_build", "generate_examples.py"]

# numpydoc: the docstrings are NumPy style; no automatic member tables
numpydoc_show_class_members = False
numpydoc_xref_param_type = False

# myst-nb: render the notebook from its saved outputs (it needs R to execute)
nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

html_theme = "pydata_sphinx_theme"
html_title = f"heatwave3 {release}"
html_logo = "../../logo.png"
html_theme_options = {
    "github_url": "https://github.com/robwschlegel/heatwave3",
    "icon_links": [
        {"name": "PyPI", "url": "https://pypi.org/project/heatwave3/",
         "icon": "fa-brands fa-python"},
        {"name": "R package", "url": "https://robwschlegel.github.io/heatwave3/",
         "icon": "fa-brands fa-r-project"},
    ],
    "navigation_with_keys": False,
}

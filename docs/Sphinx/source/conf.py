# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
import os
import re
import subprocess

pybuc_path = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../python")
)
sys.path.insert(0, pybuc_path)

with open(os.path.join(pybuc_path, "version.hpp")) as f:
    for line in f:
        vline = re.search(
            '#define BOBKIT_VERSION "(?P<value>[0-9dev.-]+)"', line
        )  # noqa 501
        if vline is not None:
            VERSION = re.sub(
                '.+"([0-9]+.[0-9]+.[0-9]+)(-dev)?"', "\\1\\2", vline.group()
            )
            break

project = "BOBKit"
copyright = "2024, Soon Wen Hoh"
author = "Soon Wen Hoh"
release = VERSION

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "breathe",
    "sphinx_rtd_theme",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- Breathe configuration ---------------------------------------------------
# -- Generating Doxygen xml files --------------------------------------------
read_the_docs_build = os.environ.get("READTHEDOCS", "") == "True"
if read_the_docs_build:
    subprocess.call(
        "cd ../../Doxygen ; doxygen Doxyfile-bobkit.cfg", shell=True
    )  # noqa501

breathe_projects = {"buccaneer": "../../Doxygen/generated_docs/xml"}
breathe_default_project = "BOBKit"

primary_domain = "cpp"
highlight_language = "cpp"

# autosummary
autosummary_generate = True

# autodoc patch sys.path
# sys.path.insert(0, os.path.abspath("../../../python"))
autodoc_default_options = {
    "members": None,
    "imported-members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import sys
import os
import re
from datetime import datetime

# import subprocess

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
author = "Soon Wen Hoh"
copyright = f"{datetime.now().year}, {author}"
release = VERSION

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # "breathe",
    "sphinx_rtd_theme",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.viewcode",
    # "myst_parser",
    # "sphinx_immaterial.apidoc.python.apigen",
    # "sphinx_autodoc_typehints",
    # "autoapi.extension",
]

# python_apigen_modules = {
#    "bobkit": "api/bobkit.",
# }
#
# python_apigen_default_groups = [
#    (".*", "Public members"),
#    ("class:.*", "Classes"),
#    (r".*\.__(init|new)__", "Constructors"),
#    (r".*\.__(str|repr)__", "String representation"),
# ]
#
# python_apigen_rst_prolog = """
# .. default-role:: py:obj
#
# .. default-literal-role:: python
#
# .. highlight:: python
# """
templates_path = ["_templates"]
exclude_patterns = []
# html_theme_options = {
#    "font": False,
# }

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
# html_theme = "sphinx_immaterial"
# html_static_path = ["_static"]

# -- Breathe configuration ---------------------------------------------------
# -- Generating Doxygen xml files --------------------------------------------
# read_the_docs_build = os.environ.get("READTHEDOCS", "") == "True"
# dox_path = os.path.abspath(
#    os.path.join(os.path.dirname(os.path.abspath(__file__)), "./../../Doxygen")
# )
# sys.path.insert(0, dox_path)
# if read_the_docs_build:
#    subprocess.call("cd ../../Doxygen ; doxygen", shell=True)  # noqa501
#    # dox_path = os.path.abspath(
#    #    os.path.join(os.path.dirname(os.path.abspath(__file__)), "./../../Doxygen")
#    # )
#    # sys.path.insert(0, dox_path)
#
# doxygen_xml_path = os.path.join(dox_path, "generated_docs/xml")
# breathe_projects = {"buccaneer": doxygen_xml_path}
# breathe_default_project = "BOBKit"

# primary_domain = "cpp"
# highlight_language = "cpp"

# autosummary
autosummary_generate = True

# autodoc patch sys.path
# sys.path.insert(0, os.path.abspath("../../../python"))
# sys.path.insert(0, os.path.abspath("../../../bobkit-stubs"))
# stubspath = os.path.join(pybuc_path, "bobkit")
# sys.path.insert(0, stubspath)
# autodoc_mock_imports = ["bobkit"]
# autodoc_mock_impor
autodoc_default_options = {
    "members": None,
    "imported-members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# autoapi directives
# pystubs_path = os.path.abspath(
#    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../../bobkit-stubs")
# )
# sys.path.insert(0, pystubs_path)
# autoapi_dirs = [pystubs_path]

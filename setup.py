
# from glob import glob
# from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# ext_modules = [
#     #Pybind11Extension("buildkit",["python/buildkit.cpp"], define_macros=[('VERSION', __version__)], ), # noqa: E501
#     Pybind11Extension("buildkit", sorted(glob("python/*.cpp")), define_macros=[('VERSION', __version__)], ), # noqa: E501
# ]

if __name__ == "__main__":
    setup()
# setup(
#    name="buildkit",
#    version=__version__,
#    author="Soon Wen Hoh",
#    author_email="soonwen.hoh@york.ac.uk",
#    description="Library for model building.",
#    long_description="",
#    ext_modules=ext_modules,
#    extras_require={"test": "pytest"},
#    cmdclass={"build_ext": build_ext},
#    python_requires=">=3.7")

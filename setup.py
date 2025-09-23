# from glob import glob
# from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension, find_packages
import glob


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


ext_modules = [
    CMakeExtension("bobkit_ext.clipper"),
    CMakeExtension("bobkit_ext.buccaneer"),
    CMakeExtension("bobkit_ext.protein_db"),
    CMakeExtension("bobkit_ext.util"),
    #    "_bobkit",
    #    glob.glob("include/_bobkit/*.cpp") + glob.glob("python/*.cpp"),
    #    include_dirs=["include/buccaneer/", "python/"],
    #    language="c++",
    #    extra_compile_args=["-std=c++11"],
    #    define_macros=[("DOCTEST_CONFIG_DIABLE", None)],
    # )
]

# ext_modules = [
#     #Pybind11Extension("buildkit",["python/buildkit.cpp"], define_macros=[('VERSION', __version__)], ), # noqa: E501
#     Pybind11Extension("buildkit", sorted(glob("python/*.cpp")), define_macros=[('VERSION', __version__)], ), # noqa: E501
# ]
setup(name="bobkit", packages=find_packages(), ext_modules=ext_modules)

# if __name__ == "__main__":
#    setup()
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

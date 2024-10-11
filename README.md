# BOBKit - **B**ioMacromolecular M**o**del **B**uilding **Kit**
## **Buccaneer**
Buccaneer is a C++ software for statistical protein chain tracing. It identifies connected alpha-carbon positions using a likelihood-based density target.

## **BOBKit**
BOBKit is the python bindings for Buccaneer. We have used pybind11 for this project for now.
Individual steps from Buccaneer can be run as python methods with the appropriate arguments.
Currently, it requires GEMMI to be built using the same compiler for the python bindings to work properly. This will 

## **Installing**
1. Clone the repository to your local machine.

   `git clone https://github.com/scotthoh/BOBKit.git`

2. Go into cloned directory: BOBKit

   `cd BOBKit`

3. Create a virtual environment with the venv python module. Currently tested on python3.8

   `python3 -m venv bobkit_env`

4. Activate virtual environment

   `source bobkit_env/bin/activate`

5. Installing BOBKit modules using pip.* This will first download, compile and install required dependencies (FFTW3 v3.3.4, GEMMI v0.6.4, and a branch of Clipper library from CCP4 repository) before installing the BOBKit python bindings.

   `python3 -m pip install .`

   *It is recommended to install BOBKit using pip to achieve best desired results. Building it manually with CMake might also work.

   It is also possible to skip the installation of FFTW3 and use the version of FFTW3 (>=3.3.4) on your system. Requirements for this is pip>=24.2 to use the option --config-settings (seems there is a bug for pip 23 [mentioned in this issue](https://github.com/pypa/pip/issues/11325#issuecomment-1474703878))
   
   `python3 -m pip install . --config-settings=cmake.args=-DSKIP_INSTALL_FFTW3=ON`

More documentation are being added, work in progres...
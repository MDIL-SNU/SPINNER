# SPINNER
SPINNER(Structure Prediction of Inorganic crystals using Neural Network potentials with Evolutionary and Random searches)

If you use SPINNER, please cite this article: 

Here do we describe minimal instruction to run the example of SPINNER
If you want more information such as tuning parameters, please visit our online manual(https://spinner-csp.readthedocs.io)

## Installation

### Download SPINNER
```
git clone https://github.com/MDIL-SNU/SPINNER.git
```

### Python requirements

SPINNER supports Python :code:`3.6` or higher version. SPINNER utilizes Python modules in the following:

  - mpi4py (https://bitbucket.org/mpi4py/mpi4py/)
  - PyYAML (https://pypi.org/project/PyYAML/)
  - numpy (https://numpy.org/)
  - pybind11 (https://github.com/pybind/pybind11/)

You can download these modules by

```

  pip3 install mpi4py PyYAML numpy pybind11
```

or directly downloading the source codes from the webpages.

### C++ compilers

SPINNER supports CMAKE 3.20 or higher version and gcc version 7.3.1 or higer version.

### Python binding of randSpg and LAMMPS (important)
SPINNER utilizes python binding of randSpg and LAMMPS code. However, those codes are slightly modified to be incorporated with SPINNER, so you need to comiple with source codes provided with SPINNER, not with the original ones.

### Install randSpg
To install randSpg, do follows

```

  cd /SPINNER-directory/randSpg-vspinner/
  mkdir build
  cd build
  cmake ..
  make -j3
```

to bind randSpg with python, do follows

```

  cd /SPINNER-directory/randSpg-vspinner/python
  mkdir build
  cd build
  cmake ..
  make –j3
  cp pyrandspg.cpython* /directory-where-your-python-is/lib/python3/site-packages/
```

### Install LAMMPS
This is the most tricky part of the installation steps. You can install python-lammps by follows but it may not work depending on your system. Please look into LAMMPS forum (https://www.lammps.org/forum.html) or LAMMPS manual Ch. 2 (https://docs.lammps.org/Python_run.html) for detailed discussion.

```

  cd /SPINNER-directory/lammps-vspinner
  cd src
  make yes-python
  make XXX mode=shlib
  make install-python
```

Here, XXX is the name of the make file (Makefile.XXX in src/MAKE directory). Note that LAMMPS have to be installed with serial version (not mpi version). The optimization is recommended for your system but default serial option is also possible.

## Usage
To use SPINNER, 1 file (XXX.yaml) and 2 directories (input directory and src) are required.

### input file (XXX.yaml; XXX is named by the user)
Parameter list to control SPINNER code is listed in XXX.yaml. 
The simplest form of input.yaml is described below:
```YAML
# XXX.yaml

  input_dir:        input_directory_name
  output_dir:       output_directory_name
  initial_volume:   326.0
  structure:
      generation: 400
  material:
      Li:  4
      Hf:  2
      N:   4
```

### input directory
In input directory (input_dir: in input file), LAMMPS potential file should be located. (Potential file have to be named potential.) Potential should be SIMPLE-NN format (https://github.com/MDIL-SNU/SIMPLE-NN).

### src directory
src directory should be in the running directory. You can copy and paste the src directory to each running directory or you can run multiple calculations in one running folder.


## Running code

```

  cd src
  mpirun -np core_number python3 main.py XXX.yaml
```

## Example
See examples directory.

```

  sh run.sh
```

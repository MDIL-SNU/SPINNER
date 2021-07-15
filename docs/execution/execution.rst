=======================
2. Execution of SPINNER
=======================

2.1 Required input files
========================
Running SPINNER requires one :code:`input file` (yaml format), one :code:`input directory` (potential file should be in it), and :code:`src` directory. All of these files and directories should be in the running directory.

2.1.1 Input file
----------------
SPINNNER uses YAML style configuration file. All setting parameters used in SPINNER can be controlled in “XXX.yaml” while XXX is the specific name that user defined. These are the mandatory setting:

::

  input_dir:        input_directory_name
  output_dir:       output_directory_name
  initial_volume:   326.0
  structure:
      generation: 400
  material:
      Li:  4
      Hf:  2
      N:   4

See other tags for :doc:`/tags/tags` section. 

2.1.2 Input directory
---------------------
In input directory (:code:`input_dir` in input file), LAMMPS potential file :code:`potential` should be located. (Potential file have to be named :code:`potential`.) Potential should be SIMPLE-NN format (https://github.com/MDIL-SNU/SIMPLE-NN).

2.1.3 src directory
-------------------
:code:`src` directory should be in the running directory. You can copy and paste the src directory to each running directory or you can run multiple calculations in one running folder.

2.2 Running SPINNER
===================

You can execute SPINNER by follows

::

  cd src
  mpirun -np core_number python3 main.py XXX.yaml

where XXX.yaml is the input file.

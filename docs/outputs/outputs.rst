===============
4. Output files
===============

4.1 Energy and volume data
==========================

 - :code:`totalGen`: information such as energy, volume, and generation methods of total population.
 - :code:`totalbest`: information such as energy, volume, and generation methods of the best structure in each generation.
 - :code:`BestResults`: information of best results that are survived in the current generation (updated at each generation)
 - :code:`best_history`: history of BestResults file

 Energy 10000.0 eV in these files means that the structure violates the geometrical constraints or have the same structure as the other structure in the pool.


4.2 Structure files
===================

 - :code:`POSCARs`: gathered structures before relaxation
 - :code:`CONTCARs`: relaxed structures


4.3 Time log
============
 - :code:`timelog`: log file containing information of time
 - :code:`specific_time`: log file containing time information of each core used in the calculation


4.4 Input information
=====================

 - :code:`input.yaml`: total input used in the calculation
 - :code:`potential`: potential lfile used in the calculation


4.5 Etc output files
====================

 - :code:`ave_atomic_e`: average atomic enegy of each element (used in the crossover algorithm)
 - :code:`distance_info`: similarity distance information: 
   format: (generation) (population1) (population2) (distance)
 - :code:`random_structure_log`: information of generation time and spacegroup (space group of initial structure not the final structure) of random structure.


4.6 Directories and files for calculation (does not contain information)
========================================================================

 - :code:`popXX` (XX=integer)(directories): lammps calculations run in this folder
 - :code:`random_structure` (directory): random structure generations run in this folder
 - :code:`potential10`: potential file to generate rdf (generated automatically)



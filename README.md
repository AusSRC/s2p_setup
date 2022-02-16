# s2p_setup

The purpose of the `s2p_setup.py` script is to assist with setting up
SoFiA-X runs. SoFiA-X is the parallel framework around the SoFiA 2
extragalactic HI source finding pipeline.


## Prerequisites

* SoFiA 2: https://github.com/SoFiA-Admin/SoFiA-2
* SoFiA-X: https://github.com/AusSRC/SoFiAX


## Usage     

    s2p_setup.py <s2p_template> <data_cube> <template_par_file> <unique_name> <output_directory>


## Arguments 

* `<s2p_template>`       Template `s2p_setup.ini` file.

* `<data_cube>`          Input FITS data cube on which SoFiA 2 is to be
                         run. Note that only the header of that file will
                         be extracted and no significant amount of memory
                         will be needed.

* `<template_par_file>`  SoFiA 2 template parameter file from which the
                         settings to be used for all regions will be
                         extracted.

* `<unique_name>`        Unique name to be used when creating the results
                         tables in the database. This must be unique to
                         ensure that any existing tables from previous
                         runs are not overwritten. It will also be used
                         to name SoFiA's output products and catalogues.

* `<output_directory>`   Output products directory for SoFiA 2 output products.

## Purpose   

This script can be used to automatically partition a data cube into
conveniently sized regions with a certain amount of overlap that can then
be fed into the parallel SoFiA-X framework. A few basic settings can be
adjusted in the `s2p_setup.ini` file. The script will also create the
necessary SoFiA-X setup files, and SoFiA-X can be launched by simply
calling the auto-generated `run_sofiax.sh` script. It is strongly
recommended to visually check and, if necessary, manually correct all
generated setup files before executing the `run_sofiax.sh` script to
launch SoFiA-X.

Note that this script does not actually cut up the data cube, but it
generates the required number of SoFiA parameter files for reading in and
processing sub-regions of the cube.


## Copyright and licence

Copyright (C) 2020 Tobias Westmeier

s2p_setup is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see http://www.gnu.org/licenses/.

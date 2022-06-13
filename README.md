# s2p_setup

The purpose of the `s2p_setup.py` script is to assist with setting up
SoFiA-X runs. SoFiA-X is the parallel framework around the SoFiA 2
extragalactic HI source finding pipeline.


## Prerequisites

* SoFiA 2: https://github.com/SoFiA-Admin/SoFiA-2
* SoFiA-X: https://github.com/AusSRC/SoFiAX


## Usage

```
s2p_setup.py \
    --config <s2p_template> \
    --image_cube <image_cube> \
    --run_name <run_name> \
    --sofia_template <template_parameter_file> \
    --output_dir <output_directory> \
    --products_dir <products_directory>
```


## Arguments

Required arguments:

| Argument | Description | 
| --- | --- |
| `s2p_template` | Template `s2p_setup.ini` file. | 
| `image_cube` | Input FITS data cube on which SoFiA 2 is to be run. Note that only the header of that file will be extracted and no significant amount of memory will be needed. | 
| `run_name` | Name of the pipeline run. | 
| `template_par_file` | SoFiA 2 template parameter file from which the settings to be used for all regions will be extracted. |
| `output_directory` | Output products directory for SoFiA 2 output products. |
| `products_directory` | Sub-cube products sub-directory for SoFiA 2 output products. |
| `region*` | RA and declination values (RA min, RA max, Dec min, Dec max) for constraining region of image cube to extract. If not provided this will default to the pixel region defined in the configuration file. |

*Argument is optional

## Purpose   

This script can be used to automatically partition a data cube into
conveniently sized regions with a certain amount of overlap that can then
be fed into the parallel SoFiA-X framework. A few basic settings can be
adjusted in the `s2p_setup.ini` file.

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

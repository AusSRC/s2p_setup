#!/usr/bin/env python

"""
## ____________________________________________________________________ ##
##                                                                      ##
## s2p_setup.py - Setup script for the SoFiAX parallel framework        ##
## Copyright (C) 2020 Tobias Westmeier                                  ##
## ____________________________________________________________________ ##
##                                                                      ##
## Address:  Tobias Westmeier                                           ##
##           ICRAR M468                                                 ##
##           The University of Western Australia                        ##
##           35 Stirling Highway                                        ##
##           Crawley WA 6009                                            ##
##           Australia                                                  ##
##                                                                      ##
## E-mail:   tobias.westmeier [at] uwa.edu.au                           ##
## ____________________________________________________________________ ##
##                                                                      ##
## This program is free software: you can redistribute it and/or modify ##
## it under the terms of the GNU General Public License as published by ##
## the Free Software Foundation, either version 3 of the License, or    ##
## (at your option) any later version.                                  ##
##                                                                      ##
## This program is distributed in the hope that it will be useful,      ##
## but WITHOUT ANY WARRANTY without even the implied warranty of       ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         ##
## GNU General Public License for more details.                         ##
##                                                                      ##
## You should have received a copy of the GNU General Public License    ##
## along with this program. If not, see http://www.gnu.org/licenses/.   ##
## ____________________________________________________________________ ##
##                                                                      ##
"""


import sys
import math
import configparser
from astropy.io import fits


# Function to search for substring in list
def substr_search(input_list, substr):
    for i, s in enumerate(input_list):
        if substr in s:
            return i
    return -1


# Check command-line arguments
argv = sys.argv
argc = len(argv)

if (argc != 7):
    sys.stderr.write("\n\033[1736m Usage     \033[0m\n\n")
    sys.stderr.write("   \033[1ms2p_setup.py\033[0m \033[3m<setup_config> <data_cube> <sofia_par_file> <run_name> <node_size> <output_directory>\033[0m\n\n")  # noqa
    sys.stderr.write("\033[1736m Arguments \033[0m\n\n")
    sys.stderr.write("   \033[3m<setup_config>\033[0m       Template s2p_setup.ini file with default parameter values.\n\n")
    sys.stderr.write("   \033[3m<data_cube>\033[0m          Input FITS data cube on which SoFiA 2 is to be run.\n")
    sys.stderr.write("                        Note that only the header of that file will be ex-\n")
    sys.stderr.write("                        tracted and no significant amount of memory will be\n")
    sys.stderr.write("                        needed.\n\n")
    sys.stderr.write("   \033[3m<sofia_par_file>\033[0m     SoFiA 2 template parameter file from which the set-\n")
    sys.stderr.write("                        tings to be used for all regions will be extracted.\n\n")
    sys.stderr.write("   \033[3m<run_name>\033[0m           Unique name to be used when creating the results ta-\n")
    sys.stderr.write("                        bles in the database. This must be unique to ensure\n")
    sys.stderr.write("                        that any existing tables from previous runs are not\n")
    sys.stderr.write("                        overwritten. It will also be used to name SoFiA's\n")
    sys.stderr.write("                        output products and catalogues.\n\n")
    sys.stderr.write("   \033[3m<node_size>\033[0m          Size of the computing worker nodes available to the user.\n")
    sys.stderr.write("                        Optimal sub-cube split will be determined from this value.\n\n")
    sys.stderr.write("   \033[3m<output_directory>\033[0m   Directory where output parameter file will be written.\n\n")
    sys.stderr.write("\033[1736m Purpose   \033[0m\n\n")
    sys.stderr.write("   This script can be used to automatically partition a data cube into\n")
    sys.stderr.write("   conveniently sized regions with a certain amount of overlap that can\n")
    sys.stderr.write("   then be fed into the parallel \033[1mSoFiAX\033[0m framework. A few basic settings\n")
    sys.stderr.write("   can be adjusted in the \033[3ms2p_setup.ini\033[0m file. The script will also create\n")
    sys.stderr.write("   the necessary SoFiAX setup files, and SoFiAX can be launched by simply\n")
    sys.stderr.write("   calling the auto-generated \033[3mrun_sofiax.sh\033[0m script.\n\n")
    sys.stderr.write("   Note that this script does not actually cut up the data cube, but it\n")
    sys.stderr.write("   generates the required number of SoFiA parameter files for reading in\n")
    sys.stderr.write("   and processing sub-regions of the cube.\n\n")
    sys.exit(1)

s2p_setup_config = str(argv[1])     # s2p_setup template values
input_filename = str(argv[2])       # FITS data cube
template_file = str(argv[3])        # SoFiA 2 template parameter file to be used for all subregions
db_run_name = str(argv[4])          # Unique name of run (used for database table)
node_size = int(argv[5])            # Size of the computing nodes
output_dir = str(argv[6])           # Output directory for parameter files


# Read settings from configuration file
config = configparser.ConfigParser()
success = config.read(s2p_setup_config)

if(len(success) == 0):
    sys.stderr.write("Error: Failed to read config file: s2p_setup.ini\n")
    sys.stderr.write("Please ensure that a copy of that file is in the current directory.\n")
    sys.exit(1)

x_min = int(config["boundary"]["x_min"])
x_max = int(config["boundary"]["x_max"])
y_min = int(config["boundary"]["y_min"])
y_max = int(config["boundary"]["y_max"])
z_min = int(config["boundary"]["z_min"])
z_max = int(config["boundary"]["z_max"])
overlap_spat = int(config["region"]["overlap_spat"])
overlap_spec = int(config["region"]["overlap_spec"])
n_cpu_cores = int(config["pipeline"]["cpu_cores"])
tolerance_pos = int(config["pipeline"]["tolerance_pos"])
tolerance_flux = int(config["pipeline"]["tolerance_flux"])
tolerance_spat = (config["pipeline"]["tolerance_spat"]).split(",")
tolerance_spec = (config["pipeline"]["tolerance_spec"]).split(",")
tolerance_spat = [int(x) for x in tolerance_spat]
tolerance_spec = [int(x) for x in tolerance_spec]


# Calculate appropriate region size from node_size and estimated SoFiA overhead
region_size = int(node_size / 2.3)


# Try to open FITS file
try:
    hdu = fits.open(input_filename)
except Exception:
    sys.stderr.write("Error: Failed to open FITS file: {}\n".format(input_filename))
    sys.exit(1)


# Extract header information
header = hdu[0].header

bitpix = int(header["BITPIX"])
word_size = int(abs(bitpix) / 8)
if (bitpix > 0):
    wordsize = 4  # Assume 32-bit for integer arrays
    sys.stderr.write("Warning: Data cube is of integer type assuming 32 bit per pixel.\n")

naxis = int(header["NAXIS"])
if (naxis < 3 or naxis > 4):
    sys.stderr.write("Error: Data cube is not three-dimensional.\n")
    hdu.close()
    sys.exit(1)

nx = int(header["NAXIS1"])
ny = int(header["NAXIS2"])
nz = int(header["NAXIS3"])

if (nz == 1 and naxis == 4):
    nz = int(header["NAXIS4"])
    sys.stderr.write("Warning: Swapping 3rd and 4th axis of 4D cube.\n")


# Close FITS file again
hdu.close()


# Update region boundaries if necessary
if (x_min < 0):
    x_min = 0
if (y_min < 0):
    y_min = 0
if (z_min < 0):
    z_min = 0
if (x_max >= nx):
    x_max = nx - 1
if (y_max >= ny):
    y_max = ny - 1
if (z_max >= nz):
    z_max = nz - 1

sys.stdout.write("\nInput:\n")
sys.stdout.write("  Input range:     {0}-{1}, {2}-{3}, {4}-{5}\n".format(x_min, x_max, y_min, y_max, z_min, z_max))


# Determine nominal size of subregions in each dimension
size = int(math.floor((float(1024 * 1024 * 1024 * region_size) / float(word_size)) ** (1.0 / 3.0)))


# Determine practical region size
n_reg_x = int(math.ceil(float(x_max - x_min + 1) / float(size)))
n_reg_y = int(math.ceil(float(y_max - y_min + 1) / float(size)))
n_reg_z = int(math.ceil(float(z_max - z_min + 1) / float(size)))

if (n_reg_x * n_reg_y * n_reg_z > 999):
    sys.stderr.write(
        """
        Error: The setup would create {0:d} parallel runs!\n
        Please reconsider your settings to reduce that number.\n
        """.format(n_reg_x * n_reg_y * n_reg_z))
    sys.exit(1)

size_x = int(math.floor(float(x_max - x_min + 1 + (n_reg_x - 1) * overlap_spat) / float(n_reg_x)))
size_y = int(math.floor(float(y_max - y_min + 1 + (n_reg_y - 1) * overlap_spat) / float(n_reg_y)))
size_z = int(math.floor(float(z_max - z_min + 1 + (n_reg_z - 1) * overlap_spec) / float(n_reg_z)))

ram_per_region = size_x * size_y * size_z * word_size / (1024.0 * 1024.0 * 1024.0)  # GB
ram_per_node = int(math.ceil(2.3 * ram_per_region))  # GB

sys.stdout.write("\nOutput:\n")
sys.stdout.write("  No. of regions:  {0} x {1} x {2}\n".format(n_reg_x, n_reg_y, n_reg_z))
sys.stdout.write("  Region size:     {0} x {1} x {2}\n".format(size_x, size_y, size_z))
sys.stdout.write("  RAM per region:  {0:.2f} GB (x 2.3)\n".format(ram_per_region))


# Read template parameter file
try:
    with open(template_file) as par_file:
        template_par = par_file.readlines()
except Exception:
    sys.stderr.write("Error: Failed to read template parameter file sofia.par.\n")
    sys.exit(1)


# Create set of parameter files
sys.stdout.write("\nCreating parameter files:\n")

for z in range(n_reg_z):
    for y in range(n_reg_y):
        for x in range(n_reg_x):
            index = x + n_reg_x * (y + n_reg_y * z) + 1

            par = template_par[:]

            x1 = x_min + int(math.floor(x * (size_x - overlap_spat)))
            y1 = y_min + int(math.floor(y * (size_y - overlap_spat)))
            z1 = z_min + int(math.floor(z * (size_z - overlap_spec)))
            x2 = x1 + size_x
            y2 = y1 + size_y
            z2 = z1 + size_z

            if (x1 < 0):
                x1 = 0
            if (y1 < 0):
                y1 = 0
            if (z1 < 0):
                z1 = 0
            if (x2 > x_max or x == n_reg_x - 1):
                x2 = x_max
            if (y2 > y_max or y == n_reg_y - 1):
                y2 = y_max
            if (z2 > z_max or z == n_reg_z - 1):
                z2 = z_max

            i = substr_search(par, "input.region")
            if (i < 0):
                par.append("input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z1, z2))
            else:
                par[i] = "input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z1, z2)

            i = substr_search(par, "output.filename")
            if (i < 0):
                par.append("output.filename  =  {0}_{1:03d}\n".format(db_run_name, index))
            else:
                par[i] = "output.filename  =  {0}_{1:03d}\n".format(db_run_name, index)

            # Dump parameters into new file
            filename = "{0}/sofia_{1:03d}.par".format(output_dir, index)
            sys.stdout.write("  {0}\n".format(filename))
            try:
                with open(filename, "w") as par_file:
                    for item in par:
                        par_file.write("{0}".format(item))
            except Exception:
                sys.stderr.write("Error: Failed to write output parameter file: {}\n".format(filename))
                sys.exit(1)


# Create config files
sys.stdout.write("\nCreating config files and scripts:\n")


# config.ini
sys.stdout.write("  config.ini\n")
content = []
content.append("[SoFiAX]\n")
content.append("db_hostname={0}\n".format(config["database"]["db_host"]))
content.append("db_name={0}\n".format(config["database"]["db_name"]))
content.append("db_username={0}\n".format(config["database"]["db_user"]))
content.append("db_password={0}\n\n".format(config["database"]["db_pass"]))
content.append("sofia_execute=0\n")
content.append("sofia_path={0}\n".format(config["path"]["sofia2_executable"]))
content.append("sofia_processes={0:d}\n\n".format(n_reg_x * n_reg_y * n_reg_z))
content.append("run_name={0}\n".format(db_run_name))
content.append("spatial_extent={0:d},{1:d}\n".format(tolerance_spat[0], tolerance_spat[1]))
content.append("spectral_extent={0:d},{1:d}\n".format(tolerance_spec[0], tolerance_spec[1]))
content.append("flux={0:d}\n".format(tolerance_flux))
content.append("uncertainty_sigma={0:d}\n".format(tolerance_pos))

try:
    with open("{0}/config.ini".format(output_dir), "w") as config_file:
        for item in content:
            config_file.write("{0}".format(item))
except Exception:
    sys.stderr.write("Error: Failed to write file: config.ini\n")
    sys.exit(1)


# sofiax.sh
sys.stdout.write("  sofiax.sh\n")
content = []
content.append("#!/bin/bash\n\n")
content.append("#SBATCH --job-name=sofiax\n")
content.append("#SBATCH --output={0}/logs/sofiax_output_%j.log\n".format(output_dir))
content.append("#SBATCH --error={0}/logs/sofiax_error_%j.log\n".format(output_dir))
content.append("#SBATCH -N 1 # nodes\n")
content.append("#SBATCH -n 1 # tasks\n")
content.append("#SBATCH -c {0:d} # CPUs per node\n".format(n_cpu_cores))
content.append("#SBATCH --mem={0:d}G\n\n".format(ram_per_node))
content.append("module load openssl/default\n")
content.append("module load python/3.7.4\n\n")
content.append(
    "singularity exec -B /mnt:/mnt /mnt/shared/wallaby/apps/singularity/SoFiAX/sofiax.sif sofiax -c $1 -p $2\n"
)

try:
    with open("{0}/sofiax.sh".format(output_dir), "w") as config_file:
        for item in content:
            config_file.write("{0}".format(item))
except Exception:
    sys.stderr.write("Error: Failed to write file: sofiax.sh\n")
    sys.exit(1)


# run_sofiax.sh
sys.stdout.write("  run_sofiax.sh\n")
content = []
content.append("#!/bin/bash\n")
content.append("param=( ")
for i in range(n_reg_x * n_reg_y * n_reg_z):
    content.append("{0:03d} ".format(i + 1))
content.append(")\n")
content.append("for i in \"${param[@]}\"\n")
content.append("do\n")
content.append("    sbatch {0}/sofiax.sh {0}/config.ini {0}/sofia_$i.par\n".format(output_dir))
content.append("done\n")

try:
    with open("{0}/run_sofiax.sh".format(output_dir), "w") as config_file:
        for item in content:
            config_file.write("{0}".format(item))
except Exception:
    sys.stderr.write("Error: Failed to write file: run_sofiax.sh\n")
    sys.exit(1)


sys.stdout.write(
    "\nPlease check all output files before executing the run_sofiax.sh script\nto launch the SoFiAX run.\n\n"
)

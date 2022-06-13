#!/usr/bin/env python3

import sys
import math
import argparse
import configparser
from astropy.io import fits
from astropy.wcs import WCS


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to template s2p_setup configuration file"
    )
    parser.add_argument(
        "--image_cube",
        type=str,
        required=True,
        help="Path to image cube file"
    )
    parser.add_argument(
        "--region",
        type=str,
        required=False,
        help="RA and Dec boundaries for the region of the cube that will be processed",
    )
    parser.add_argument(
        "--run_name",
        type=str,
        required=True,
        help="Name for source finding run"
    )
    parser.add_argument(
        "--sofia_template",
        type=str,
        required=True,
        help="Source finding (sofia) parameter file template"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        required=True,
        help="Output directory for cubelet parameter files",
    )
    parser.add_argument(
        "--products_dir",
        type=str,
        required=True,
        help="Output directory for cubelet catalog and product files",
    )
    args = parser.parse_args(argv)
    return args


def substr_search(input_list, substr):
    """Function to search for substring in list

    """
    for i, s in enumerate(input_list):
        if substr in s:
            return i
    return -1


def get_boundary(wcs, region_str, freq):
    """Get the pixel coordinates corresponding to a RA/Dec boundary
    Expect region as "RA_min, RA_max, Dec_min, Dec_max" string
    Returns x_min, x_max, y_min, y_max

    """
    region = [float(s) for s in region_str.strip(' ').split(',')]
    p_min = wcs.all_world2pix(region[0], region[2], 1.0, freq, 0)
    p_max = wcs.all_world2pix(region[1], region[3], 1.0, freq, 0)

    # TODO(austin): check these are correct
    x_min = math.ceil(p_max[0])
    x_max = math.ceil(p_min[0])
    y_min = math.ceil(p_min[1])
    y_max = math.ceil(p_max[1])
    return x_min, x_max, y_min, y_max


def subcube_split_dimension(w_total, w, min_overlap=256):
    """Determine number of sub-cubes (N) and overlap (overlap) for optimal split.

    """
    # Determine optimal number of sub-cubes
    N = 1
    width = N * (w - min_overlap) + min_overlap
    if width > w_total:
        raise Exception("Selected region smaller than minimum sub-cube size in one dimension.")
    while width < w_total:
        N += 1
        width = N * (w - min_overlap) + min_overlap

    # Solve for minimum overlap
    overlap = math.floor((N * w - w_total) / (N - 1))
    return N, overlap


def main(argv):
    """Script to generate cubelet parameter files for source finding run.

    """
    args = parse_args(argv)

    # Read settings from configuration file
    config = configparser.ConfigParser()
    success = config.read(args.config)

    if (len(success) == 0):
        sys.stderr.write("Error: Failed to read config file: s2p_setup.ini\n")
        sys.stderr.write("Please ensure that a copy of that file is in the current directory.\n")
        sys.exit(1)

    # Default region sizes
    x_min = int(config["boundary"]["x_min"])
    x_max = int(config["boundary"]["x_max"])
    y_min = int(config["boundary"]["y_min"])
    y_max = int(config["boundary"]["y_max"])
    z_min = int(config["boundary"]["z_min"])
    z_max = int(config["boundary"]["z_max"])
    w_x, w_y, w_z = tuple([int(v) for v in config["region"]["subcube_shape"].split(" ")])
    min_overlap_spat = int(config["region"]["min_overlap_spat"])
    min_overlap_spec = int(config["region"]["min_overlap_spec"])
    tolerance_pos = int(config["pipeline"]["tolerance_pos"])
    tolerance_flux = int(config["pipeline"]["tolerance_flux"])
    tolerance_spat = (config["pipeline"]["tolerance_spat"]).split(",")
    tolerance_spec = (config["pipeline"]["tolerance_spec"]).split(",")
    tolerance_spat = [int(x) for x in tolerance_spat]
    tolerance_spec = [int(x) for x in tolerance_spec]

    # Try to open FITS file
    with fits.open(args.image_cube) as hdu:
        # Extract header information
        header = hdu[0].header

        # Get boundary
        freq = header['CRVAL4']
        wcs = WCS(header)

        if args.region:
            sys.stdout.write("Using provided RA/Dec region")
            x_min, x_max, y_min, y_max = get_boundary(wcs, args.region, freq)
            sys.stdout.write(f"RA/Dec region {args.region} in pixels: {(x_min, x_max, y_min, y_max)}")

        bitpix = int(header["BITPIX"])
        word_size = int(abs(bitpix) / 8)
        if (bitpix > 0):
            word_size = 4  # Assume 32-bit for integer arrays
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

    n_reg_x, overlap_x = subcube_split_dimension(x_max - x_min, w_x, min_overlap_spat)
    n_reg_y, overlap_y = subcube_split_dimension(y_max - y_min, w_y, min_overlap_spat)
    n_reg_z, overlap_z = subcube_split_dimension(z_max - z_min, w_z, min_overlap_spec)

    if (n_reg_x * n_reg_y * n_reg_z > 999):
        sys.stderr.write(
            """
            Error: The setup would create {0:d} parallel runs!\n
            Please reconsider your settings to reduce that number.\n
            """.format(n_reg_x * n_reg_y * n_reg_z))
        sys.exit(1)

    ram_per_region = w_x * w_y * w_z * word_size / (1024.0 * 1024.0 * 1024.0)  # GB
    sys.stdout.write("\nOutput:\n")
    sys.stdout.write("  No. of regions:  {0} x {1} x {2}\n".format(n_reg_x, n_reg_y, n_reg_z))
    sys.stdout.write("  Region size:     {0} x {1} x {2}\n".format(w_x, w_y, w_z))
    sys.stdout.write("  RAM per region:  {0:.2f} GB (x 2.3)\n".format(ram_per_region))

    # Read template parameter file
    try:
        with open(args.sofia_template) as par_file:
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

                x1 = max((w_x - overlap_x) * x + x_min, 0)
                y1 = max((w_y - overlap_y) * y + y_min, 0)
                z1 = max((w_z - overlap_z) * z + z_min, 0)
                x2 = min(x1 + w_x, x_max)
                y2 = min(y1 + w_y, y_max)
                z2 = min(z1 + w_z, z_max)

                i = substr_search(par, "input.data")
                if (i < 0):
                    par.append("input.data  =  {0}\n".format(args.image_cube))
                else:
                    par[i] = "input.data  =  {0}\n".format(args.image_cube)

                i = substr_search(par, "output.directory")
                if (i < 0):
                    par.append("output.directory  =  {0}\n".format(args.products_dir))
                else:
                    par[i] = "output.directory  =  {0}\n".format(args.products_dir)

                i = substr_search(par, "input.region")
                if (i < 0):
                    par.append("input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z1, z2))
                else:
                    par[i] = "input.region  =  {0:d},{1:d},{2:d},{3:d},{4:d},{5:d}\n".format(x1, x2, y1, y2, z1, z2)

                i = substr_search(par, "output.filename")
                if (i < 0):
                    par.append("output.filename  =  {0}_{1:03d}\n".format(args.run_name, index))
                else:
                    par[i] = "output.filename  =  {0}_{1:03d}\n".format(args.run_name, index)

                # Dump parameters into new file
                filename = "{0}/sofia_{1:03d}.par".format(args.output_dir, index)
                sys.stdout.write("  {0}\n".format(filename))
                try:
                    with open(filename, "w") as par_file:
                        for item in par:
                            par_file.write("{0}".format(item))
                except Exception:
                    sys.stderr.write("Error: Failed to write output parameter file: {}\n".format(filename))
                    sys.exit(1)

    # config.ini
    sys.stdout.write("\nCreating config files and scripts:\n")
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
    content.append("run_name={0}\n".format(args.run_name))
    content.append("spatial_extent={0:d},{1:d}\n".format(tolerance_spat[0], tolerance_spat[1]))
    content.append("spectral_extent={0:d},{1:d}\n".format(tolerance_spec[0], tolerance_spec[1]))
    content.append("flux={0:d}\n".format(tolerance_flux))
    content.append("uncertainty_sigma={0:d}\n".format(tolerance_pos))

    try:
        with open("{0}/config.ini".format(args.output_dir), "w") as config_file:
            for item in content:
                config_file.write("{0}".format(item))
    except Exception:
        sys.stderr.write("Error: Failed to write file: config.ini\n")
        sys.exit(1)


if __name__ == '__main__':
    main(sys.argv[1:])

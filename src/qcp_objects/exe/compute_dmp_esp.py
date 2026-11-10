#!/usr/bin/env python
from qcp_objects.operations.calculate_esp_dmp import calculate_esp_dmp

def main(grid, dmp, ranks, ofile=None):
    calculate_esp_dmp(grid=grid, multipoles=dmp, ranks=ranks, ofile=ofile)


def main_par():
    # Create parser
    description="""
    This script calculates the ESP energy for a given moment file and grid file
    """
    epilog=None
    prog=None
    import argparse as ap ; par=ap.ArgumentParser(description=description, prog=prog, epilog=epilog, formatter_class=ap.RawDescriptionHelpFormatter)
    # add 
    add=par.add_argument
    add("--mom",     help="CamCasp/ORIENT/JSON moment file")
    add("--grid",    help="CamCasp/ORIENT/JSON grid file")
    add("--ranks",   help="Which ranks to consider")
    add("-o",'--ofile',   help="Name of outputfile")
    # Parse
    args=par.parse_args()
    MOM_fis=args.mom
    GRD_fis=args.grid
    ranks=args.ranks
    ofile=args.ofile

    main(GRD_fis, MOM_fis, ranks, ofile=ofile)

if __name__ == '__main__':
    main_par()


#!/usr/bin/env python
from qc_objects.operations.calculate_esp_dmp import calculate_esp_dmp

def main(grid, dmp, ranks):
    calculate_esp_dmp(grid=grid, multipoles=dmp, ranks=ranks)

if __name__ == '__main__':
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
    # Parse
    args=par.parse_args()
    MOM_fis=args.mom
    GRD_fis=args.grid
    ranks=args.ranks

    main(GRD_fis, MOM_fis, ranks)


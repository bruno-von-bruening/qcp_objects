#!/usr/bin/env python


from qcp_objects.operations.compare_esp import compare_esp

def main(try_map, ref_map):
    compare_esp(try_map, ref_map)

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
    add("--trial",  type=str,   help="Map that is trial")
    add("--ref", type=str,   help=f"Map that is the refernces")
    # Parse
    args=par.parse_args()
    try_map=args.trial
    ref_map=args.ref

    main(try_map, ref_map)

if __name__ == '__main__':
    main_par()


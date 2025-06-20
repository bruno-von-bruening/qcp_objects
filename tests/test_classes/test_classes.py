#!/usr/bin/env python
import os, glob

# create classes instances (from dummy data)
qcprop_home=os.path.realpath('../..')
assert os.path.isdir(qcprop_home)

qcprop_test_data=os.path.join(qcprop_home,'tests', 'data')
assert os.path.isdir(qcprop_test_data)

def get_dummy_file(file_name, source=qcprop_test_data):
    dummy_file=os.path.join(source, file_name)
    assert os.path.isfile(dummy_file), f"Dummy file cannot be found: {dummy_file}"
    return dummy_file

def create_multipoles():
    from qc_objects.objects.properties import multipoles
    dummy_file=get_dummy_file('H2.mom')
    mom=multipoles(dummy_file)

def create_solution():
    from qc_objects.objects.basis import molecular_radial_basis
    dummy_file=get_dummy_file('H2_sol.json')
    molecular_radial_basis(dummy_file)

if __name__=="__main__":
    def do_test(script, args, tag):
        try:
            if isinstance(args, type(None)):
                script()
            else:
                raise Exception()
        except Exception as ex:
            raise Exception(f"Failure in testing {tag}: {ex}")
    test=[
        (create_multipoles, None, 'multipoles'),
        (create_solution, None, 'solution'),
    ]
    for t in test:
        do_test(*t)
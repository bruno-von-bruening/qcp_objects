
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

from qc_objects.operations.calculate_esp_dmp import calculate_esp_dmp
def test_dmp_esp():
    mom_file=get_dummy_file('H2.mom')
    grid_file=get_dummy_file('H2_idsurf.json')
    calculate_esp_dmp(multipoles=mom_file, grid=grid_file)


def do_test(script, args, tag):
    try:
        if isinstance(args, type(None)):
            script()
        else:
            raise Exception()
    except Exception as ex:
        raise Exception(f"Failure in testing {tag}: {ex}")
if __name__=="__main__":
    test=[
        (test_dmp_esp, None, 'dmp_esp'),
    ]
    for t in test:
        do_test(*t)
from . import *

from typing import List
from qcp_objects.objects.properties import multipoles as mom_cl

# import sys,os; sys.path.insert(1,os.environ["PY_MODULES"])
# import itertools as it, copy, re, shutil, pickle, json
# # Special                                                                                                                         
# import mod_import as m_imp; imp=m_imp.imp()
# import mod_objects as m_obj
# names=["numpy", "pandas"]
# np, pd  = imp.pkg(names)
# # Custom                                                                                                                          
# import mod_database as m_db, mod_io as m_io, mod_encoder as m_enc, mod_geometry as m_geo

@validate_call
def calculate_esp_dmp(multipoles: file, grid: file, ranks: List[int]|str='all_available',
                      ofile=None):
    """"""

    # Read into py-objects 
    mom=mom_cl(multipoles)
    grid=isod_grid_cl(grid)
    # Evaluate the moments
    total_grid, individual_grids=mom.eval_on_grid(grid, ranks=ranks)
    final_file=total_grid.to_file(onam='test.map.json' if ofile is None else ofile,format='json')

    results=dict(
        files={
            'map':final_file
        },
        results={
            'statistical_indicators':total_grid.gen_statistical_indicators()
        }
    )
    # Drop the results
    results_file='results.yaml'
    with open(results_file, 'w') as wr:
        yaml.safe_dump(results,wr)
    print(f"Dropped results to {results_file}")
    return results_file

    # return_dis=[]
    # for k in esp_map_di.keys():
    #     esp_map_di[k].to_file(onam={'rank':k, 'lead':lead_tag})
    #     esp_map_di[k].gen_statistical_indicators()
    #     di={    'file':esp_map_di[k].output_file,
    #             'iso_density_value': esp_map_di[k].iso_density_value,
    #             'rank': k,
    #     }
    #     return_dis.append(di)
    # return return_dis

def make_esp_dmp_file_di(MOM_fis, GRD_fis, mol_nam=None, db_nam=None, dir=None, lmax=None):
    ''' Takes moment files and grid files and will generate the ESP-DMP files.
    These files are saved in a json object (list of dictionaries) whose name is returned.
    If mol_nam is provided this entry will be added to each individual dicionary.
    If db_nam is provided the database will be named as such. Otherwise a default name will be given.
    '''
    # Get a list of dictionaries each describing the calaculation
    if type(MOM_fis)==dict:
        tuple_list=[ (k,v) for k,v, in MOM_fis.items() ]
    elif type(MOM_fis) == list:
        tuple_list=[ (None, v) for v in MOM_fis ]
    else:
        raise Exception(f"File type handling for MOM_fis not implemented: {type(MOM_fis)}")
    file_dis=[]
    for method_key, mom_fi in tuple_list:
        for grid_fi in GRD_fis:
            sub_file_dis= \
                calculate_esp_dmp(grid=grid_fi, dmp=mom_fi, lmax=lmax, lead_tag=None) # lead_tag will be read from mopment
            if type(method_key)!=type(None):
                [ f.update({'method':method_key,}) for f in sub_file_dis]
            file_dis+=sub_file_dis

    # add the molecule name if its there
    if type(mol_nam)!=type(None):
        [ f.update({'molecule_name': mol_nam, }) for f in file_dis ]

    # save as object and return object file name
    if type(db_nam)==type(None):
        db_nam='ESP-DMP_file_table.json'
    if type(dir)!=type(None):
        os.path.join(dir, db_nam)
    with open(db_nam, 'w') as wr:
        json.dump( file_dis, wr)
    return db_nam





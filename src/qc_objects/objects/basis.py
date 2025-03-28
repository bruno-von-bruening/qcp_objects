#from . import *
from . import *
import json
import numpy as np
from qc_objects.file_readers.camcasp_reader import read_camcasp_solution, read_camcasp_basis

class radial_atomic_basis(BaseModel):
    exponent_scales         : List[List[float]]
    exponent_orders         : List[List[float]]
    normalizations          : List[List[float]]
    length_of_contracted    : int|None=None
    #
    center                  : List[Tuple[float,float,float]]|None=None
    coefficients            : List[float]|None=None
    initial_coefficients    : List[float]|None=None

    #angular_functions       : None=None
    def check(self): # should check if all arrays have the same length
        raise Exception(f"Not implemented yet")

class radial_atomic_basis(radial_atomic_basis):
    """ Basis set for an individual site
    or more generically an specific atom type 
    """
    @validate_call
    def __init__(self, inp:str|dict|radial_atomic_basis, inp_type:str|None=None):
        # Make dict and then initialize
        if isinstance(inp,dict):
            dic=inp
        elif isinstance(inp,radial_atomic_basis):
            dic=inp.model_dump()
        elif isinstance(inp_type, str):
            if inp_type.upper()=='CAMCASP':
                inp_dic=read_camcasp_basis(inp).dic
                update_keys=[ 'exponent_scales', 'normalizations' ]
                dic={}
                for k in update_keys:
                    if not k in inp_dic.keys():
                        raise Exception(f"Key \'{k}\' is not in avaialable keys returned by {read_camcasp_basis}: {list(inp_dic.keys())}")
                    else:
                        # We need to reshape to accomodate contracted functions
                        dic.update({k:[]})
                        assert isinstance(inp_dic[k], list), f"Expected list for {k} in {inp}, got: {inp_dic[k]}"
                        for x in inp_dic[k]:
                            if isinstance(x, (float,int)):
                                x=[x]
                            elif isinstance(x, (list, tuple)):
                                pass
                            else: raise Exception(f"Found unkown type in {k} when reading {inp} via {read_camcasp_basis}: {x} ({type(x)})")
                            dic[k].append(x)
                # Camcasp does not give exponent orders but all basis sets are gaussians
                dic.update({'exponent_orders': [[2] for x in dic['exponent_scales']]})
            else:
                raise Exception(f"Type {inp_type} not implemented")
        else:
            raise Exception(f"")
        super().__init__(**dic)
        self.length_of_contracted=len(self.exponent_orders)
    
    def generate_horton_basis(self):
        """ Generate dictionary needed for horton basis set definition """
        # Do not need normalization (implicit in horton) and no coefficients (since basis set)
        horton_mapper={
            'exponents':self.exponent_scales,
            'exponents_orders':self.exponent_orders,
            'initials':np.array( [ 1./len(self.exponent_scales) for x in self.exponent_scales ] ).tolist(),
        }
        return horton_mapper
        

#class radial_basis(base_basis):
#    def __init__(self, inp, inp_type=None):
#        super().__init__(inp, inp_type=inp_type)
#    def to_dict(self,keys=None, format=None):
#        return self.model_dump()
        # at={
        #     'exponent_scales': self.exponent_scales,
        #     'exponent_orders': self.exponent_orders,
#        #     'coefficients': self.coefficients,
#        #     'initial_coefficients':self.initial_coefficients,
#        #     'normalizations': self.normalizations,
        # }
        # new_at={}
        # if type(keys)==type(None):
        #     new_at=at
        # else:
        #     for key in keys:
        #         if key not in at.keys():
        #             raise Exception(f"Requestest key {key} not avaiable")
        #         else:
        #             new_at.update({key: at[key]})
        # if isinstance(format, type(None)):
        #     pass

        # cleaned_dic={}
        # for key in new_at.keys():
        #     if isinstance(new_at[key], type(np.array([])) ):
        #         cleaned_dic.update({key:new_at[key].tolist()})
        #     else:
        #         cleaned_dic.update({key:new_at[key]})
        # return cleaned_dic

class molecular_basis(BaseModel):
    atom_types      : List[str]|None=None
    basis_expansions: List[radial_atomic_basis]
    centers         : List[Tuple[float,float,float]]
    num_sites       : int=-1
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
        self.num_sites=len(self.basis_expansions)
        self.check_integrity()

    def check_integrity(self):
        try:
            assert self.num_sites>0, f"Number of sites is not positive but {self.num_sites}"

            if not self.atom_types is None:
                assert len(self.atom_types)==self.num_sites
            assert len(self.basis_expansions)==self.num_sites
            if not self.centers is None:
                assert len(self.centers)==self.num_sites

        except Exception as ex:
            raise Exception(f"Failed to validate integrity of {molecular_basis}: {ex}")
class molecular_basis(molecular_basis):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.check_integrity()
    def to_json(self, output_file=None):
        dic=self.to_dict()
        cleaned_dic={}
        for key in dic.keys():
            if isinstance(dic[key], type(np.array([])) ):
                cleaned_dic.update({key:dic[key].tolist()})
            else:
                cleaned_dic.update({key:dic[key]})
        class_name=self.__class__.__name__
        if type(output_file)==type(None):
            output_file=class_name+f".json"
        with open(output_file, 'w') as wr:
            json.dump(cleaned_dic, wr)
        print(f"Printed {class_name} to {output_file}")
        return output_file
class molecular_radial_basis(molecular_basis):
    def __init__(self,inp):
        if isinstance(inp, str):
            if inp.endswith('.ISA'):
                dic=self.isa_to_dict(inp)
            elif inp.upper().endswith('.JSON'):
                dic=self.read_json(inp)
            else:
                raise Exception(f"Unkown file: {inp}")
        elif isinstance(inp, dict):
            pass
        else:
            raise Exception()
        super().__init__(**dic)
    def isa_to_dict(self, isa_file):
        dic=read_camcasp_solution(isa_file).dic 
        atom_types=dic['atom_types']
        the_coordinates=dic['coordinates'] # coordinates interferes with method in inherited class

        atoms=dic['atoms']
        basis_centers=[]
        for at in atoms:
            dic={}
            [ dic.update({k:at[k]}) for k in ['exponent_scales','normalizations','exponent_orders','coefficients'] ]
            atom=radial_atomic_basis(dic)
            basis_centers.append(atom)

        return dict(
            basis_expansions=basis_centers,
            atom_types=atom_types,
            centers=the_coordinates,
        ) 
    def read_json(self, json_file):
        try:
            with open(json_file,'r') as rd:
                data=json.load(rd)
            coor_key='coordinates'
            atty_key='atom_types'
            shape_key='shape_functions'
            #num_keys=['exponent_scales','exponent_orders','normalizations','coefficients']
    

            coordinates=data[coor_key]
            atom_types=data[atty_key]
            shape=data[shape_key]

            # Initially this were just plain lists but we want list of list to enable contracted types
            for s in shape:
                keys=["exponent_orders", "exponent_scales",'normalizations']
                for k in keys:
                    if not isinstance(s[k][0], list):
                        s[k]=[ [ss] for ss in s[k]]
                initials=s['initials']
                del s['initials']
                s.update({'initial_coefficients':initials})


            return dict(
                centers=coordinates,
                atom_types=atom_types,
                basis_expansions=[ radial_atomic_basis(x) for x in shape],
            )
        except Exception as ex:
            raise Exception(f"Cannot process json file {json_file}: {str(ex)}")

    #def from_dict(self,dic):
    #    mandatory_keys=['atom_types', 'coordinates', 'basis_expansions']
    #    available_keys=mandatory_keys
    #    provided_keys=list(dic.keys())
    #    self.basis_expansions=[]
    #    for key in provided_keys:
    #        if not key in available_keys:
    #            raise Exception(f"Provided unkown key {key}")
    #        else:
    #            if key in ['atom_types', 'coordinates']:
    #                if key=='coordinates':
    #                    the_key='centers'
    #                else:
    #                    the_key=key
    #                setattr(self,the_key,dic[key])
    #            elif key in ['shape_functions']:
    #                for shape_func in dic[key]:
    #                    self.basis_expansions.append( radial_atomic_basis(shape_func))
    #    raise Exception(self)
    def to_dict(self):
        """ Save in Json format """
        dic={} ; atoms=[]
        for i,b in enumerate(self.basis_expansions):
            at_dic=b.to_dict()
            atoms.append(at_dic)
        dic={'shape_functions':atoms, 'coordinates':self.coordinates, 'atom_types':self.atom_types}
        cleaned_dic={}
        for key in dic.keys():
            if isinstance(dic[key], type(np.array([])) ):
                cleaned_dic.update({key:dic[key].tolist()})
            else:
                cleaned_dic.update({key:dic[key]})
        return cleaned_dic
    
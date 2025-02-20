#from . import *
import json
import numpy as np

class base_basis():
    """ Basis set for an atom type """
    exponent_scales=None
    exponent_orders=None
    normalizations=None
    angular_functions=None
    coefficients=None
    initial_coefficients=None
    def __init__(self, inp, inp_type=None):
        if isinstance(inp,dict):
            self.read_from_dict(inp)
        elif isinstance(inp_type, str):
            if inp_type.upper()=='CAMCASP':
                inp_dic=read_camcasp_basis(inp).dic
                update_keys=[ 'exponent_scales', 'exponent_orders', 'normalizations' ]
                for k in update_keys:
                    setattr(self, k, inp_dic[k])
            else:
                raise Exception(f"Type {inp_type} not implemented")
        else:
            raise Exception(f"Have to provide string for \'inp_type\' atm {type(inp)}")
    def read_from_dict(self,inp):
        necessary_keys=['exponent_scales','normalizations','exponent_orders','coefficients']
        available_keys=necessary_keys+['initials','initial_coefficients']
        provided_keys=list(inp.keys())
        line=""
        for key in provided_keys:
            if not key in available_keys:
                line+=f"key \'{key}\' is not in available keys\n"
        if line!='':
            raise Exception(line+f"Available keys are {available_keys}")
        for key in necessary_keys:
            if not key in provided_keys:
                raise Exception(f"key {key} not under keys of input dict {provided_keys}")
            else:
                setattr(self, key, inp[key])
class radial_basis(base_basis):
    def __init__(self, inp, inp_type=None):
        super().__init__(inp, inp_type=inp_type)
    def to_dict(self,keys=None, format=None):
        at={
            'exponent_scales': self.exponent_scales,
            'exponent_orders': self.exponent_orders,
            'coefficients': self.coefficients,
            'initial_coefficients':self.initial_coefficients,
            'normalizations': self.normalizations,
        }
        new_at={}
        if type(keys)==type(None):
            new_at=at
        else:
            for key in keys:
                if key not in at.keys():
                    raise Exception(f"Requestest key {key} not avaiable")
                else:
                    new_at.update({key: at[key]})
        if isinstance(format, type(None)):
            pass

        cleaned_dic={}
        for key in new_at.keys():
            if isinstance(new_at[key], type(np.array([])) ):
                cleaned_dic.update({key:new_at[key].tolist()})
            else:
                cleaned_dic.update({key:new_at[key]})
        return cleaned_dic

class molecular_basis():
    atom_types=None
    basis_centers=None
    _coordinates=None
    num_sites=None
    @property
    def coordinates(self):
       return self._coordinates
    @coordinates.setter
    def coordinates(self,value):
        value=np.array(value)
        shape=value.shape
        if len(shape)!=2:
            raise Exception()
        if shape[1]!=3:
            raise Exception(f"List of 3D coordinates, got:{value}")
        self._coordinates=value
    def __init__(self):
        pass
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
        super().__init__()
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
        self.from_dict(dic)
    def isa_to_dict(self, isa_file):
        dic=read_camcasp_solution(inp).dic 
        atom_types=dic['atom_types']
        teh_coordinates=dic['coordinates'] # coordinates interferes with method in inherited class

        atoms=dic['atoms']
        basis_centers=[]
        for at in atoms:
            dic={}
            [ dic.update({k:at[k]}) for k in ['exponent_scales','normalizations','exponent_orders','coefficients'] ]
            atom=radial_basis(dic)
            basis_centers.append(atom)

        return dict(
            basis_centers=basis_centers,
            atom_types=atom_types,
            coordinates=coordinates,
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
            return dict(
                coordinates=coordinates,
                atom_types=atom_types,
                shape_functions=shape,
            )
        except Exception as ex:
            raise Exception(f"Cannot process json file {json_file}: {ex}")

    def from_dict(self,dic):
        mandatory_keys=['atom_types', 'coordinates', 'shape_functions']
        available_keys=mandatory_keys
        provided_keys=list(dic.keys())
        self.basis_centers=[]
        for key in provided_keys:
            if not key in available_keys:
                raise Exception(f"Provided unkown key {key}")
            else:
                if key in ['atom_types', 'coordinates']:
                    setattr(self,key,dic[key])
                elif key in ['shape_functions']:
                    for shape_func in dic[key]:
                        self.basis_centers.append( radial_basis(shape_func))
    def to_dict(self):
        """ Save in Json format """
        dic={} ; atoms=[]
        for i,b in enumerate(self.basis_centers):
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
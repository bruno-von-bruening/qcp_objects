#from . import *
from . import *
import json
import numpy as np
from qc_objects.file_readers.camcasp_reader import read_camcasp_solution, read_camcasp_basis


class radial_atomic_basis(BaseModel):
    exponent_scales         : List[List[float]]
    exponent_orders         : List[List[float]]
    normalizations          : List[List[float]]
    #length_of_contracted    : int|None=None
    #
    center                  : List[Tuple[float,float,float]]|None=None
    coefficients            : List[float]|None=None
    initial_coefficients    : List[float]|None=None

    #angular_functions       : None=None
    def check(self): # should check if all arrays have the same length
        raise Exception(f"Not implemented yet")
    
    @property
    def length_of_contracted(self):
        length_keys=['exponent_scales','exponent_orders','normalizations']
        if self.coefficients  is not None:
            length_keys+=['coefficients']
        if self.initial_coefficients is not None:
            length_keys+=['initial_coefficients']
        
        lengths=[ len(getattr(self,k)) for k in length_keys]
        assert len(set(lengths))==1, f"Multiple different lengths {lengths} {length_keys}"
        return lengths[0]

class radial_atomic_basis(radial_atomic_basis):
    """ Basis set for an individual site
    or more generically an specific atom type 
    """
    @validate_call
    def __init__(self, inp:str|dict|radial_atomic_basis=None, inp_type:str|None=None, **kwargs):
        # Make dict and then initialize
        if isinstance(inp,dict):
            dic=inp
            assert len(kwargs)==0
        elif isinstance(inp,radial_atomic_basis):
            dic=inp.model_dump()
            assert len(kwargs)==0
        elif isinstance(inp_type, str):
            assert len(kwargs)==0
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
        elif inp is None:
            dic=kwargs
        else:
            raise Exception(f"")
        super().__init__(**dic)
        #self.length_of_contracted=len(self.exponent_orders)
    
    def generate_horton_basis(self):
        """ Generate dictionary needed for horton basis set definition """
        # Do not need normalization (implicit in horton) and no coefficients (since basis set)
        horton_mapper={
            'exponents':self.exponent_scales,
            'exponents_orders':self.exponent_orders,
            'initials':np.array( [ 1./len(self.exponent_scales) for x in self.exponent_scales ] ).tolist(),
        }
        return horton_mapper
class radial_atomic_grid(BaseModel):
    start: float=0
    end: float=10
    units:  Literal['au','ang','angstrom']='au'
    number_of_points: int=Field(200, gt=1)
    @property
    def grid_points(self):
        return np.linspace( self.start, self.end, num=self.number_of_points)
# Add evaluation
class radial_atomic_basis(radial_atomic_basis):
    @val_call
    def evaluate_basis_functions(self, 
        grid_points: np.ndarray|List[float]|radial_atomic_basis=radial_atomic_grid()
    ):
        ''' Takes a shapefunction and returns radial for given grid_points '''
        ''' Optionally evaluates only for given indices (allows evaluation of pGTO) '''

        # Helper
        @val_call
        def eval_gaussian(
            grid_points:np.ndarray, alpha:float,Norm:float=1., order:float=2,  coeff:float=1.
        ):
            exponents= -alpha * np.power(grid_points, order)
            thres=20
            return [ coeff*Norm * np.exp(x) if x<thres else 0 for x in exponents ]
        @val_call
        def eval_contracted_gaussian(
            grid_points:np.ndarray, 
            alphas:Annotated[List[float],Len(min_length=1)], 
            Norms:List[float], 
            orders:List[float]|None=None, 
            coeff:float=1.
        ):
            assert len(alphas)==len(Norms)
            if orders is None:
                orders=[ 2 for i in range(len(alphas))]
            else:
                assert len(orders)==len(alphas)
            
            contracted_function=np.sum([ 
                eval_gaussian(grid_points, alpha=alpha, Norm=Norm, order=order) for alpha, Norm, order in zip(alphas, Norms, orders )
            ], axis=0)
            return contracted_function
        
        # Grid points should be numpy arrayz
        if isinstance(grid_points, radial_atomic_grid): grid_points=grid_points.grid_points
        grid_points=np.array(grid_points) 
        
        length=self.length_of_contracted
        def get_parameters(index):
            keys=['exponent_scales','normalizations','exponent_orders']#,'coeffiecients']
            return [ getattr(self, k)[index] for k in keys ]
        input=[ get_parameters(i) for i in range(length) ]
        
        # Evaluation of coefficientless gaussian for every function
        basis_functions= np.array([ eval_contracted_gaussian(grid_points, *inp) for inp in input ])
        return basis_functions

    @val_call
    def evaluate(self,
        grid_points: np.ndarray|List[float]|radial_atomic_basis=radial_atomic_grid()
    ):
        coefficients=self.coefficients
        if coefficients is None:
            raise AssertionError(f"Requested evaluation of shape functions without coefficients being set")

        # Grid points should be numpy arrayz
        if isinstance(grid_points, radial_atomic_grid): grid_points=grid_points.grid_points
        grid_points=np.array(grid_points) 

        # Get basis function and multiply them with the coefficients
        basis_functions=self.evaluate_basis_functions(grid_points)
        evaluated_function=np.sum(np.array(coefficients)[:,None]*basis_functions, axis=0)
        return grid_points.tolist(),evaluated_function.tolist()
        

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
    def __init__(self,inp=None, **kwargs):
        if isinstance(inp, str):
            assert len(kwargs)==0
            if inp.endswith('.ISA'):
                dic=self.isa_to_dict(inp)
            elif inp.upper().endswith('.JSON'):
                dic=self.read_json(inp)
            else:
                raise Exception(f"Unkown file: {inp}")
        elif isinstance(inp, dict):
            assert len(kwargs)==0
            pass
        elif inp is None:
            assert len(kwargs)>0
            dic=kwargs
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
    
class molecular_radial_basis(molecular_radial_basis):
    def evaluate_shape_functions(self):
        return [ x.evaluate() for x in self.basis_expansions]
    
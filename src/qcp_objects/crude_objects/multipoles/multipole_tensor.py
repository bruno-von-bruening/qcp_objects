from . import *
from qcp_objects.math_objects import allowed_representation
from qcp_objects.math_objects.multipoles import  num_elements_per_rank
from qcp_objects.math_objects.multipoles import key_mapper, triangular_index, symmetrize_tensor

class multipoles_base(myBaseModel):
    """ Delivers a single tensor of given rank"""

    rank: int
    tensor_elements: List[float]
    symmetrized: bool = False
    representation: allowed_representation

    def get_shape(self, rank, symm:bool|None=None):
        if symm is None:
            symm=self.symmetrized
        
        if not symm:
            shape= [3]*rank
        else:
            shape= [num_elements_per_rank(rank,self.representation, symm=True)]*rank
        return shape

    @property
    def matrix(self):
        return np.array(self.tensor_elements).reshape( self.get_shape(self.rank) )

    @val_call
    def from_dict(self,di:Dict[str,float], symm:bool=False):
        """ Assumes a dict with """
        r=list(set( len(k) for k in di.keys() ))
        assert len(r)==1, f"from_dict only works for single rank multipole tensors got ranks {r}:\n{di.keys()}"
        rank=r[0]
        
        if symm:
            di=symmetrize_tensor(di)
            mat=np.zeros(num_elements_per_rank(rank,'C', symm=True))
            for k,v in di.items():
                idx=key_mapper(k)
                tri_idx=triangular_index(*sorted( [ idx[:],idx[::-1] ] , reverse=True)[0]) # get large to small (lower triangular!)
                mat[tri_idx]=v
        else:
            mat=np.zeros(self.get_shape(rank,symm=False))
            for k,v in di.items():
                idx=key_mapper(k)
                mat[tuple(idx)]=v
        mat=mat.flatten().tolist()
        
        return dict(
            rank=rank,
            tensor_elements=mat,
            symmetrized=symm,
            representation='C'
        )

    def __init__(self, *args, **kwargs):
        if len(args)>0:


            assert len(args)==1, "Only single positional argument allowed"
            di=args[0]
            assert isinstance(di,Dict), "If positional argument is given, it must be a dict"
            
            assert not 'tensor_elements' in kwargs, "If positional dict argument is given, tensor_elements cannot be given as kwarg"
            symm=kwargs.get('symmetrized', False)
            load_di=self.from_dict(di, symm=symm)
            for k,v in load_di.items():
                if k in kwargs.keys():
                    assert kwargs[k]==v, f"Conflict in loading from dict and kwarg for {k}: {kwargs[k]} vs {v}"
                kwargs.update( {k:v} )

        super().__init__(**kwargs )

class molecular_multipoles_base(myBaseModel):
    ranks: List[int]
    tensor_elements: List[float]
    symmetrized: List[bool]|bool = False
    representation: List[allowed_representation]|allowed_representation
    
    def from_dict(self,di:Dict[str,float], symmetrize:bool=True):
        """ Assumes a dict with shape 'xx':'value' """

        outdi=dict(
            ranks=[],
            tensor_elements=[],
            symmetrized=[],
            representation=[]
        )

        # gather by rank:
        ranks=list(set( len(k) for k in di.keys() ))
        for r in ranks:
            # gather all keys of this rank
            subdi={ k:v for k,v in di.items() if len(k)==r }
            # build tensor for this rank
            mom=multipoles_base(subdi, rank=r, symmetrized=symmetrize)
            outdi['ranks'].append(r)
            outdi['tensor_elements']+=mom.tensor_elements
            outdi['symmetrized'].append(mom.symmetrized)
            outdi['representation'].append(mom.representation)

        return outdi

    def check_num_els(self):
        excpected_num_els=[
            num_elements_per_rank(r, rep, symm=symm) for r,rep,symm in zip(
                self.ranks,
                self.representation if isinstance(self.representation,list) else [self.representation]*len(self.ranks),
                self.symmetrized if isinstance(self.symmetrized,list) else [self.symmetrized]*len(self.ranks)
            )
        ]
        assert len(self.tensor_elements)==sum(excpected_num_els), f"Number of tensor elements {len(self.tensor_elements)} does not match expected {sum(excpected_num_els)} for ranks {self.ranks} with symmetrization {self.symmetrized} and representation {self.representation}"  
        

    def __init__(self, *args, **kwargs):
        if len(args)>0:
            assert len(args)==1, "Only single positional argument allowed"
            di=args[0]
            assert isinstance(di,Dict), "If positional argument is given, it must be a dict"
            
            assert not 'tensor_elements' in kwargs, "If positional argument is given, tensor_elements cannot be given as kwarg"
            outdi=self.from_dict(di, symmetrize=True)

            for k,v in outdi.items():
                assert not k in kwargs, f"If positional argument is given, {k} cannot be given as kwarg"
                kwargs.update( {k:v} )

        super().__init__(**kwargs )
        for k in ['symmetrized','representation']:
            if isinstance( getattr(self,k), list):
                unique_els=set( getattr(self,k) )
                if len(unique_els)==1:
                    setattr(self,k, unique_els.pop() )
        
        self.check_num_els()










    def symmetrize(self):
        raise NotImplementedError("Symmetrization not implemented yet")
        symmetrized = True

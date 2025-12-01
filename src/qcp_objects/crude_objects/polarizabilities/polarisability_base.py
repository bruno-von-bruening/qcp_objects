from ...objects.polarisabilities import *
from . import *

from qcp_objects.math_objects.multipoles import num_elements_per_rank, triangular_index
from itertools import permutations

allowed_representation=Literal['S','C']


class polarizability_single(myBaseModel):
    """
    single moment like dipole quadrupole, etc.
    """
    matrix: np.ndarray # where the values are stored
    # what if perturbation is not traceless (I will need matrix applied before)
class polarizability_single(polarizability_single):
    def traceless(self) -> polarizability_single:

        # make matrix traceless (in case not dipolar)
        mat=self.matrix
        ndim=len(mat.shape)
        if ndim==1:
            new_mat=mat
        elif ndim==2:
            if mat.shape[0]==1:
                new_mat=mat
            elif mat.shape[0]==mat.shape[1]:
                new_mat=mat - np.eye(mat.shape[0]) * np.trace(mat)/3
            else: raise Exception()
        
        # create new instance with updated matrix
        data=self.model_dump()
        data['matrix']=mat
        return polarizability_tensor(**data)

    def model_dump(self) -> dict:
        data=self.model_dump(exclude=['matrix'])
        data.update(matrix=self.matrix.tolist())
        return data
            



# v0.1.04 added direct provision of AC shift

class polarizability_tensor(myBaseModel):
    """
    Representation: Either list specifying each rank or string indicating all ranks
    PHYSICAL MEANING: Matrix of Induced mutlipole elements for external multipolar field elements
    """
    # plain array with all matrix elements
    tensor_elements: List[float] # plain elements

    # Formatting
    induced_ranks : List[int]
    induced_representation: allowed_representation
    field_ranks : List[int]
    field_representation: allowed_representation
    symmetrized : bool=True

    number_of_elements: int|None=None # for conveneince

@val_call
def polarizabilities_from_dict(
    inp:Dict[str,float], symmetrize=True
) -> dict:

    def key_mapper(key):
        if key=='0':
            return []
        else:
            mapper=dict( 
                (k,i) for i,k in enumerate(['x','y','z'])
            )
            return [ mapper[k] for k in key ]

    @val_call
    def flat_key(key:str, ranks:List[int])->int:
        """ identifies the position of key e.g. (xx) within a flatt and concatenated multipole tensor"""
        try:
            ints=key_mapper(key)

            rank=len(ints)
            if symmetrize and len(ints)>1:
                local_ind = triangular_index(*reversed(ints)) # provided as upper triangular but I store as lower!
            else:
                local_ind = sum(ints[i] * 3**(rank-i-1) for i in range(rank))  # first should be small!
            ind=sum(
                num_elements_per_rank(i,rep='C',symm=symmetrize) for i  in [ 
                    r for r in ranks if r<rank 
                ]
            )+local_ind
        except Exception as ex:
            raise Exception(f"Could not make flat key for {key} with ranks {ranks}, {ints}: {ex}") from ex

        return ind


    def cartesian_description(key)->Tuple[str,str]:
        reg=r"[x-z,0]+@[x-z,0]+"
        assert len(key)>2, f"Expected at least three key {reg}, got {key}"
        assert re.fullmatch(reg, key), f"Got key {key}, expected something of form {reg}"
        return key.split('@')

    pairs=sorted([ (cartesian_description(k),v) for k,v in inp.items()], key=lambda x: x[0] )


    if symmetrize:
        blocks=[]
        for (a,b),v in pairs:
            blocks+=[ dict(
                ind=len(key_mapper(a)), 
                ind_val=a,
                field=len(key_mapper(b)),
                field_val=b,
                val=v,
            )]
        # first clear all duplicates in the blocks
        def small_symmetrize(blocks, inner, outer):
            new_blocks=[]
            thres=1.e-8
            for j in sorted(list(set([ b[inner] for b in blocks]))):
                # continue if nothing to symmetrize
                sel_inner=[ b for b in blocks if b[inner]==j]
                if j<2:
                    new_blocks+=sel_inner
                else:
                    for i in sorted(list(set([ b[outer] for b in sel_inner]))):
                        sel=[ b for b in sel_inner if b[outer]==i]

                        for outer_val in list(set([ the_sel[f"{outer}_val"] for the_sel in sel ])):
                            uniques={}
                            small_sel=[ the_sel for the_sel in sel if the_sel[f"{outer}_val"]==outer_val]
                            for s in small_sel:
                                #choices=sorted( [s[f"{inner}_val"], ''.join(reversed(list(s[f"{inner}_val"])))], key=lambda x:list(x) )
                                base = s[f"{inner}_val"]
                                # enumerate all unique permutations that have the same counts of x,y,z
                                choices = sorted(
                                    {''.join(p) for p in permutations(base)},
                                    key=lambda x: list(x)
                                )
                                identified=choices[0]


                                if not identified in uniques.keys():
                                    found=[ ss for ss in small_sel if ss[f"{inner}_val"] in choices ]
                                    if len(found)>0:
                                        vals=[ ss['val'] for ss in found ]
                                        avg=sum(vals)/len(vals)
                                        difs=[ abs(v-avg) for v in vals ]
                                        assert all( d<thres for d in difs), f"Difference larger {thres}: {difs}"
                                        val=avg
                                    
                                    else: raise Exception(found)

                                    uniques.update({ identified : s})
                                    uniques[identified].update({
                                        'val':val, f"{inner}_val":identified,
                                    })
                            new_blocks+=list( uniques.values() )
                
            return new_blocks
        blocks=small_symmetrize(blocks, 'field','ind')
        blocks=small_symmetrize(blocks, 'ind','field')
        pairs=[
            ((b['ind_val'],b['field_val']), b['val']) for b in blocks
        ]





    ranks_ind, ranks_field=[
        list(set(len(p[0][i]) for p in pairs)) for i in range(2)
    ]
    shape=[
        sum( 
            num_elements_per_rank(i,rep='C',symm=symmetrize) for i  in ranks
        )
        for ranks in [ranks_ind,ranks_field]
    ]


    pairs_in_int=[]



    for p,v in pairs:
        try:
            pairs_in_int.append((
                tuple([ flat_key(p[i], ranks) 
                        for i,ranks in enumerate([ranks_ind, ranks_field])
                ]), v)
            )
        except Exception as ex:
            raise Exception(f"Could not process pair {p} with value {v}: {ex}\n{pairs}") from ex


    try:
        mat=np.zeros(shape)
        
        for (i,j),v in pairs_in_int:
            if not (j,i) in [ x[0] for x in pairs_in_int]:
                mat[i,j]=v
            #    pass
            #else:
            #    found=[ x for x in pairs_in_int if x[0]==(j,i)]
            #    assert len(found)==1
            #    v=(v+found[0][1])/2
            mat[i,j]=v
    except Exception as ex: raise Exception(f"Could not build matrix of shape {shape}: {ex}") from ex




        # symmetrize
    


    # @val_call
    # def symmetrize(self,thres=1.e-10) -> polarizability_tensor:
    #     mat=self.mat

    mat_flat=mat.reshape(-1).tolist()

    return dict( # format of polarisability tensor
        tensor_elements=mat_flat,
        symmetrized=symmetrize,
        induced_ranks=ranks_ind,
        induced_representation='C',
        field_ranks=ranks_field,
        field_representation='C',
    )
class polarizability_tensor(polarizability_tensor):
    def __init__(self, *args,**kwargs):
        if len(args)==1:
            if isinstance(args[0],dict):
                assert len(kwargs)==0, f"Provided dictionary as input but also additional kwargs!"
                kwargs=polarizabilities_from_dict(args[0])
        elif len(args)>0:
            raise Exception(f"More than 1 positional argument!")

        super().__init__(**kwargs)
        self.number_of_elements=int(np.sum( np.array(kwargs['tensor_elements']).reshape(-1).size ))
    def elements_per_rank(self,representation, rank,symm):
        return num_elements_per_rank(rank, rep=representation, symm=symm)
    def representation_vector_length(self,representation, ranks, symm):
        return np.sum([ self.elements_per_rank(representation, n, symm) for n in ranks ])
    @property
    def expected_shape(self):
        induced_length, field_length=[
            self.representation_vector_length(rep, rank, self.symmetrized)
            for rep,rank in [
                (self.induced_representation, self.induced_ranks) , (self.field_representation, self.field_ranks)
            ]
        ]
        return (induced_length, field_length)


    @property
    def matrix(self):
        induced_length, field_length=self.expected_shape
        return np.array(self.tensor_elements).reshape(induced_length, field_length)

    def check_validity(self)->bool:
        """ Check that the number of elements is consistent with ranks and representations """
        expected_shape=self.expected_shape
        expected_elements= expected_shape[0]*expected_shape[1]
        if expected_elements!=self.number_of_elements:
            raise Exception(f"Inconsistent number of elements: Expected {expected_elements} but got {self.number_of_elements}")
        return True
    def make_traceless(self):
        # go over blocks
        mat=self.matrix
        for j in range(mat.shape[1]):
            vec=mat[:,j]
            start=0
            for r in self.induced_ranks:
                length=self.representation_vector_length(self.induced_representation, [r], self.symmetrized)
                if r<2:
                    pass
                elif r==2:
                    indices=[ triangular_index(i,i) for i in range(3) ]
                    trace=np.sum([vec[ind+start] for ind in indices])
                    for ind in indices:
                        vec[ind+start]-=trace/3
                elif r==3:
                    warn("Traceless for rank 3 skipped for now")
                    # for ind in np.indices( 3*[length]):
                    #     a,b,c=ind
                    #     if a<=b and b<=c:
                    #         pass
                    #     else:
                    #         ind=triangular_index(a,b,c)
                    #         # Formula is 
                    #         vec[ind+start]-= np.linalg.norm([a,b,c])**2 * (
                    #             (c if a==b else 0)+
                    #             (b if a==c else 0)+
                    #             (a if b==c else 0) 
                    #         )/5
                elif r==4:
                    warn("Traceless for rank 4 skipped for now")
                else:
                    raise NotImplementedError(f"Traceless for rank {r} not implemented yet")
                start+=length
            mat[:,j]=vec

        mod=self.model_dump()
        mod.update(tensor_elements=mat.reshape(-1).tolist())
        return polarizability_tensor( **mod)

    
class list_of_polarizability_tensors(myBaseModel):
    freqsq: List[float]|None=None
    tensor_elements: List[List[float]]

    induced_ranks : List[List[int]]
    induced_representation: List[allowed_representation]|allowed_representation
    field_ranks : List[List[int]]
    field_representation: List[allowed_representation]|allowed_representation
    freqsq: List[float]|None=None

    number_of_elements: int|None=None# for conveneince
    @property
    def available_frequencies(self):
        return sorted(list(self.freqsq), reverse=True)
    @val_call
    def filter_freqsq(self, filter_values: float | List[float]):
        if not isinstance(filter_values, list): filter_values=[ filter_values]
        length=len(self.tensor_elements)

        not_found=[ x for x in filter_values if x not in self.available_frequencies]
        assert len(not_found)==0, f"Following requested frequencies are not vailable {not_found}"

        keys=[ x for x in list_of_polarizability_tensors.model_fields if not x in ['number_of_elements'] ]
        selection=[ i for i in range(length) if self.freqsq[i] in filter_values ]
        reduced=dict([ (k,[]) for k in keys])

        data=self.model_dump()
        for k in ['induced_representation','field_representation']:
            if isinstance(data[k], str): data[k]=[ data[k] for i in range(length)]
        for k in keys:
            try:
                vals= data[k]
                reduced[k]+= [vals[i] for i in selection]
            except Exception as ex: raise Exception(f"{k}: {ex}")
        
        return list_of_polarizability_tensors(**reduced)
    def model_dump(self,**kwargs):
        dd=super().model_dump(exclude='number_of_elements')
        dd.update({'number_of_elements': int(self.number_of_elements)})
        return dd
    def __init__(self, *args,**kwargs):
        super().__init__(*args,**kwargs)
        self.number_of_elements=int(np.sum( np.array(self.tensor_elements).reshape(-1).size ))
        if len(set(self.induced_representation))==1: self.induced_representation=self.induced_representation[0]
        if len(set(self.field_representation))==1: self.field_representation=self.field_representation[0]

    @classmethod
    @val_call
    def from_list(cls, tensors: List[polarizability_tensor], **kwargs):
        if not tensors:
            raise ValueError("No polarizability_tensor objects provided.")
        # Collect all fields except number_of_elements
        keys = [k for k in tensors[0].model_fields if k != 'number_of_elements']
        data = {k: [] for k in keys}
        for t in tensors:
            for k in keys:
                data[k].append(getattr(t, k))
        # Handle representation fields: flatten if all are the same
        for k in ['induced_representation', 'field_representation']:
            if k in data and len(set(data[k])) == 1:
                data[k] = data[k][0]

        for k,v in data.items():
            assert k not in kwargs.keys(), f"Found key {k} both in provided tensors and kwargs"

        return cls(**data,**kwargs)
    @val_call
    def split_to_centers(self) -> List[polarizability_tensor]:
        """ Split into a list of tensors each describing a site"""
        keys=[ k for k in polarizability_tensor.model_fields if k not in ['number_of_elements'] ]
        data=self.model_dump(include=keys)
        length=len(self.tensor_elements)
        for k in ['induced_representation','field_representation']:
            if isinstance(data[k], str): data[k]=[ data[k] for i in range(length)]

        centers=[]
        for i in range(length):
            centers+= [ polarizability_tensor(**dict( (k,data[k][i]) for k in keys)) ]
        return centers




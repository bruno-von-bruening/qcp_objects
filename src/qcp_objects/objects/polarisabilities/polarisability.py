from . import *

from qcp_objects.file_readers.file_readers import pot_file_reader

from ...crude_objects.polarizabilities.polarisability_base import allowed_representation, polarizability_tensor

class Molecular_Polarisability():
    pass


# read C12.pot
@val_call
def read_pot_file(file:pdtc_file):
    reader = pot_file_reader(file_path=file)
    return reader


from .. import *



class Pairwise_Polarizabilities(BaseModel):
    """
    Representation: Either list specifying each rank or string indicating all ranks
    PHYSICAL MEANING: Matrix of Induced mutlipole elements for external multipolar field elements
    """
    # plain array with all matrix elements
    tensor_elements: List[List[float]] # plain elements

    # Formatting
    induced_ranks : List[List[int]]
    induced_representation: List[allowed_representation]|allowed_representation
    field_ranks : List[List[int]]
    field_representation: List[allowed_representation]|allowed_representation
    freq2: List[float]|None=None
    site_indices:        List[Tuple[int,int]]

    number_of_elements: int|None=None# for conveneince


from qcp_objects.file_readers.camcasp_reader import read_isapol_pair_file
class Pairwise_Polarizabilities(Pairwise_Polarizabilities):
    def __init__(self, *args,**kwargs):
        if len(args)>0:
            assert len(args)==1, f"Only one positional element!"
            the_file=args[0]
            assert os.path.isfile(the_file)

            reader=read_isapol_pair_file(the_file)
            dic=reader.model_dump()

            # Check for collisions with kwargs
            for k in dic.keys():
                if k in kwargs.keys():
                    if kwargs[k]==dic[k]: del dic[k]
                    else: raise Exception(f"found key {k} both from kwargs and from {the_file} with different vavlues: {kwargs[k]} {dic[k]}")
            kwargs.update(**dic) 
        super().__init__(**kwargs, number_of_elements=0)
        self.number_of_elements=np.sum( np.array(kwargs['tensor_elements']).reshape(-1).size )
    @val_call
    def split_to_centers(self) -> List[polarizability_tensor]:
        """ Split into a list of tensors each describing a site"""
        keys=[ k for k in polarizability_tensor.model_fields if k not in ['number_of_elements'] ]
        data=self.model_dump(include=keys)
        length=len(self.site_indices)
        for k in ['induced_representation','field_representation']:
            if isinstance(data[k], str): data[k]=[ data[k] for i in range(length)]

        centers=[]
        for i in range(length):
            centers+= [ polarizability_tensor(**dict( (k,data[k][i]) for k in keys)) ]
        return centers

    @property
    def available_frequencies(self):
        return sorted(list(self.freq2), reverse=True)
    @val_call
    def filter_freq2(self, filter_values:float|List[float]) -> Pairwise_Polarizabilities:
        if not isinstance(filter_values, list): filter_values=[ filter_values]
        length=len(self.site_indices)

        not_found=[ x for x in filter_values if x not in self.available_frequencies]
        assert len(not_found)==0, f"Following requested frequencies are not vailable {not_found}"

        keys=[ x for x in self.model_fields if not x in ['number_of_elements'] ]
        selection=[ i for i in range(length) if self.freq2[i] in filter_values ]
        reduced=dict([ (k,[]) for k in keys])

        data=self.model_dump()
        for k in ['induced_representation','field_representation']:
            if isinstance(data[k], str): data[k]=[ data[k] for i in range(length)]
        for k in keys:
            try:
                vals= data[k]
                reduced[k]+= [vals[i] for i in selection]
            except Exception as ex: raise Exception(f"{k}: {ex}")
        
        return Pairwise_Polarizabilities(**reduced)
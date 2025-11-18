from . import *
from .polarisability_base import list_of_polarizability_tensors, polarizability_tensor

class Pairwise_Polarizabilities_Base(list_of_polarizability_tensors, myModel):
    #  should execute the init of myModel):
    # 
    """
    Representation: Either list specifying each rank or string indicating all ranks
    PHYSICAL MEANING: Matrix of Induced mutlipole elements for external multipolar field elements
    """
    site_indices:        List[Tuple[int,int]] # Pairs of sites
    
    def __init__(self, *args, **kwargs):
        myModel.__init__(self, *args, **kwargs)
        list_of_polarizability_tensors.__init__(self, *args, **kwargs)
    
# TODO make for both localised and pairwise!
    @classmethod
    @val_call
    def from_list(cls, tensors: list_of_polarizability_tensors|List[polarizability_tensor],**kwargs):
        if not isinstance(tensors, list_of_polarizability_tensors):
            return super().from_list(tensors, **kwargs)
        else:
            return super().__init__(**tensors.model_dump(),**kwargs)

class Localized_Distributed_Polarizabilities_Base(list_of_polarizability_tensors ,myModel):
    """ """
    site_tags: List[str]

    def __init__(self, *args, **kwargs):
        myModel.__init__(self, *args, **kwargs)
        list_of_polarizability_tensors.__init__(self, *args, **kwargs)
    @classmethod
    @val_call
    def from_list(cls, tensors: list_of_polarizability_tensors|List[polarizability_tensor],**kwargs):
        if not isinstance(tensors, list_of_polarizability_tensors):
            return super().from_list(tensors, **kwargs)
        else:
            return super().__init__(**tensors.model_dump(exclude=['number_of_elements']),**kwargs)

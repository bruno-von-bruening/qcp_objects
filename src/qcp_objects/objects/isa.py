from . import *

class Radial_Weights(BaseModel):
    coordinates_of_sites    : List[Tuple(float,float,float)]
    types_of_sites          : List[str]

    ## what about contractions?
    functions_per_site  : List[inr] # how many functions belong to each site
    decay_factors       : List[float] # exponent factor 
    decay_orders        : List[float] # exponent order (normally 1 or 2)
    normalizations      : List[float]

    coefficients        : List[float]   | None =None
    initial_coefficients : Listr[float] | None =None

    #def __init__(self, *args, **kwargs):
    #    try:
    #        for k, i in kwargs.items():
    #            if isinstance(i, (list, np.ndarray)):
    #                kwargs[k]=drop_to_string(i)
    #        super().__init__(*args,**kwargs)
    #    except Exception as ex:
    #        raise Exception(f"Problem in initializting \'ISA_Weights_Base\': {ex}")
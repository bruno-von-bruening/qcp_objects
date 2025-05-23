
#from . import *
import json
import os,sys, re
import numpy as np

from pydantic import validate_call, BaseModel, Field, PlainSerializer, BeforeValidator; val_call=validate_call(config=dict(arbitrary_types_allowed=True))
from typing import List, Annotated, Tuple, Literal
from annotated_types import Len
def nd_array_before_validator(x):
    # custom before validation logic
    if isinstance(x, List):
        x = np.array(x)
    return x

def nd_array_serializer(x):
    # custom serialization logic
    return x.tolist()
    # return np.array2string(x,separator=',', threshold=sys.maxsize)

from pydantic import ConfigDict, ValidationError
mynpar = Annotated[   np.ndarray,
                            BeforeValidator(nd_array_before_validator),
                            PlainSerializer(nd_array_serializer, return_type=List),
                        ]
from qc_utils.file_handling import load_json_or_yaml
from qc_global_utilities.encoding_and_conversion.encoding import nuclear_charge_to_element_symbol, element_symbol_to_nuclear_charge
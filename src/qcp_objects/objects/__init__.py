
#from . import *
import json
import os,sys, re
import numpy as np
from logging import getLogger
log = getLogger(__name__)

from qcp_global_utils.encoding_and_conversion.constants import ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM
from qcp_global_utils.encoding_and_conversion.encoding import element_symbol_to_nuclear_charge, nuclear_charge_to_element_symbol

from pydantic import validate_call, BaseModel, Field, PlainSerializer, BeforeValidator; val_call=validate_call(config=dict(arbitrary_types_allowed=True))
from typing import List, Annotated, Tuple, Literal, Union
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

from pydantic import ConfigDict, ValidationError, validate_call; val_call=validate_call(config=dict(arbitrary_types_allowed=True))
mynpar = Annotated[   np.ndarray,
                            BeforeValidator(nd_array_before_validator),
                            PlainSerializer(nd_array_serializer, return_type=List),
                        ]
from qcp_global_utils.environment.file_handling import load_json_or_yaml
from qcp_global_utils.encoding_and_conversion.encoding import nuclear_charge_to_element_symbol, element_symbol_to_nuclear_charge

allowed_length_units=['BOHR','ANGSTROM']
class my_units(BaseModel):
    LENGTH: Literal[tuple(allowed_length_units)]='BOHR'

from qcp_objects.file_readers.file_readers import fchk_reader , read_xyz_file
from qcp_objects.file_readers.file_writers import xyz_writer, txyz_writer

file=str
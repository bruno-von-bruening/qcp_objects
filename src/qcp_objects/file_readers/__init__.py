import re,os
import numpy as np

from pydantic import validate_call, BeforeValidator
from typing import List, Tuple

from pydantic import Field, BaseModel
from qc_global_utilities.encoding_and_conversion.encoding import element_symbol_to_nuclear_charge, nuclear_charge_to_element_symbol
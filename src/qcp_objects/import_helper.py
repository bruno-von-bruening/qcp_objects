from pydantic import validate_call, BaseModel;
val_call=validate_call(config=dict(arbitrary_types_allowed=True))
from typing import Literal, Union, Tuple, List

import numpy as np, os, sys
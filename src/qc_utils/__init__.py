from pydantic import BeforeValidator, BaseModel, validate_call
from typing import Annotated, List, Union, Tuple
from enum import Enum
import os, subprocess as sp
import json, yaml
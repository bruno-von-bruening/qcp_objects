try: # rich only acts on printout so its entirely options
    from rich.traceback import install
    install(suppress=[ 'pydantic', 'pydantic_core' ])
except Exception as ex: pass 

import re

from pydantic import validate_call, BaseModel, ValidationError;
val_call=validate_call(config=dict(arbitrary_types_allowed=True))
from typing import Literal, Union, Tuple, List, Dict

import numpy as np, os, sys

from qcp_global_utils.pydantic.pydantic import pdtc_file
from qcp_global_utils.environment.file_handling import load_json_or_yaml, is_json_or_yaml 

class myBaseModel(BaseModel):
    class Config:
        arbitrary_types_allowed = True
        extra = 'forbid'

import logging
log=logging.getLogger(__name__)
warn=log.warning

class myModel(myBaseModel):
    def __init__(self, *args, **kwargs):
        if len(args)>1:
            raise Exception(f"Provided more than one positional argument to {self.__class__.__name__} but no such arguments are expected!")
        elif len(args)>0:
            if not hasattr(self, 'read_file'):
                raise Exception(f"Coding error, class {self.__class__.__name__} does not implement read_file method but received positional arguments")
            else:
                assert callable(self.read_file), f"Coding error, class {self.__class__.__name__} does not implement read_file method but received positional arguments"
                new_kwargs=self.read_file(args[0])

                for k,v in new_kwargs.items():
                    # If not there add it
                    if not k in kwargs.keys():
                        kwargs.update({k:v})
                    else:
                        # If there do nothing is its the same else raise warning that doubly define
                        if kwargs[k]== v: pass
                        else: 
                            warn(f"Key {k} (for {type(self).__name__}) provided both from file {args[0]} and from kwargs with different values: {v} and {kwargs[k]}. Using value from kwargs.")
        super().__init__(**kwargs)

from enum import Enum



import logging

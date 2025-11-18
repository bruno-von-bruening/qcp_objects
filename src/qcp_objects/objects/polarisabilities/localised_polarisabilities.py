from . import *

from qcp_objects.file_readers.polarisabilities import read_isapol_pol_file
from qcp_objects.crude_objects.polarizabilities.polarizability_crude import Localized_Distributed_Polarizabilities_Base

class Localized_Distributed_Polarizabilities(Localized_Distributed_Polarizabilities_Base):
    """ """
    # C 6 to C 12 coeffiecinets (Tensors)
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    @val_call
    def read_file(self, file_name:pdtc_file):
        file_reader=read_isapol_pol_file(file_name)
        obj=file_reader.read_file()
        assert issubclass(Localized_Distributed_Polarizabilities, type(obj)), f"Expected to read Localized_Distributed_Polarizabilities but got {type(obj)}"
        return obj.model_dump()



    
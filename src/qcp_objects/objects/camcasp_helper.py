from qcp_objects.util.import_helper import *
from qcp_objects.file_readers.camcasp_reader import read_camcasp_sites
class camcasp_sites(BaseModel):
    site_types: List[str]
    site_tags: List[str]
    coordinates: List[Tuple[float,float,float]]
    def __init__(self,*args,**kwargs):
        if len(args)==1:
            input=args[0]
        elif len(args)>1: raise Exception(f"More than one positional argument: {args}")
        else: input=None

        if not input is None:
            assert len(kwargs)==0
            kwargs=self.read_camcasp_input(input)

        super().__init__(**kwargs)
    def read_camcasp_input(self, fi):
        """Read camcasp sitenames input file"""
        return read_camcasp_sites(fi).dic

        

        
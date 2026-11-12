
from . import *
from qcp_global_utils.environment.file_handling import temporary_file

from qcp_objects.util.import_helper import myBaseModel

enumline=Tuple[int, str]
class enumline():
    no: int
    line: str
    def __init__(self, no:int, line:str):
        self.no=no
        self.line=line

class file_reader_2(myBaseModel):
    file_path: pdtc_file|None = None # input file
    comtag: str|None = None # 
    lines: List[enumline]|None = None # list of lines
    #_lineno_nocom: list[int]|None = None # List of number of lines without comments
    _lines_nocom: List[enumline]|None = None # List of number of lines without comments
    _lines_no_nocom: List[int]|None = None # List of number of lines without comments
    @property
    def lines_nocom(self):
        """ Get all lines without comments """
        if self._lines_nocom is None:
            if self.comtag is None:
                self._lines_nocom=self.lines
            else:
                self._lines_nocom = [ enumline(i,l) for i,l in enumerate(self.lines) if l.strip().startswith(self.comtag) ]
        return self._lines_nocom
    @property
    def lines_no_nocom(self):
        """ Get all line numbers without comments """
        if self._lines_no_nocom is None:
            self._lines_no_nocom = [ v.no for v in self.lines_nocom ]
        return self._lines_no_nocom
    
    @val_call
    def __init__(self, inp_fi:pdtc_file|None=None, **kwargs):
        if inp_fi is not None:
            assert not 'file_path' in kwargs, "Expected either inp_fi or file_path keyword argument but not both"
            kwargs.update( file_path=inp_fi )
            assert not 'lines' in kwargs, "Expected either inp_fi or lines keyword argument but not both"

            if is_json_or_yaml(inp_fi):
                raise Exception("Json or yaml files are not supported: {inp_fi}")

            # TODO
            # Take care about double definition between kwargs and file!
            with open(inp_fi) as rd:
                kwargs.update( lines=[ enumline(i,l.strip('\n')) for i,l in enumerate(rd.readlines()) ] )

            if len(kwargs['lines'])<1:
                raise Exception(f"Empty File: \"{inp_fi}\"")
        super().__init__(**kwargs)
class file_reader_2(file_reader_2):
    @val_call
    def get_line_no(self, match=None|str, fullmatch=False, single=False, in_range=None, case_sensitive=True):
        # Check that type is a string
        args= ( [] if case_sensitive else [re.IGNORECASE] )
        
        op=re.search if not fullmatch else re.fullmatch
        results=[ v.no for v in self.lines_nocom if op( match, v.line, *args ) ]

        # Get range
        if in_range is None:
            pass
        else:
            results=[ x for x in results if x>=in_range[0] and x<=in_range[1]]

        # Check number
        if single:
            if len(results)!=1:
                raise Exception(self.file_path, results, match)
            return results[0]
        else:
            return results
    @val_call
    def get_block_indices(self, start_tag:str, end_tag:str, case_sensitive:bool=True) -> List[Tuple[int,int]]:
        """ Start and end matches"""

        re_args= ( [] if case_sensitive else [re.IGNORECASE] )
        def is_start(l:str): return bool(re.match(start_tag, l, *re_args))
        def is_end(l:str): return bool(re.match(end_tag,l, *re_args))


        # Get all possible start and ends
        block_starts=[]
        block_ends=[]
        for li in self.lines_nocom:
            i,l=li.no, li.line
            if is_start(l):
                block_starts+=[i]
            if is_end(l):  # Do not do elif since appearently not everything has an end statement (at least in orient pol files)
                block_ends+=[i]
        
        # Match the starts with the ends (next larger i than start)
        block_indices=[]
        for start_index in block_starts:
            larger_ends=[x for x in block_ends if x>start_index]
            assert len(larger_ends)>0, f"Did not find end statement after line (index: {start_index}): {self.lines_nocom[start_index]}"
            end_index=min(larger_ends)
            block_indices+=[(start_index, end_index)]
        if len(block_indices)==0:
            raise Exception(f"Did not find blocks for start_tag: {start_tag} and end_tag: {end_tag} in file {self.file_path}")
        return block_indices
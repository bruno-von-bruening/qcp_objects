from . import *
from .file_reader_base import file_reader_2, enumline

from qcp_objects.crude_objects.polarizabilities.polarisability_base import (
    allowed_representation, polarizability_tensor , list_of_polarizability_tensors,
)
from qcp_objects.crude_objects.polarizabilities.polarizability_crude import (
    Pairwise_Polarizabilities_Base, Localized_Distributed_Polarizabilities_Base,
)


igncase=re.IGNORECASE

class site_helper(myBaseModel):
    index: int|None=None
    tag: str=None
    site_names: Tuple[str,str] 
    site_indices : Tuple[int,int]|int|None=None
    ranks_induced   : List[int]
    ranks_field     : List[int]
    freq: float
    represenation: allowed_representation = 'S'
    def __init__(self,format:Literal['pairwise','localised'], *args,**kwargs):
        if len(args)>0:
            assert len(args)==1, f"Only one positional element!"
            line=args[0]
            assert isinstance(line,str), f"Expected a string line but got {type(line)}"
            dic=self.read_block_header(line, format=format)
            kwargs.update(**dic)
        try:
            super().__init__(**kwargs)
        except ValidationError as ve:
            raise Exception(f"Could not fill {type(self).__name__}: {ve}")
    def read_block_header(self, line:str, format:Literal['pairwise','localised']):
        ar = list(enumerate(line.split()))[1:] # Assume the first element is the block tag

        # --- Tag identifier functions ---
        def idf_index(tag):
            return re.match(r'INDEX',tag, re.IGNORECASE)
        def idf_site_names(tag):
            return re.search(r'SITE-(NAMES|LABELS)', tag, re.IGNORECASE)
        def idf_site_indices(tag):
            return re.match(r'SITE-INDICES', tag, re.IGNORECASE)
        def idf_rank(tag):
            return re.match(r'RANK', tag, re.IGNORECASE)
        def idf_cart_spher(tag):
            return re.match(r'CARTSPHER', tag, re.IGNORECASE)
        def idf_freq2(tag):
            return re.match(r'FREQ(2|SQ)', tag, re.IGNORECASE)

        # --- Getter functions ---
        def get_index(tag, *args):
            if len(args) != 1:
                raise ValueError(f"INDEX expects 1 argument, got {len(args)}: {args}")
            val = args[0]
            def index_convention(num):
                """CamCASP has fortran indexing (starting from 1)"""
                return num-1
            try:
                return dict(index=index_convention(int(val)))
            except Exception as ex:
                raise Exception(f"Cannot interprete as int: {val}\n {line}") from ex
        def get_site_names(tag, *args):
            if len(args) != 2:
                raise ValueError(f"SITE-NAMES expects 2 arguments, got {len(args)}: {args}")
            site_one, site_two = args
            if re.match(tag,r'SITE-INDICES', igncase):
                try:
                    site_one=int(site_one)
                    site_two=int(site_two)
                except Exception as ex:
                    raise Exception(f"Cannot interprete as int: {site_one}, {site_two}\n {line}") from ex
            return dict(site_names=(site_one, site_two))
        def get_site_indices(tag, *args):
            if len(args) != 1:
                raise ValueError(f"SITE-INDICES expects 1 arguments, got {len(args)}: {args}")
            site_one = args
            try:
                site_one=int(site_one)
            except Exception as ex:
                raise Exception(f"Cannot interprete as int: {site_one}\n {line}") from ex
            return dict(site_indices=site_one)
        def get_rank(tag, *args) -> dict:
            vals = list(args)
            if len(vals)==7:
                assert format=='pairwise', f"RANK with 7 arguments only valid for pairwise format, got {format}"
                assert vals[1]==':' and vals[3]=='BY' and vals[5]==':', f"Expected ' : ',  ' BY ', ' : ' keywords in rank definition but got {vals}"
                ints=[vals[0], vals[2], vals[4], vals[6]]
                try:
                    rank_ind_start=int(ints[0])
                    rank_ind_end=int(ints[1])
                    rank_field_start=int(ints[2])
                    rank_field_end=int(ints[3])
                    return dict(
                        ranks_induced=list(range(rank_ind_start, rank_ind_end+1)),
                        ranks_field=list(range(rank_field_start, rank_field_end+1))
                    )
                except Exception as ex: raise Exception(f"Cannot interprete as int: {ints}\n {line}") from ex
            elif len(vals)==3:
                assert re.match(r'TO', vals[1], igncase), f"Expected 'TO' keyword as second argument of ranks but got {vals[1]}"
                assert format=='localised', f"RANK with 3 arguments only valid for localised format, got {format}"
                ints=[vals[0], vals[2]]
                try:
                    rank_start=int(ints[0])
                    rank_end=int(ints[1])
                    return dict(
                        ranks_induced=list(range(rank_start, rank_end+1)),
                        ranks_field=list(range(rank_start, rank_end+1))
                    )
                except Exception as ex: raise Exception(f"Cannot interprete as int: {ints}\n {line}") from ex
            else:
                raise ValueError(f"RANK expects 3 or 7 arguments, got {len(vals)}: {vals}")
        def get_representation(tag, *args):
            if len(args) != 1:
                raise ValueError(f"CARTSPHER expects 1 argument, got {len(args)}: {args}")
            return dict(representation=args[0])
        def get_freq(tag, *args):
            if len(args) != 1:
                raise ValueError(f"FREQ2 expects 1 argument, got {len(args)}: {args}")
            val = args[0]
            try:
                return dict(freq=float(val))
            except Exception as ex:
                raise Exception(f"Cannot interprete as float: {val}\n {line}") from ex

        # --- Pairing: (idf_func, get_func, attribute_name) ---
        tag_map = [
            (idf_index, get_index, 'index'),
            (idf_site_names, get_site_names, 'site_names'),
            (idf_site_indices, get_site_indices, 'site_indices'),
            (idf_rank, get_rank, 'ranks'),
            (idf_cart_spher, get_representation, 'represenation'),
            (idf_freq2, get_freq, 'freq'),
        ]

        # Find all tag positions and their mapping
        tag_positions = []
        for idx, word in ar:
            for idf_func, get_func, attr in tag_map:
                if idf_func(word):
                    tag_positions.append((idx, word, get_func, attr))
                    break

        tag_positions.sort(key=lambda x: x[0])


        kwargs = {}
        try:
            for i, (start_idx, tag, get_func, attr) in enumerate(tag_positions):
                end_idx = tag_positions[i + 1][0] if i + 1 < len(tag_positions) else len(ar)+1
                elements = [x[1] for x in ar if start_idx <= x[0] < end_idx]
                # Remove the tag itself from elements
                elements = elements[1:]
                add_kwargs = get_func(tag, *elements)
                kwargs.update(add_kwargs)
        except Exception as ex:
            raise Exception(f"Error processing line: {line}\n (identified tags were {[ar[x[0]-1][1] for x in tag_positions ]}).\nTried to fill {add_kwargs}") from ex

        start_of_args= min([ x[0] for x in tag_positions ])
        if start_of_args==2:
            kwargs.update( tag=ar[0][1] )


        return kwargs
    @val_call
    def make_pol_tensors(self, values:List[float]) -> polarizability_tensor:
        pol_tens= polarizability_tensor(
            tensor_elements=values,
            induced_ranks=self.ranks_induced,
            field_ranks=self.ranks_field,
            induced_representation=self.represenation,
            field_representation=self.represenation,
        )
        pol_tens.check_validity()
        return pol_tens


allowed_representation=Literal['S','C']
class read_isapol_pol_file(myBaseModel):
    _comtag=['!']
    # tensor_elements: List[List[float]]    # 2D <- sites x plain matrix
    # induced_ranks: List[List[int]]             
    # field_ranks: List[List[int]]
    # induced_representation: List[allowed_representation]|allowed_representation
    # field_representation: List[allowed_representation]|allowed_representation
    # site_indices: List[Tuple[int,int]]
    # freq2: List[float]

formats=Literal['pairwise','localised']
class formats(Enum):
    pairwise='pairwise'
    localised='localised'

class read_isapol_pol_file(file_reader_2):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
    
    def read_file(self, expected_format:formats|None=None):

        # There are two file formats for localised and not localised for some reason
        
        if expected_format is None:
            endfile_tags = self.get_line_no(match=r'ENDFILE')
            if len(endfile_tags)==0:
                format=formats.pairwise
            elif len(endfile_tags)>=1: # For some reasons somtimes these files use endfile instead of end
                format=formats.localised

        if format==formats.pairwise:
            start_tag, end_tag=r'POL' , r'END'
        else:
            start_tag, end_tag=r'ALPHA', r'(ALPHA|ENDFILE)'
        block_indices = self.get_block_indices(start_tag, end_tag,case_sensitive=False)

        @val_call
        def get_nums(content:List[enumline]) -> List[float]:
            nums=[]
            for c in content:
                try:
                    nums+= [ float(x) for x in c.line.split() ]
                except Exception:
                    raise Exception(f"Cannot interprete as float (on line {c.no}): {c.line}")
            return nums

        # Format
        pol_tens=[]
        freqs=[]
        sites=[]
        for bstart, bend in block_indices:
            site=site_helper( 'localised' if not format==formats.pairwise else 'pairwise', self.lines_nocom[bstart].line )
            sites+=[ site ]

            # Get the numbers
            content=self.lines_nocom[bstart+1 : bend ]
            nums=get_nums(content)

            pol_tens+=[ site.make_pol_tensors(nums) ]
            freqs+=[ site.freq ]

        list_of_pol_tens=list_of_polarizability_tensors.from_list(pol_tens, freqsq=freqs)
        if format==formats.pairwise:
            obj=Pairwise_Polarizabilities_Base.from_list(pol_tens, freqsq=freqs,site_indices=[ x.site_indices for x in sites],)
        else:
            obj=Localized_Distributed_Polarizabilities_Base.from_list(pol_tens, freqsq=freqs,site_tags=[ x.tag for x in sites],)

        return obj
        

        
        # Now make List of polarisability tensors
        # How to put to high-level?



            
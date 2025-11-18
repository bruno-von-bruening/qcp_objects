from . import *
# Constructor that only includes the Tinker like frame definition
# Actual frame which gets the actual cartesian local axis frames from the definition and an geometry

from enum import Enum

# Copy from local_axis_frame@$OWN_PACKS/multipole_scripts/modules/mod_database.py
class AXIS_TYPE(str,Enum):
    Z_ONLY      ="Z_ONLY"
    BISECTOR    ="BISECTOR"
    Z_THEN_X    ="Z_THEN_X"
    TRISECTOR   ="TRISECTOR"
    Z_THEN_BISECTOR="Z_THEN_BISECTOR"

tinker_frame_types = Union[Tuple[int,int,int],Tuple[int,int,int,int]]

def read_sulley_local_axis_file(file)->List[ tinker_frame_types ]:
    tuples = []
    with open(file, "r") as f:
        lines=f.readlines()
    for line in lines:
        parts = line.strip().split()
        if len(parts)>0:
            assert len(parts) in [3,4], f"Expected 3 or 4 entries per line but found {len(parts)} in line {line}"
            try:
                int_parts = tuple([int(x) for x in parts])
            except:
                raise Exception(f"Expected integer entries but found {parts} in line {line}")
            finally:
                tuples.append(int_parts)
    return tuples

class local_axis_frame_constructor(BaseModel):
    definition: List[ tinker_frame_types ] # atom_index, axis_type, axis_atom_index

    def __init__(self,*args,**kwargs):

        # TAKE Care about Input
        if len(args)==1:
            input=args[0]
        elif len(args)>1:
            raise Exception("local_axis_frame_constructor expects either a single argument (the input definition) or keyword arguments")
        else:
            input=None
        if 'input' in kwargs:
            input=kwargs.pop('input')

        if input is not None:
            assert isinstance(input, str), f"Expected string input but found {type(input)}"
            assert os.path.isfile(input), f"Expected a file but found {input}"

            assert not 'definition' in kwargs, "Expected either input or definition keyword argument but not both"
            # Read the file
            if any([input.endswith(x) for x in ['yaml','json','yml']]):
                kwargs.update(load_json_or_yaml(input))
            elif input.endswith('.txt'):
                tuples = read_sulley_local_axis_file(input)
                kwargs.update(dict(definition=tuples))
        super().__init__(**kwargs)

    @classmethod
    def interprete_tinker(cls,index1,index2, index3, index4=None):
        """ Interpret the Tinker local axis frame definition (up to three numbers)"""
        a=AXIS_TYPE
        if index4 is not None:
            if all(np.array([index2,index3,index4])<0):
                ty= a.TRISECTOR
            else:
                ty= a.Z_THEN_BISECTOR
        else:
            assert (index1>0), f"Expected one of {[a.Z_ONLY,a.BISECTOR,a.Z_THEN_X]} but found {index1},{index2},{index3} (first two not positive)"
            if index3 == 0:
                ty=a.Z_ONLY
            elif index3 < 0:
                assert index2<0, f"Found third index negative but second not negative (expected for {a.BISECTOR.value}): {index1},{index2},{index3}"
                ty=a.BISECTOR
            else:
                ty=a.Z_THEN_X
        return ty
    @classmethod
    def orient_single_axis(cls, axis_type: AXIS_TYPE, indices: Tuple[str,str,str,str|None]) -> str:
        """ Given the axis type, print the orientation in ORIENT format
        Expects incdices in CamCASP index (e.g. H<index>)
        """
        # SEE orient manual page 10 (but I think The stone TInker format gives local axis frames orient does not know)
        def direct(the_ax:Literal['x','y','z'], i1,i2, reverse=False):
            """ """
            return f"{the_ax} From "+ (f"{i1} to {i2}" if not reverse else f"{i2} to {i1}")
        def global_ax(the_ax:Literal['x','y','z'], the_ax2:Literal['x','y','z'], reverse=False):
            """ """
            return f"{the_ax} Global {'-' if reverse else ' '}{the_ax2}"
        def bisector(the_ax:Literal['x','y','z'], i1,i2,i3):
            return f"{the_ax} From {i1} between {i2} and {i3}"


        def orient_z_only(i1,i2):
            z=direct('z',i1,i2)
            x=global_ax('x','x')
            return [z,x]
        def orient_bisector(i1,i2,i3):
            # Z should be negative bisector but orient does not allow that
            log.warning("Orienting with bisector should be negative (-Z)")
            z=bisector('z',i1,i2,i3)
            x=direct('x',i1,i3)
            return [z,x]
        def orient_z_the_x(i1,i2,i3):
            z=direct('z',i1,i2, reverse=True)
            x=direct('x',i1,i3)
            return [z,x]
        def orient_z_then_bisector(a,b,c,d):
            z=direct('z',a,b, reverse=True)
            x=bisector('x',a,c,d)
            return [z,x]
        def orient_trisector(a,b,c,d):
            log.warning("Orienting with trisector is not following a precise definition in orient")
            x=direct('x',b,c)
            y=direct('y',c,d)
            z=direct('z',a,b)
            return [x,y,z]


        if axis_type == AXIS_TYPE.Z_ONLY:
            axs=orient_z_only(*indices[:2])
        elif axis_type == AXIS_TYPE.BISECTOR:
            axs=orient_bisector(*indices[:3])
        elif axis_type == AXIS_TYPE.Z_THEN_X:
            axs=orient_z_the_x(*indices[:3])
        elif axis_type == AXIS_TYPE.Z_THEN_BISECTOR:
            axs=orient_z_then_bisector(*indices)
        elif axis_type == AXIS_TYPE.TRISECTOR:
            axs=orient_trisector(*indices)
        else:   
            raise Exception(f"Forgot implementin axis type {axis_type.value}")

        orient_line=f"{indices[0]} "+' '.join(axs)
        return orient_line
    def orient_block(self, camcasp_indices:List[str]):
        def orient_all(self, camcasp_indices=List[str]):
            orient_lines=[]
            for ax in self.definition:
                orient_lines.append( self.orient_single_axis(
                    self.interprete_tinker(*ax),
                    tuple([(camcasp_indices[abs(i)-1] if i is not None else None) for i in ax])
                ))
            return orient_lines
        orient_block='\n'.join([
            'Axes', 
            *[ 4*' '+x for x in orient_all(self, camcasp_indices) ],
            'End'
        ])+'\n'
        return orient_block




#   Internal object accordint to Tinker definition
# class local_axis_frame(BaseModel):
#     atom_index: int
#     description: str|None=None
from . import *


coor_ty=List[Tuple[float,float,float]]
at_ty=List[str]
con_ty=List[Tuple[int,int,int]] # index implicit which atom connected to, bon_order

class writer(BaseModel):
    output_file: str
    def __init__(self,*args,**kwargs):
        if len(args)==1:
            kwargs.update({'output_file':args[0]})
        elif len(args)>1:
            raise Exception(f"Too many positional arguments: {args}")
        super().__init__(**kwargs)
    def write_to_file(self,lines: list):
        for l in lines:
            if type(l)!=str:
                raise Exception(f"Expected string but found {type(l)}: {l}")
        if type(self.output_file)==type(None):
            raise Exception(f"You need to provide a file name")
        with open(self.output_file, 'w') as wr:
            wr.writelines(''.join([x.strip('\n')+"\n" for x in lines]))



def write_xyz_file(output_file:str, coordinates:coor_ty, atom_types:at_ty, comment:str|None=None) -> List[str]:

    assert len(coordinates)==len(atom_types), f"Number of coordinates ({len(coordinates)}) must match number of atom types ({len(atom_types)})" 
    num_sites=len(coordinates)
    
    lines=[f"{num_sites}"] # first line in xyz is number of atoms
    lines+=[ '' if comment is None else comment]

    for coor, at_ty in zip(coordinates, atom_types):
        coor_form=[ f"{x:.8f}" for x in coor ]
        new_line=f"  {at_ty:<2}"+4*' '+ ' '.join([f"{x:>12}" for x  in coor_form])
        lines.append(new_line)
    
    lines=[l.strip('\n')+'\n' for l in lines]
    return lines

def write_txyz_file(output_file:str, coordinates:coor_ty, atom_types:at_ty, connectivity:con_ty, comment:str|None=None) -> List[str]:
    # assert len(connectivity)==len(coordinates), f"Number of connectivity entries ({len(connectivity)}) must match number of coordinates ({len(coordinates)})"

    # Check
    connected_indices=[ ]    
    for c in connectivity:
        connected_indices+=[c[0],c[1]]
    not_included=set([ x  for x in connected_indices if not x in list(range(len(coordinates)))])
    if len(not_included)>0:
        raise Exception(f"Connectivity indices must be between 0 and {len(coordinates)-1} (thats the coordinate range), but found {not_included}")


    crude_lines=write_xyz_file(output_file, coordinates, atom_types, comment)
    crude_lines=[ l.strip('\n') for l in crude_lines ]
    max_len_crude_line=max([len(l) for l in crude_lines])

    def get_conneciton(index):
        """ Consider that txyz starts from 1 not from zero!"""
        connections=[]
        for a,b,bond in connectivity:
            if bond>0:
                if (a==index):
                    add_me=b
                elif(b==index):
                    add_me=a
                else: add_me=None
                if not add_me is None and not add_me in connections:
                    connections.append(add_me)
                
        connections=[ x+1 for x in connections ] # txyz starts from 1 not from zero!
        return connections


    for i in range(2,len(crude_lines)):
        index=i-2
        connections=get_conneciton(index)
        crude_lines[i]=' '.join([
            f"{index+1:>4}", crude_lines[i], f"{atom_types[index]}_{index}", ' '.join([f"{x}" for x in connections])
        ])
    return crude_lines



class xyz_writer(writer):
    coordinates     : coor_ty
    atom_types      : List[str]
    comment         : str|None=None
    write           : bool=True
    def __init__(self, output_file, **kwargs):
        super().__init__(output_file, **kwargs)
        if self.write:
            lines=self.create_lines()
            self.write_to_file(lines)
    def create_lines(self):
        return write_xyz_file(self.output_file, self.coordinates, self.atom_types, self.comment)
class txyz_writer(xyz_writer):
    connectivity    : List[Tuple[int,int,int]]
    def create_lines(self):
        return write_txyz_file(self.output_file, self.coordinates, self.atom_types, self.connectivity, self.comment)

from . import *

class writer():
    output_file=None
    def __init__(self,output_file):
        if type(output_file)!=str:
            raise Exception()
        else:
            self.output_file=output_file
    def write_to_file(self,lines: list):
        for l in lines:
            if type(l)!=str:
                raise Exception(f"Expected string but found {type(l)}: {l}")
        if type(self.output_file)==type(None):
            raise Exception(f"You need to provide a file name")
        with open(self.output_file, 'w') as wr:
            wr.writelines(lines)
class xyz_writer(writer):
    coordinates=None
    atom_types=None
    comment=None
    write=True
    def __init__(self, output_file, write=True, coordinates=None, atom_types=None, comment=None):
        super().__init__(output_file)
        self.write=write
        self.coordinates=np.array(coordinates)
        if type(atom_types)==type(None):
            raise Exception()
        self.atom_types=tuple(atom_types)
        self.comment=comment
        if self.write:
            lines=self.create_lines()
            self.write_to_file(lines)
    def create_lines(self):
        if type(self.coordinates)!=type(np.array([])):
            raise Exception()
        if type(self.atom_types)!=tuple:
            raise Exception()
        
        num_sites=len(self.coordinates)
        lines=[f"{num_sites}"]
        if type(self.comment)==type(None):
            lines+=[' ']
        else:
            lines+=[ self.comment ]

        for coor, at_ty in zip(self.coordinates, self.atom_types):
            coor_form=[ f"{x:.8f}" for x in coor ]
            new_line=f"  {at_ty:<2}"+4*' '+ ' '.join([f"{x:>12}" for x  in coor_form])
            lines.append(new_line)
        
        lines=[l.strip('\n')+'\n' for l in lines]
        return lines
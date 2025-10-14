from . import *

coor_ty=List[Tuple[float,float,float]]

from rdkit.Chem import rdmolfiles, Mol, rdDetermineBonds, Mol as RDMol, rdmolops, MolFromXYZBlock

class geometry_hidden(BaseModel):
    ifi: str|None=None
    units: my_units # Controls only input/output units
    atom_types: List[str]|None=None
    multiplicity:int|None=Field(default=None, ge=1 )
    multiplicity_ion: int|None=Field(default=None, ge=1)
    charge:int|None=Field(default=None)
    adjacency_matrix:List[List[int]]|None=None
    #
    _internal_units: my_units=my_units( LENGTH='BOHR' ) # Controls only input/output units
    _coordinates: List[Tuple[float,float,float]]|None=None # allways save in Bohr
    @property
    def coordinates(self):
        ''' Coorindates are allways saved in Bohr, but are returned in current unit '''
        conversion=self.internal_to_external_unit('LENGTH', self.units.LENGTH)
        return self._coordinates*conversion
    
    def read_from_file(self,ifi):
        raise Exception(f"This is supposed to be done on higher level (since other objects (e.g. multipoles) may be provided)")
    @val_call 
    def __init__(self,
        *args,
        **kwargs,
    ):
        def read_input_file(inp):
            kwargs=self.read_from_file(inp)
            return kwargs

        def check_input(inp):
            if inp is None:
                if not ('coordinates' in kwargs.keys() and 'atom_types' in kwargs.keys()):
                    raise Exception(f"Need to provide either input file or coordinates/atom types. You provided input_file (or input file was None)")
            assert isinstance(inp, str),f"Input file must be a string, but found {type(inp)}"
            assert os.path.isfile(inp), f"Input file must be a file, but found {inp}"

        input_provided=False
        if 'input' in kwargs.keys():
            inp=kwargs.pop('input')
            input_provided=True
        elif len(args)>0:
            if len(args)==1:
                inp=args[0]
                input_provided=True
            elif len(args)>1:
                raise Exception(f"Too many positional arguments: {args}")


        if input_provided:
            assert inp is not None, f"Input file is None"
            check_input(inp)
            kwargs.update(dict(ifi=inp))
            kwargs=read_input_file(inp)
        else:
            if not ('coordinates' in kwargs.keys() or 'atom_types' in kwargs.keys()):
                raise Exception(f"No input file provided buat also no coordinates/atom types. Provided: {kwargs.keys()}")
            if not  ('length_units' in kwargs.keys()):
                raise Exception(f"Need to provide length unit when coordinates are read in over options")
            else:
                length_units=kwargs['length_units']
                assert length_units in allowed_length_units
                kwargs.update({'units':{'LENGTH':length_units}})

            atom_types=kwargs['atom_types']
            if not atom_types is None:
                tmp_at_ty=[]
                for x in atom_types:
                    if isinstance(x, str):
                        tmp_at_ty.append(x)
                    elif isinstance(x, int):
                        tmp_at_ty.append(nuclear_charge_to_element_symbol(x))
                    else:
                        raise Exception(f"Unkown type {x}")
                kwargs['atom_types']=tmp_at_ty
        
        coordinates=kwargs['coordinates']
        del kwargs['coordinates']
        super().__init__(**kwargs)
        self._coordinates=np.array(coordinates)/self.internal_to_external_unit('LENGTH',self.units.LENGTH)
    def model_dump(self):
        dic=super().model_dump()
        dic['coordinates']=self.coordinates.tolist()
        dic['length_units']=self.units.LENGTH
        del dic['units']
        del dic['_internal_units']
        del dic['_coordinates']
        return dic

class geometry_hidden(geometry_hidden):
        
    def internal_to_external_unit(self, key, input_unit=None):
        ''' Return conversion factor for given unit '''
        
        internal_units=self._internal_units.model_dump()
        external_units=self.units.model_dump()
        if not key in internal_units.keys():
            raise Exception()
        if not key in external_units.keys():
            raise Exception()
        internal_unit, output_unit= internal_units[key], external_units[key]
        
        #setattr(self.units,key,None) # Set this explicitly before requesting output, to make things clear!
        
        def unit_val(unit_key):
            # Values in au
            if unit_key.upper() in ["BOHR"]:
                ret_unit=1.
            elif unit_key.upper() in ['ANGSTROM','ANG']:
                ret_unit=ANGSTROM_TO_BOHR
            else:
                raise Exception(f"Unkown unit key: {unit_key}")
            return ret_unit

        internal_unit_val=unit_val(internal_unit)
        output_unit_val=unit_val(output_unit)
    
        conversion_factor= internal_unit_val/output_unit_val
        return conversion_factor

    def connec(self):
        def get_rdkit_mol():
            if type(self.ifi) != type(None):
                raw_mol = rdmolfiles.MolFromXYZFile(self.ifi)
                rdmol = RDMol(raw_mol)
            else:
                self.units.LENGTH='ANGSTROM'
                mol_block=f"{len(self.coordinates)}\n\n"
                for ty,coor in zip( self.atom_types, self.coordinates):
                    bohr_to_angstrom=0.529177249
                    coor=[f"{float(x)}" for x in coor]
                    mol_block+=f"{ty} {' '.join(coor)}\n"
                rdmol=MolFromXYZBlock(mol_block)
            return rdmol
        rdmol=get_rdkit_mol()
        rdDetermineBonds.DetermineBonds(rdmol, charge=0)
        self.adjacency_matrix=np.array(rdmolops.GetAdjacencyMatrix(rdmol, useBO=True))
    @property
    def camcasp_atom_tags(self):
        camcasp_atom_tags=[]
        at_count={}
        if type(self.atom_types)==type(None):
            raise Exception(f"Atom types are not defined")
        for at_ty in self.atom_types:
            if not at_ty in at_count.keys():
                at_count.update({at_ty:1})
            else:
                at_count[at_ty]+=1
            camcasp_atom_tags.append(f"{at_ty}_{at_count[at_ty]:02d}")
        return camcasp_atom_tags
    @property
    def camcasp_sites(self):
        coords=self.coordinates
        camcasp_atom_tags=self.camcasp_atom_tags
        atom_types=self.atom_types
        nuclear_charges=self.nuclear_charges
        return [
            dict(coordinate=c, tag=tag, atom_type=atty, charge=charge)                for
            c,tag,atty, charge  in zip( coords, camcasp_atom_tags, atom_types, nuclear_charges)
        ]
    @property
    def camcasp_sites_string(self):
        lines=[]
        for site in self.camcasp_sites:
            coor, tag, at_ty, Z = list(site.values())
            line=' '.join(['{:<6} {:>5}']+3*["{:>12}"]+['Type {:<2}'])
            line=line.format(tag,f"{Z:.1f}", *[f"{x:.8f}" for x in coor], at_ty)
            lines.append(line)
        return lines
    @property
    def camcasp_coordinates_string(self):
        lines=[]
        self.units.LENGTH='BOHR'
        lines+=[ "units {}".format(self.units.LENGTH) ]
        lines+=self.camcasp_sites_string

        return '\n'.join(lines)
    @property
    def nuclear_charges(self):
        return [ element_symbol_to_nuclear_charge(x) for x in self.atom_types]
    def print_out(self, output_name=None, format='xyz'):
        if format=='xyz':
            self.units.LENGTH='ANGSTROM'
            writer=xyz_writer(output_name, write=True, coordinates=self.coordinates, 
                    atom_types=self.atom_types, comment='Generated automatically')
        elif format=='txyz':
            self.units.LENGTH='ANGSTROM'
            connec_list=[]
            for i in range(len(self.adjacency_matrix)):
                for j in range(len(self.adjacency_matrix)):
                    if self.adjacency_matrix[i][j]>0:
                        connec_list.append( (i,j,self.adjacency_matrix[i][j]) )
            
            writer=txyz_writer(output_name, write=True, coordinates=self.coordinates, 
                    atom_types=self.atom_types, comment=None,
                    #f"Charge={self.charge} Multiplicity={self.multiplicity}",
                    connectivity=connec_list)
        return output_name
    def to_xyz_file(self, output_name=None):
        return self.print_out(output_name=output_name, format='xyz')
    def to_txyz_file(self, output_name=None):
        if self.adjacency_matrix is None:
            self.connec()
        return self.print_out(output_name=output_name, format='txyz')
    @validate_call
    def com(self, weighted:bool=False) -> Tuple[float, float, float]:
        coords=np.array(self.coordinates)
        charges=np.array(self.nuclear_charges)
        if weighted: charges=np.sqrt(charges)
        return np.sum([   [c/np.sum(charges)*coord] for c,coord in zip(charges,coords)   ], axis=0)[0]
        #    , axis=0
    def to_dict(self):
        dic={
            'coordinates':self.coordinates.tolist(),
            'atom_types': self.atom_types,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
            'length_units': self.units.LENGTH,
        }
        return dic
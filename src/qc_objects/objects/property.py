from . import *
from .grids import molecular_isodensity_grid
import copy
from qc_objects.file_readers.file_readers import fchk_reader
from qc_global_utilities.encoding_and_conversion.constants import ANGSTROM_TO_BOHR, BOHR_TO_ANGSTROM
from typing import Literal
from qc_objects.file_readers.camcasp_reader import read_mom_file

class my_units(BaseModel):
    LENGTH: Literal['BOHR','ANGSTROM']
    #energy: Literal['HARTREE','kJ']
class geometry(BaseModel):
    ifi: str|None=None
    coordinates: List[List[float]]|None=None # allways save in Bohr
    units: my_units # Controls only input/output units
    _internal_units: my_units=my_units( LENGTH='BOHR' ) # Controls only input/output units
    atom_types: List[str]|None=None
    multiplicity:int=Field(default=None, ge=1 )
    multiplicity_ion: int|None=Field(default=None, ge=1)
    charge:int=0
    adjacency_matrix:List[List[int]]|None=None

class geometry(geometry):
    @property
    def coordinates(self):
        ''' Coorindates are allways saved in Bohr, but are returned in current unit '''
        conversion=self.unit_handler('LENGTH')
        return self._coordinates*conversion
    @coordinates.setter
    def coordinates(self,value):
        """ Saves coordinates in Bohr """
        coor=None
        if type(value)==list:
            coor=np.array(value)
        elif isinstance(value, type(np.array([])) ):
            coor=value
        else:
            raise Exception(f"Unkown format: {type(value)}")

        shape=coor.shape
        if len(shape)!=2:
            raise Exception()
        elif shape[1]!=3:
            raise Exception()

        conversion=self.unit_handler('LENGTH')
        self._coordinates=coor/conversion
    def unit_handler(self, key):
        ''' Return conversion factor for given unit '''
        if not key in self._internal_units.keys():
            raise Exception()
        if not key in self.units.keys():
            raise Exception()
        internal_unit, output_unit= self._internal_units[key], self.units[key]
        self.units[key]=None # Set this explicitly before requesting output, to make things clear!
        
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
    
        conversion_factor= output_unit_val/internal_unit_val
        return conversion_factor
        
    def __init__(self,inp=None,coordinates=None, atom_types=None, length_units=None):
        if type(inp)==str:
            if inp.endswith('.xyz'):
                info_dict=read_xyz_file(inp)
                self.units['LENGTH']='ANGSTROM'
            elif inp.endswith('.sdf'):
                self.units['LENGTH']='ANGSTROM'
                info_dict=read_sdf_file(inp)
            elif inp.endswith('.fchk'):
                kwargs=fchk_reader(inp).make_geometry()
                kwargs.update({'units':{'LENGTH':'BOHR'}})
            else:
                raise Exception(f"File format not handeled yet: {inp}")
            # self.coordinates=info_dict.coordinates
            # self.atom_types=info_dict.atom_types
            # self.charge=info_dict.charge
            # self.multiplicty=info_dict.multiplicity
            super().__init__(**kwargs)
        elif isinstance(inp, type(None)):
            if type(coordinates) != type(None):
                if type(length_units)==type(None):
                    raise Exception(f"Need to provide length unit when coordinates are read in over options")
                self.units['LENGTH']=length_units
                self.coordinates=coordinates
            else:
                raise Exception(f"No input file provided, need coordinates then")
            if type(atom_types) != type(None) :
                tmp_at_ty=[]
                for x in atom_types:
                    if isinstance(x, str):
                        tmp_at_ty.append(x)
                    elif isinstance(x, int):
                        tmp_at_ty.append(atomic_charge_to_atom_type(x))
                    else:
                        raise Exception(f"Unkown type {x}")
                self.atom_types=tmp_at_ty

            else:
                raise Exception(f"No input file provided, need atom types then")
        else:
            raise Exception(inp)
    def connec(self):
        #from rdkit.Chem import AllChem,rdDetermineBonds
        if type(self.ifi) != type(None):
            raw_mol = rdchem.rdmolfiles.MolFromXYZFile(self.ifi)
            rdmol = rdchem.Mol(raw_mol)
        else:
            mol_block=f"{len(self.coordinates)}\n\n"
            for ty,coor in zip( self.atom_types, self.coordinates):
                bohr_to_angstrom=0.529177249
                coor=[f"{float(x)*bohr_to_angstrom}" for x in coor]
                mol_block+=f"{ty} {' '.join(coor)}\n"
            rdmol=rdchem.MolFromXYZBlock(mol_block)
            #mol_block = rdchem.MolToMolBlock(rdmol)
        rdDetermineBonds.DetermineBonds(rdmol, charge=0)
        self.adjacency_matrix=np.array(rdmolops.GetAdjacencyMatrix(rdmol, useBO=True))
        #m_io.out().warn(f"Implement check for adjacency matrix, ie total valence electrons, number of bonds per species")
        #self.adj_mat=rdchem.rdmolops.GetAdjacencyMatrix(rdmol, useBO=True)
    #@classmethod
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
    def camcasp_coordinates_string(self):
        lines=[]
        self.units.LENGTH='BOHR'
        lines+=[ "units {}".format(self.units.LENGTH) ]

        coords=self.coordinates
        camcasp_atom_tags=self.camcasp_atom_tags
        atom_types=self.atom_types
        nuclear_charges=self.nuclear_charges
        for coor, tag, at_ty, Z in zip( coords, camcasp_atom_tags, atom_types, nuclear_charges):
            line=' '.join(['{:<6} {:>5}']+3*["{:>12}"]+['Type {:<2}'])
            line=line.format(tag,f"{Z:.1f}", *[f"{x:.8f}" for x in coor], at_ty)
            lines.append(line)
        return '\n'.join(lines)
    @property
    def nuclear_charges(self):
        return [ element_symbol_to_nuclear_charge(x) for x in self.atom_types]
    def print_out(self, output_name=None, format='xyz'):
        if format=='xyz':
            self.units['LENGTH']='ANGSTROM'
            writer=xyz_writer(output_name, write=True, coordinates=self.coordinates, 
                    atom_types=self.atom_types, comment='Generated automatically')
    def to_dict(self):
        dic={
            'coordinates':self.coordinates.tolist(),
            'atom_types': self.atom_types,
            'charge': self.charge,
            'multiplicity': self.multiplicity,
        }
        return dic

class rotate_multipoles():
    """ Cartesian is in from xx,yy,zz,xy,xz,yz
    Handle multipoles internally in BOHR!
    """
    # Slow version for a quick one 
    def __init__(self, input_format=None, output_format=None):
        # implement test function
        pass
    def rotate_quick(self, array):

        BOHR_TO_ANGSTROM = 1./0.52917721067
        matrices=self.get_mat_sph_to_crt_stone()
        assert any([len(array)==x**2 for x in range(1,10)])
        max_rank, ranges=self.get_min_rank(array)
        assert len(matrices)>=max_rank, f"{len(matrices)}, {max_rank}"

        return_array=[]
        for i in range(max_rank):
            #lw,up=ranges[i]
            sub_vec=np.array([ array[x] for x in range(*ranges[i]) ])
            mat=matrices[i]
            assert mat.shape[1]==len(sub_vec), f"{mat.shape} {len(sub_vec)}"
            cart_vec=np.matmul(mat, sub_vec)*BOHR_TO_ANGSTROM**i
            return_array.append(cart_vec)
        flat_array=[]
        for x in return_array:
            flat_array+=[float(y) for y in  x ]
        return flat_array


    def get_min_rank(self, array):
        levels=[ 2*n+1 for n in range(10) ]
        lev_ind=[ sum(levels[:i+1]) for i,lev in enumerate(levels) ]
        ranges=[ (lev_ind[i],lev_ind[i+1]) for i in range(len(levels)-1)]
        ranges=[ (x[0],x[1]) for x in ranges]
        ranges=[ [0,1] ]+ranges
        found=False
        for i, lev in enumerate(lev_ind):
            if lev==len(array):
                found=True
                return i, ranges[:i+1]
        if not found:
            raise Exception(len(array), lev_ind)

        
    def get_mat_sph_to_crt_stone(self):
        """ Provides Matrices to convert spherical representation to cartesian representation in stones definition for a given multipolar rank
        see Theory of Intermolecular Forces, Stone, Second Edition (2013), p.287
        return a list of matrices, one for each rank from monopoles to hexadecapoles

        The transformtation is achieved by implementing the equations in previously mentioned book chapter.
        Spherical Moments are represented as vector indicating their order (as ordered on the previous pages)
        the cartesian moments are then expressed as vector describing the components, e.g. for qxx:
            q_xx = -1/2 q20  +  sqrt(3)/2 q22c
            q_xx = ( -1/2   0   0   sqrt(3)/2   0 ) dot (q20 q21c q21s q22c q22s)  # latter is the vector representation of observed multipoles 
        """

        def make_vector_for_index(index):
            """ Find i=index-2n+1 (n as large as possible)
            return vector ( #2(n+1)+1 )[i]
            Purpose is describing multipole levels as matrix rows (square matrix for every rank)"""
            # Size per level
            levels=[ 2*n+1 for n in range(10) ]
            # End of every level
            lev_starts=[0]+[ sum(levels[:i+1]) for i,lev in enumerate(levels) ]

            index_lev=None
            for i, level_length in enumerate(levels):
                level_end   =lev_starts[i+1]
                level_start =lev_starts[i]
                if level_end>index: # If the level shatters the index, its contained
                    sub_index= index - (level_start)
                    assert sub_index>=0 and sub_index<level_length, f"Something wrong in this function"
                    index_lev= (sub_index , level_length)
                    break
            if isinstance(index_lev, type(None)):
                raise Exception(f"Did not find the represenation for {index_lev}")
            else:
                sub_index,level_length=index_lev
                represenation=np.zeros(level_length)
                represenation[sub_index]=1.
            return represenation

        # Transform to vector
        q00   = make_vector_for_index(0)
        #Dip
        q10   = make_vector_for_index(1)
        q11c  = make_vector_for_index(2)
        q11s  = make_vector_for_index(3)
        # Quad
        q20   = make_vector_for_index(4)
        q21c  = make_vector_for_index(5)
        q21s  = make_vector_for_index(6)
        q22c  = make_vector_for_index(7)
        q22s  = make_vector_for_index(8)
        q30   = make_vector_for_index(9)
        # Octo
        q31c  = make_vector_for_index(10)
        q31s  = make_vector_for_index(11)
        q32c  = make_vector_for_index(12)
        q32s  = make_vector_for_index(13)
        q33c  = make_vector_for_index(14)
        q33s  = make_vector_for_index(15)
        # Hexa
        q40   = make_vector_for_index(16)
        q41c  = make_vector_for_index(17)
        q41s  = make_vector_for_index(18)
        q42c  = make_vector_for_index(19)
        q42s  = make_vector_for_index(20)
        q43c  = make_vector_for_index(21)
        q43s  = make_vector_for_index(22)
        q44c  = make_vector_for_index(23)
        q44s  = make_vector_for_index(24)
        
        # Helper for better readability
        sq2  =np.sqrt(2.)
        sq3  =np.sqrt(3.)
        sq5  =np.sqrt(5.)
        sq8  =np.sqrt(8.)
        sq10 =np.sqrt(10.)
        sq12 =np.sqrt(12.)
        sq24 =np.sqrt(24.)
        sq35 =np.sqrt(35.)
        sq70 =np.sqrt(70.)

        # Dipoles
        x = q11c
        y = q11s
        z = q10
        # Quadrupoles
        xx = - 1/2   * q20   + sq3/2 * q22c
        yy = - 1/2   * q20   - sq3/2 * q22c
        zz =           q20
        xy =   sq3/2 * q22s
        xz =   sq3/2 * q21c
        yz =   sq3/2 * q21s
        # Octupoles
        xxx =   sq5/sq8  * q33c  -  sq3/sq8 * q31c
        xxy =   sq5/sq8  * q33s  -  1/sq24  * q31s
        xyy = - sq5/sq8  * q33c  -  1/sq24  * q31c
        yyy = - sq5/sq8  * q33s  -  sq3/sq8 * q31s
        xxz =   sq5/sq12 * q32c  -  1/2     * q30
        xyz =   sq5/sq12 * q32s
        yyz = - sq5/sq12 * q32c  -  1/2     * q30
        xzz =   sq2/sq3  * q31c
        yzz =   sq2/sq3  * q31s
        zzz =              q30
        ## # hexadecapoles
        xxxx =    3/8   * q40  -  sq5/4  * q42c + sq35/8*q44c
        xxxy =    1/8   * (  -sq5    * q42s  +  sq35 * q44s  )
        xxyy =    1/8   * q40   - sq35/8 * q44c
        xyyy =  - 1/8   * (  sq5     * q42s  +  sq35 * q44s  )
        yyyy =    3/8   * q40   + sq5/4*q42c + sq35/8*q44c
        xxxz =    1/16  * (  -3*sq10 * q41c  +  sq70 * q43c  )
        xxyz =    1/16  * (  -sq10   * q41s  +  sq70 * q43s  )
        xyyz =  - 1/16  * (  sq10    * q41c  +  sq70 * q43c  )
        yyyz =  - 1/16  * (  3*sq10  * q41s  +  sq70 * q43s  )
        xxzz =  - 1/2   * q40   +  sq5/4  * q42c
        xyzz =  sq5/4   * q42s
        yyzz =  - 1/2   * q40   -  sq5/4  * q42c
        xzzz =  sq5/sq8 * q41c
        yzzz =  sq5/sq8 * q41s
        zzzz =            q40

        monos=np.array([
            q00
        ])
        dipos=np.array([
            x, y, z,
        ])
        quads=np.array([ 
        xx, yy, zz, xy, xz, yz, 
        ])
        octos=np.array([
            xxx, xxy, xyy, yyy, xxz, 
            xyz, yyz, xzz, yzz, zzz,
        ])
        hexas=np.array([
            xxxx, xxxy, xxyy, xyyy, yyyy, 
            xxxz, xxyz, xyyz, yyyz, xxzz,
            xyzz, yyzz, xzzz, yzzz, zzzz,
        ])

        return [
            monos, dipos, quads, octos, hexas
        ]


class multipoles_data(BaseModel):
    atom_types: List[str|int]
    atom_coordinates: mynpar
    site_dmp_di: List
    tag_to_rank: dict
    output_file: str = None
    input_file: str
    _au_to_kj_mol=2625.49963947
    class Config:
        arbitrary_types_allowed=True
    @property
    def ranks_list(self):
        if type(self.site_dmp_di)==type(None):
            raise Exception(f"Asked for rank, but moments not initialized")
        tag_to_rank=self.tag_to_rank
        ranks_list=[]
        for mom in self.site_dmp_di:
            tags=mom.keys()
            if any( [t for t in tags if t not in tag_to_rank.values() ]):
                raise Exception(f"Moment features tags not known:\n{tags} {tag_to_rank}")
            ranks=[]
            for t in tags:
                rank=[ k for k in tag_to_rank.keys() if t==tag_to_rank[k] ]
                if len(rank)!=1:
                    raise Exception(f"Did not found rank for key {t} in {tag_to_rank}")
                ranks.append(rank[0])
            ranks_list.append(ranks)
        return ranks_list
    @property
    def max_rank(self):
        """ Return the maximum rank that occurs in all moments """
        finished=False
        max_rank=-1
        ranks_list=self.ranks_list
        for i in range(100):
            if sum([ 1 for x in ranks_list if not i in x ])>0:
                break
            else:
                max_rank=i
        if max_rank<0:
            raise Exception(f"Maximum rank smaller than 0, something is off")
        return max_rank

class multipoles(multipoles_data):
    class Config:
        arbitrary_types_allowed=True
    def flatten_moment(self,mom):
        if isinstance(mom, dict):
            new_mom=[]
            for k in mom.keys():
                new_mom+=[ float(x) for x in list(mom[k]) ]
        else:
            raise Exception(f"Unkown type")
        return new_mom
    #
    def __init__(self,inp_fi):
        if type(inp_fi)!=str:
            raise Exception(f"Expected input as string, got {type(inp_fi)} ({inp_fi})")
        else:
            if inp_fi.endswith('.mom'):
                input_file=inp_fi
                multipoles_raw=read_mom_file(inp_fi)
                atom_types=multipoles_raw.atom_types
                atom_coordinates=multipoles_raw.atom_coordinates
                site_dmp_di=multipoles_raw.site_dmp_di
                tag_to_rank=multipoles_raw.tag_to_rank
            else:
                raise Exception(f"Error in creating multipoles")
        super().__init__(input_file=input_file, atom_types=atom_types, atom_coordinates=atom_coordinates, site_dmp_di=site_dmp_di, tag_to_rank=tag_to_rank)
    def return_geom(self, connectivity=False):
        geom=geometry(atom_types=self.atom_types, coordinates=self.atom_coordinates)
        geom.connec()
        if connectivity:
            geom.connec()
        return geom
    @validate_call
    def which_ranks_evaluate(self,lmax: int | None = None,ranks =None ) -> List[int]:
        """ Determine list of ranks that will be evaluated """
        do_ranks=[]

        assert [ isinstance(x, type(None)) for x in [lmax,ranks] ].count(True) <2, f"Expect exatecly one argument others are NONE, got: lmax={lmax}, ranks={ranks}"
        # If nothing provided go with the following
        if lmax is None and ranks is None:
            do_ranks=self.ranks_list
        if ranks is not None:
            if isinstance(ranks, str):
                if bool(re.fullmatch(r"max[0-9]*",ranks)):
                    max_rank=int(ranks.replace('max',''))
                    do_ranks=[ [x for x in y if x<=max_rank] for y in self.ranks_list
                              ]
                elif ranks in [ 'all' , 'all_available' ]:
                    do_ranks=copy.deepcopy(self.ranks_list)
                else: raise Exception(f"Provided key ranks=\'{ranks}\', which is unkown")
            elif isinstance(ranks, (list, np.ndarray)):
                do_ranks=[]
                for the_set in self.rank_list:
                    if not all([x in the_set for x in ranks]):
                        print(f"Warning: requested ranks ({ranks}) are not supported by available ranks ({the_set})")
                    do_ranks.append([x for x in the_set if x in ranks])
            else:
                raise Exception(f"Unkown type for ranks: {ranks}")

        # 
        from itertools import chain
        all_ranks=sorted(list(set( list(chain(*do_ranks)) )))
        ranks_to_sites={}
        for rank in all_ranks:
            ranks_to_sites.update({rank: [k for k,v in enumerate(do_ranks) if rank in v]})

        return ranks_to_sites, do_ranks
    @validate_call
    def eval_on_grid(self, grid=None, onam=None, laf=None, lmax=None, ranks =None, change_at_cor=True, print_ranks=True):
        """
        1. Create grid from provided grid (either file or grid object)
        2. In case a global cartesian representation does not exist, get its transformation
            3. Generate an ESP map
        """
        # Check the grid
        if type(grid)==type(None):
            raise Exception(f"You have to provide a grid!")
        elif isinstance(grid, molecular_isodensity_grid):
            pass
        elif type(grid)==str:
            if os.path.isfile(grid):
                grid=molecular_isodensity_grid(grid)
            else:
                raise Exception(f"Passed string but is not valid file {grid}")
        else:
            raise Exception(f" Cannot handle type of grid :\n {type(grid)}")
        
        ranks_to_sites, ranks_to_do=self.which_ranks_evaluate(lmax=lmax, ranks=ranks)

        eval_points=grid.vertices
        
        mol_esp={} # list of every atomic contribution to grid
        for r, site_mask in ranks_to_sites.items():
            rk_tag=self.tag_to_rank[r]
            at_map=[]
            for coor, mom in [ (self.atom_coordinates[i], self.site_dmp_di[i]) for i in site_mask]:
                if rk_tag in mom.keys():
                    mom_cmp=mom[rk_tag]
                    vectors=eval_points-coor
                    T_mat=self.get_T_matrix(vectors,rank=r)
                    if mom_cmp.shape!=T_mat[:,0].shape:
                        raise Exception(F"Uneqal first dimension between moments and T matrix")

                    new_scalars=self._au_to_kj_mol*np.matmul(mom_cmp,T_mat)
                    #new_map=grid.provide_scalars( new_scalars)
                    at_map.append(new_scalars)
                else:
                    raise Exception(f"Rank \'{rk_tag}\' not in {mom}")
            total_scalars=sum(at_map)
            total_map=grid.provide_scalars(total_scalars, energy_units='kj/mol')
            total_map.gen_statistical_indicators()
            mol_esp.update({rk_tag : total_map})
        total_esp=np.sum(np.array([x.scalars for x in mol_esp.values()]), axis=0)
        total_map=grid.provide_scalars(total_esp, energy_units='kj/mol')
        total_map.gen_statistical_indicators()
        return total_map, mol_esp
    def get_T_matrix(self, vectors, rank=None,  plus_facs=True):
        f"""
        Produce T matrix in spherical tensor representation according to Stone's convention.
        See Appendix E, p. 285 in Stone's Intermolecular Interaction book.
        Then multiplies by  r^(-1-k) * (Q^k*E^k)*r^(-k)      
        ( factorial and prefactor are accounted for in spherical harmonics expansion)
        """
        x,y,z=[ vectors[:,i] for i in range(3) ]
        r=np.sqrt( x**2 + y**2 + z**2 )
        if type(rank)==type(None):
            raise Exception(f"Cannot proces rank {rank}")
        elif rank in [0,'mon']:
            T_matrix= np.array( [ np.ones(x.size)  ])
        elif rank==1:
            T_matrix=np.array( [ z, x, y])
        elif rank==2:
            q20  = 1./2.         * (3.*z**2 - r**2)
            q21c = np.sqrt(3.)    * x*z
            q21s = np.sqrt(3.)    * y*z
            q22c = np.sqrt(3.)/2. * ( x**2 - y**2)
            q22s = np.sqrt(3.)    * x*y
            T_matrix=np.array( [ q20,q21c,q21s,q22c,q22s ])
        elif rank==3:
            q30     = 1./2.                 * (5.*z**3 -3.*z*r**2)
            q31c    = np.sqrt(6.)/4.        * (5.*x*z*z - x*r**2)
            q31s    = np.sqrt(6.)/4.        * (5.*y*z*z - y*r**2)
            q32c    = np.sqrt(15)/2.        * ( x*x*z-y*y*z ) 
            q32s    = np.sqrt(15)           * x*y*z
            q33c    = np.sqrt(10.)/4.       * (x**3 - 3.*x*y*y)
            q33s    = np.sqrt(10.)/4.       * (3.*x*x*y - y**3)
            T_matrix=np.array( [q30,q31c, q31s, q32c, q32s, q33c, q33s ])
        elif rank==4:
            q40     = 1./8.                 * (35.*z**4 - 30.*z**2*r**2 +3.*r**4 )
            q41c    = np.sqrt(10.)/4.       * (7.*x*z*z*z - 3.*x*z *r**2)
            q41s    = np.sqrt(10.)/4.       * (7.*y*z*z*z - 3.*y*z *r**2)
            q42c    = np.sqrt(5.)/4.        * (7.*z**2 - r**2) * (x**2 - y**2 )
            q42s    = np.sqrt(5.)/2.        * (7.*z**2 - r**2) * ( x*y )
            q43c    = np.sqrt(70.)/4.       * x*z *  (x**2 - 3.*y**2)
            q43s    = np.sqrt(70.)/4.       * y*z *  (3.*x**2 - y**2)
            q44c    = np.sqrt(35.)/8.       * (x**4 - 6.*x*x*y*y + y**4)
            q44s    = np.sqrt(35.)/2.       * (x*x*x*y - x*y*y*y)
            T_matrix=np.array( [ q40 , q41c, q41s, q42c, q42s, q43c, q43s, q44c, q44s, ])
        else:    
            raise Exception(f"Rank \'{rank}\' not yet implemented")

        if plus_facs:
            # r^(-1-2*k)
            div=np.sqrt(np.array([ x**2+y**2+z**2 ]))**(-1-2*rank)
            T_matrix*=div #(-1.)**(rank) * div
            
        return T_matrix
    def to_dict(self):
        moments_per_site=self.site_dmp_di
        for i,m in enumerate(moments_per_site):
            for k in m.keys():
                moments_per_site[i][k]= moments_per_site[i][k].tolist()
        moment_type='spherical'
        convention='stone'
        atom_types=self.atom_types
        coordinates=np.array(self.atom_coordinates).tolist()
        ranks=self.ranks_list

        return {
                'moments':[ self.flatten_moment(mom) for mom in moments_per_site],
            'moment_type':moment_type,
            'convention':convention,
            'ranks':ranks,
            'atom_types':atom_types,
            'coordinates': coordinates,
        }
    @property
    def get_moments(self):
        convention='cart'
        if isinstance(convention, type(None)): raise Exception(f"Provide convention to this function")
        elif convention=='cart':
            new_array=[]
            for ar in self.site_dmp_di:
                ar=self.flatten_moment(ar)
                new_array.append(rotate_multipoles().rotate_quick(ar))
            return new_array

# Interaction maps
from . import *

from qc_objects.util.util import load_json

class molecular_isodensity_grid_meta(BaseModel):
    num_vertices: int=-1
    num_faces: int=-1
    iso_density_value: float

class molecular_isodensity_grid_data(BaseModel):
    class Config:
        arbitrary_types_allowed=True
    input_file: str | None = None
    output_file: str | None = None
    vertices: mynpar
    length_units: str
    faces: mynpar
    meta: molecular_isodensity_grid_meta
    def __init__(self, *args, **kwargs):
        if not 'iso_density_value' in kwargs.keys(): raise Exception(f"{kwargs.keys()}")
        kwargs.update({'meta':molecular_isodensity_grid_meta(*args,**kwargs)})
        super().__init__(*args, **kwargs)
        self.meta.num_vertices=len(self.vertices)
        self.meta.num_faces=len(self.faces)
    # What is this
    #sites: List | None=None
    #dic={}
    #double_lines_di: dict | None=None

class molecular_isodensity_grid(molecular_isodensity_grid_data):
    def __init__(self,inp=None,  mode='grid', **kwargs):
        if type(inp)==str:
            read_file=False
            if inp=='empty':
                info={}
            elif type(inp)==str:
                assert os.path.isfile(inp), f"Input {inp} is not a file"
                info=dict(input_file=inp)
                assert inp.split('.')[-1] in ['map','grid','json'], f"No valid file extension {inp}"
                read_file=True
            else:
                raise Exception()
            if read_file:
                info.update(self.read_file(inp))
        elif type(inp)==dict:
            info=self.process_dict(inp)
        elif inp is None:
            info={}
        else:
            raise Exception(f"Unkown type of input: {type(input)}: {input}")

        for key in list(set(info.keys()))+list(set(kwargs.keys())):
            if key in info.keys() and key in kwargs:
                if info[key]==kwargs[key] or any([x[key] is  None for x in [info, kwargs]]):
                    del kwargs[key]
                else:
                    raise Exception(f"Contradicting keys (for key={key}): kwargs[{key}]={kwargs[key]}, info[{key}]={info[key]}")
        super().__init__(**info, **kwargs)
    def read_file(self,inp):
        # Check extension and decide on mode
        units={}

        extension=inp.split('.')[-1]
        if extension in ['grid']:
            mode='grid'
        elif extension in ['map']:
            mode='map'
        elif extension in ['json']:
            mode='json'
        else:
            raise Exception(f"Do not know how to handle file extension {extension}: {inp}")
        
        if mode in ['grid','map']:
            info=read_map_file(inp, mode=mode)
            info=self.process_dict(info)
        elif mode in ['json']:
            try:
                info=load_json(inp)
                try:
                    info=self.process_dict(info)
                except Exception as ex:
                    raise Exception(f"Problem in processing data obtained from {inp}: {ex}")
            except Exception as ex:
                raise Exception(f"Problem when processing {inp} as json: {ex}")
        else:
            raise Exception(f"Did not implement mode {mode}")
        return info

    def process_dict(self,info):
        #assert 'data_type' in info.keys(), f"Excpectyed key \'data_type\' in {inp}"
        faces=info['faces']
        vertices=info['vertices']
        length_units = info['length_units']
        iso_density_value=info['meta']['iso_density_value']
        assert iso_density_value is not None, f"Isodensity value is None but needs to be definde: (check file: \'{inp}\')"
        
        return dict(
            vertices=vertices,
            length_units=length_units,
            iso_density_value=iso_density_value,
            faces=faces,
        )

    def read_dat_file(self,inp):
        with open(inp,'r') as rd:
            lines=[l.strip() for l in rd.readlines()]
        
        # get the index of the esp
        lines_no=[ (i,l) for i,l in enumerate(lines) ]
        label_line=[ (i,l) for i,l in lines_no if l.startswith("LABELS") ]
        if len(label_line)!=1:
            raise Exception(f"Did not find label line in {inp}")
        label_no, label_line=label_line[0]
        esp_ind=[ i-1 for i,l in enumerate(label_line.split()) if l=='E(1)elst' ]
        if len(esp_ind)!=1:
            raise Exception()
        esp_ind=esp_ind[0]

        # get scalars
        end_line=[ (i,l) for i,l in lines_no if l=="END" ]
        if len(end_line)!=1:
            raise Exception()
        end_no, end_line=end_line[0]
        try:
            esp_lines=lines[label_no+1:end_no]
            esp_vals=[ line.split()[esp_ind] for line in esp_lines ]
            esp_vals=[ float(x) for x in esp_vals ]
            esp_vals=np.array(esp_vals)
        except Exception as e:
            raise Exception(f"{e}")
        
        # Check units
        end_line=[ (i,l) for i,l in lines_no if l.startswith("ENERGY-UNITS") ]
        if len(end_line)!=1:
            raise Exception()
        eng_unit=end_line[0][1].split()[1]
        if not 'ENERGY' in self.units.keys():
            self.units.update({'ENERGY':'kj/mol'})
        if self.units['ENERGY']!='kj/mol':
            raise Exception()
        if eng_unit in ['CM-1']:
            conv=1.196266e-2
        else:
            raise Exception()
        esp_vals= conv* esp_vals 

        return esp_vals
    
    @validate_call( config=ConfigDict(arbitrary_types_allowed=True) )
    def provide_scalars(self,scalars:mynpar, energy_units: str):
        new_map=interaction_map(inp=self , scalars=scalars, energy_units=energy_units)
        return new_map

class statistical_indicators(BaseModel):
    avr: float 
    min: float 
    max: float
    rms: float
    abs: float
    min_rel_avr: float 
    max_rel_avr: float
    rms_rel_avr: float
    abs_rel_avr: float
    def __init__(self,scalars):
        avr=np.average(scalars)
        min=np.min(scalars)
        max=np.max(scalars)
        min_rel_avr=min-avr
        max_rel_avr=max-avr
        def get_rmse(vals):
            avr=np.average(vals)
            return np.sqrt( np.average ( (vals-avr)**2 ) )
        def get_abs(vals):
            return np.average(np.abs(vals))
        rms     = get_rmse(scalars)
        abs     = get_abs(scalars)
        rms_rel_avr     = get_rmse(scalars-avr)
        abs_rel_avr     = get_abs(scalars-avr)

        super().__init__(
            avr=avr,min=min, max=max, rms=rms, abs=abs, 
            min_rel_avr=min_rel_avr, max_rel_avr=max_rel_avr, rms_rel_avr=rms_rel_avr, abs_rel_avr=abs_rel_avr
        )

        #pos_dev = [ x for x  in self.scalars if x>=average ]
        #neg_dev = [ x for x in self.scalars if x<=average ]
        #rmse_pos=get_rmse(pos_dev)
        #rmse_neg=get_rmse(neg_dev)
class interaction_map_data(molecular_isodensity_grid):
    scalars: mynpar 
    statistics: statistical_indicators | None = None
    charge_violation: float|None=None
    energy_units: str
    def __init__(self,*args,  **kwargs):
        try:
            super().__init__(*args, **kwargs)
            if self.scalars.shape != self.vertices[:,0].shape:
                raise Exception(f"Shapes of scalars and grid do not match : {self.scalars.shape}, {self.vertices.shape}")
        except ValidationError as val_ex:
            raise Exception(f"Pydantic error in {val_ex.title}:\nerrors={val_ex.errors()}")
class interaction_map(interaction_map_data):
        #if isinstance(self, interaction_map):
        #    if type(scalars)==type(None):
        #        scalars=info_dict['values']
        #    if type(scalars)==type(np.a_rray([])):
        #        scalars=scalars
        #    elif type(scalars)==str:
        #        scalars=read_dat_file(scalars)
        #        units.update({'ENERGY':'kj/mol'})
        #    if scalars.shape[0]!=grid_points.shape[0]:
        #        raise Exception(f"Sizes dont match: {scalars.shape[0]} {grid_points.shape[0]}\n{inp}")
        #    dic.update(dict(scalars=scalars))
    def __init__(self,inp, scalars:mynpar|None=None, energy_units:str|None=None):#, charge_violation=None):#length_units:str|None=None, charge_violation:int=None):
        """ If instance of class has been provided just load the old class,
        other type should be packed in the correct dictinoary type"""

        # Get things in dictionary shape
        if isinstance(inp, (interaction_map, molecular_isodensity_grid) ):
            inp=inp.model_dump()
        elif isinstance(inp, dict):
            pass
        elif isinstance(inp, str):
            inp=load_json_or_yaml(inp)
        else:
            raise Exception(f"Unkown data type for input: {type(inp)}")

        # Process the dictionary    
        scalars_key='scalars'
        energy_units_key='energy_units'
        if scalars_key in inp.keys():
            assert scalars is None, f"Double definition of {scalars_key}"
            scalars=inp[scalars_key]
            del inp[scalars_key]
        if energy_units_key in inp.keys():
            assert energy_units is None, f"Double definition of {energy_units_key}"
            energy_units=inp[energy_units_key]
            del inp[energy_units_key]
        
        # Initialize class
        super().__init__(inp=inp, scalars=scalars, energy_units=energy_units)
    def gen_statistical_indicators(self):
        self.statistics=statistical_indicators(self.scalars)
        return self.statistics.model_dump()
    @validate_call
    def to_file(self, onam: str|None, map_dir: bool=True, format:str='json'):
        """
        Prints current grid to esp-map
        """
        # defaults for output name
        if isinstance(onam, type(None)):
            if not isinstance(self.output_file, type(None)):
                onam=self.output_file
            else:
                onam='dummy'
        if format.lower() in ['json']:
            onam=onam.replace('.json','')+'.json'
            dic=self.model_dump()
            with open(onam,'w') as wr:
                json.dump(dic, wr)
        else:
            
            def make_line(xyz,eng):
                line_ls    = [ f"{x:.8f}" for x in xyz ]
                line_ls    += [ f"{eng:.6f}" ]
                return line_ls
            di=self.double_lines_di
            print_keys=['mol_nam','dens_type', 'iso-dens', 'units', 'sites', 'step_size', 'num_points','num_triangles']
            avail_keys=list(di.keys())

            lines=[]
            for p_k in print_keys:
                line=''
                if p_k in ['sites','units']:
                    if p_k in ['sites']:
                        line='\n'.join(self.sites)
                    if p_k in ['units']:
                        line="{:<15}  ".format("UNITS")+' '.join(list(self.units.values() ))
                elif p_k in avail_keys:
                    line+= '{:<15}  '.format(di[p_k]['str'])
                    if not p_k in self.dic:
                        raise Exception()
                    val=self.dic[p_k]
                    if type(val) in [int, str, float]:
                        line+=f"{val}"
                    elif type(val) in [ list ]:
                        line+=f"{' '.join(val)}"
                    elif type(val) in [dict ]:
                        line+=f"{' '.join( list(val.values() ) )}"
                    else:
                        raise Exception(f"{type(val)}")
                else:
                    raise Exception(f"Abort: {p_k}")
                if line!='':
                    lines.append(line+'\n')

            
            # Print data
            lines.append('BEGIN DATA\n')
            for i, coor in enumerate(self.grid_points):
                scalar=self.scalars[i]
                lines.append("   {:>12} {:>12} {:>12} {:>12}\n".format( *[f"{x:.8f}" for x in coor], f"{scalar:.8f}") )
            for i, triang in enumerate(self.triangles): # Dont forget to follow index convention
                lines.append("  {:>5} {:>5} {:>5}\n".format( *[f"{int(x)+1}" for x in triang ]) ) 
            lines.append('END DATA\n')

            # Write
            if type(onam)==str:
                pass
            elif type(onam)==type(None):
                onam='test.map'
            elif type(onam)==dict:
                if 'rank' in onam.keys():
                    rank=f"_{onam['rank']}"
                else:
                    rank=''
                if 'lead' in onam.keys():
                    lead=f"{onam['lead']}"
                else:
                    lead='ESP-MAP'

                if 'iso-dens' in onam.keys():
                    iso=f"_iso-{onam['iso-dens']}".strip('0').replace('.','p')
                elif type(self.iso_density_value)!=type(None):
                    iso='_iso-'+f"{self.iso_density_value}".strip('0').replace('.','p')
                else:
                    iso=''
                onam=f"{lead}{rank}{iso}.map"
            else:
                raise Exception()
            if map_dir:
                m_io.mkdir("maps", mode="ifnotexists")
                if os.path.dirname(onam)=="": 
                    onam="maps/"+os.path.basename(onam)
                else:
                    onam=os.path.dirname(onam)+"/maps/"+os.path.basename(onam)
            with open(onam, 'w') as wr:
                wr.writelines(lines)
        self.output_file=os.path.realpath(onam)
        print(f"Printed map to {onam}")
        return onam
# Subtraction of maps
class interaction_map(interaction_map):
    @validate_call
    def add_maps(self, ref_map:interaction_map, fac=1., fac_self=1.):

        def check_map_compability(map_0,map_1):
            # check grid points are the same
            assert map_0.vertices.shape==map_1.vertices.shape, f"Coordinates of maps do not have same shape"
            assert map_0.length_units == map_1.length_units, f"Length units between maps are different, but could implement conversion (TODO):\n{map_0.length_units}\n{map_1.length_units}"
            for trial_point, ref_point in zip(map_0.vertices, map_1.vertices):
                thres=1.e-8
                assert np.linalg.norm( trial_point-ref_point)<thres, f"Difference between points (L2 norm) larger than {thres}:\n{trial_point.tolist()}\n{ref_point.tolist()}"
            # assert energy units are the same
            assert map_0.energy_units.lower()==map_1.energy_units.lower(), f"Energy units between maps are different, but could implement convertsion (TODO):\n{map_0.energy_units}\n{map_1.energy_units}"
            return True

        # Check that both maps are compatible
        check_map_compability(self, ref_map)

        # Calculate the new scalars
        new_scalars=self.scalars*fac_self + fac*ref_map.scalars

        # Create new map with updated scalars
        cur_map_dic=self.model_dump()
        del cur_map_dic['scalars']
        new_map=interaction_map(inp=cur_map_dic, scalars=new_scalars)
        new_map.gen_statistical_indicators()

        return new_map
    @validate_call
    def sub(self, ref_map:interaction_map):
        return self.add_maps(ref_map, fac=-1)

class grid():
    input_file=None
    grid_points=None
    values=None
    triangles=None
    statistical_indicators=None
    output_file=None
    def __init__(self,inp=None):
        if type(None)==type(inp): # Start with empty object
            pass
        elif str==type(inp): # Read a file
            # Read file according to extension
            self.input_file=inp
            if inp.endswith('.map'):
                info_dict=read_map_file(inp).dic
                self.grid_points=info_dict['points']
                self.values=info_dict['values']
                self.triangles=info_dict['triangles']
            if inp.endswith('.grid'):
                info_dict=read_grid_file.dic
            else:
                raise Exception(f"Unkown extension of file \"{inp}\"")
    def gen_statistical_indicators(self):
        def get_rmse(vals):
            avr=np.average(vals)
            return np.sqrt( np.average ( (vals-avr)**2 ) )
        average = np.average(self.values)
        minimum     = np.min(self.values)
        maximum     = np.max(self.values)
        min_dev = minimum-average
        max_dev = maximum-average
        rmse    = get_rmse(self.values)
        pos_dev = [ x for x  in self.values if x>=average ]
        neg_dev = [ x for x in self.values if x<=average ]
        rmse_pos=get_rmse(pos_dev)
        rmse_neg=get_rmse(neg_dev)
        self.statistical_indicators={
            'average':average, 'minimum':minimum, 'maximum':maximum, 
            'min_dev':min_dev, 'max_dev':max_dev, 'rmse':rmse
        }

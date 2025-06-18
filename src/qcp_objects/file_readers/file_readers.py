from . import *

class file_reader():
    lines=None
    lines_no=None
    comtag=None
    input_file=None
    dic={}
    # Init function checks if file exists reads in the line from file
    def __init__(self, inp_fi):
        from qc_global_utilities.file_handling import temporary_file

        self.input_file=inp_fi
        if self.input_file.endswith('.xz'):
            import lzma
            lines=lzma.open(self.input_file, encoding='utf-8', mode='rt').read()
        else:
            with open(self.input_file, encoding='utf-8-sig', errors='ignore') as rd:
                lines=rd.read()
        
        self.lines=[ l.rstrip('\n') for l in lines.split('\n') ]
        if len(self.lines)<1:
            raise Exception(f"Empty File: \"{self.input_file}\"")
    @property
    def lines_no(self):
        if type(self.lines)==type(None):
            raise Exception(f"Lines are not initialized")
        if type(self.lines)!=list:
            raise Exception(f"Lines are not in format list but {type(self.lines)}")
        if len(self.lines)<1:
            raise Exception(f"Lines is a list with no entries")

        self._lines_no=[ (i,l) for i,l in enumerate(self.lines) ]
        return self._lines_no
    @property
    def lines_nocom(self):
        if type(self.lines)==type(None):
            raise Exception(f"Lines are not initialized")
        if type(self.lines)!=list:
            raise Exception(f"Lines are not in format list but {type(self.lines)}")
        if len(self.lines)<1:
            raise Exception(f"Lines is a list with no entries")

        if type(self.comtag)==type(None):
            raise Exception(f"Comment tag is not defined, define as \' self.comtag \'"+
                    "\n in case you do not want to declare com tag provide empty list \'[]\'")
        else:
            self._lines_nocom=self.lines
            for comtag in self.comtag:
                self._lines_nocom=[ l.split(comtag)[0] for i,l in enumerate(self._lines_nocom) ]

        return self._lines_nocom
    @property
    def lines_num_nocom(self):
        return [ (i,l) for i,l in enumerate(self.lines_nocom) ]
    def get_line_no(self, match=None, fullmatch=False, single=False, in_range=None):
        # Check that type is a string
        if type(match)!=None:
            if type(match)!=str:
                raise Exception()
        else:
            raise Exception()
        if not fullmatch:
            results=[ v[0] for v in self.lines_num_nocom if re.search( match, v[1] ) ]
        else:
            results=[ v[0] for v in self.lines_num_nocom if re.fullmatch( match, v[1] ) ]
        # Get range
        if type(in_range)==type(None):
            pass
        else:
            results=[ x for x in results if x>=in_range[0] and x<=in_range[1]]

        # Check number
        if single:
            if len(results)!=1:
                raise Exception(self.input_file, results, match)
            return results[0]
        else:
            return results
    def get_line_item(self,l, tot_leng=None, get_item=None, data_type=None, sep=None):
        if type(l)==str:
            pass
        elif type(l)==int:
            if l>=len(self.lines) or l<0:
                raise Exception(f"Requested line number \'{l}\' is not in allowed range ([0,{len(self.lines)}]).")
            else:
                l=self.lines_nocom[l]
        else:
            raise Exception()

        if type(sep)==type(None):
            ar=l.split()
        else:
            raise Exception(f"Not implemented")

        if type(tot_leng)!=type(None):
            if type(tot_leng)!=int:
                raise Exception(f"Wrong type for option \'tot_leng\':\n {tot_leng} ({type(tot_leng)})")
            else:
                if len(ar)!=tot_leng:
                    raise Exception(f"Expected line of length {tot_leng} but got {len(ar)}:\n{l}, {ar}")

        return_single=False
        if type(get_item)==int:
            get_item=[get_item]
            return_single=True
        elif type(get_item)==str:
            if get_item.endswith('+'):
                try:
                    get_item=int(get_item[:-1])
                except Exception as e:
                    raise Exception(f"{e}")
                get_item=[ i for i in range(len(ar)) if i>=get_item]
            else:
                raise Exception()
        elif type(get_item)!=list:
            raise Exception(f"Type of number of item to get is wrong, should be integer:\n{get_item} ({type(get_item)})")

        if any( [x>len(ar) for x in get_item ] ):
            raise Exception(f"Number of get item on this line is {get_item} longer than elements on line:\n{l}")
        if any( [x<0 for x in get_item] ):
            raise Exception()
        else:
            vals=[ ar[x] for x in get_item ]

        if type(data_type)!=type(None):
            try:
                vals=[ data_type(x) for x in vals ]
            except Exception as e:
                raise Exception(f"Could not convert value {vals} for type conversion {data_type}")
        if return_single:
            vals=vals[0]
        return vals
    def read_data(self,start=None, end=None, data_type=None, object_type=None, line_items=None):
        if start<0 or start>=len(self.lines):
            raise Exception()
        if end<start+1 or end>len(self.lines):
            raise Exception()

        data=[]
        the_lines=self.lines_nocom[ start:end ]
        for l in the_lines:
            ar=l.split()

            #Check length
            if type(line_items)!=type(None):
                if type(line_items)!=int:
                    raise Exception()
                else:
                    if len(ar)!=line_items:
                        raise Exception(f"Length of line not as expects:\n{ar}\n{line_items}")
            
            if type(data_type)!=type(None):
                try:
                    ar=[ data_type(x) for x in ar ]
                except Exception as e:
                    raise Exception(f"Could not read {ar} with {data_type}:\n{e}")
            data.append(ar)
        if type(object_type)==type(None):
            return data
        elif object_type=='np':
            try:
                data=np.array(data)
            except Exception as e:
                raise Exception(f"Could not convert to numpy array: {e}")
            return data
        else:
            raise Exception()
        
class read_map_file(file_reader):
    comtag=['!']
    dic={}
    double_lines_di={'num_points':{'str':'POINTS','type':int},
            'step_size':{'str':'STEP', 'type':float},
            'num_triangles': {'str':'TRIANGLES', 'type':int},
            'iso-dens': {'str': 'ISO-DENSITY', 'type':float},
            'units':    {'str': 'UNITS', 'tot_leng':None, 'get_item':'1+', 'type':str},
            'dens_type': {'str':'DENSITY-TYPE', 'type':str},
            'mol_nam':      {'str':'MOL-NAME', 'type':str},
    }
    mode=None
    def __init__(self, inp_fi, mode='map'):
        file_reader.__init__(self,inp_fi)
        self.mode=mode
        self.extract_data(mode=mode)
    def extract_data(self, mode='map'):
        if not mode in ['map','grid']:
            raise Exception(f"Unkown keyword {mode}")

        for k in self.double_lines_di.keys():
            item=self.double_lines_di[k]
            key=item['str']
            format=item['type']
            lino=self.get_line_no(match=key, single=True)
            line=self.lines_nocom[lino]
            if not 'tot_leng' in item.keys():
                tot_leng=2
            elif type(item['tot_leng'])==type(None):
                tot_leng=item['tot_leng']
            else:
                raise Exception()
            if not 'get_item' in item.keys():
                get_item=1
            elif type(item['get_item'])!=type(None):
                get_item=item['get_item']
            else:
                raise Exception(f"Key for \'get_item\'")

            value=self.get_line_item( line, tot_leng=tot_leng, get_item=get_item, data_type=format)
            self.dic.update({k:value})

        no_pts=self.dic['num_points'] 
        # Get the data
        no_data_start=self.get_line_no(match='BEGIN DATA', single=True)
        no_data_end=self.get_line_no(match='END DATA', single=True)
        #-- Get points
        if mode in ['map']:
            line_item=4 # 3 cartesian components and a scalar
        else:
            line_item=3 # 3 carteisan components
        data=self.read_data(start=no_data_start+1, end=no_data_start+1+no_pts, data_type=float, object_type='np', 
                line_items=line_item)
        points=data[:,0:3]
        if mode in ['map']:
            values=data[:,3]
        else:
            values=None

        no_sites_start  = self.get_line_no(match='BEGIN SITES', single=True)
        no_sites_end    = self.get_line_no(match='END SITES', single=True)
        sites=self.lines_nocom[no_sites_start:no_sites_end+1]

        #-- Get the triangles
        no_triang=self.dic['num_triangles']
        triangles=self.read_data(start=no_data_start+1+no_pts, end=no_data_start+1+no_pts+no_triang, 
                data_type=float, object_type='np', line_items=3)
        triangles=triangles-1 # python index convention

        self.dic.update({'units':self.process_units(self.dic['units'])})
        self.dic.update({
            'points': points,
            'values': values,
            'triangles': triangles,
            'sites': sites
        })
    def process_units(self,list_units):
        unit_map={ 'LENGTH': ['BOHR','ANGSTROM'],
                    'ENERGY': ['AU', 'kj/mol'],
        }
        all_found=[]
        ret_di={}
        for k in unit_map.keys():
            if k=='ENERGY' and self.mode=='grid': # Dont read energy unit for pure grid, there are no scalars
                continue
            found=[ u for u in list_units if u in unit_map[k] ]
            if len(found)==0:
                warn(f"No unit provided for {k}")
            elif len(found)>1:
                m_io.out().error(f"Found more than one key corresponding to units of \'{k}\': \n{found}")
            else:
                ret_di.update({k:found[0]})
            all_found+=found
        missing=[ x for x in list_units if x not in all_found ]
        if len(missing)>0:
            raise Exception(f"Could not interprete unit keys:\n{missing}")
        return ret_di


class json_reader():
    def __init__(self, input):
        
        # Load json content
        assert isinstance(input, str)
        assert os.path.isfile(input)
        with open(input, 'r') as rd:
            content=json.load(rd)

        # Check format
        assert isinstance(content, dict)
        mand_keys=dict(
            data_type=['data_type']
        )
        key_di=dict(
            vertices=['vertices'],
            faces=['faces']
        )
    def check_keys(dic, keys):
        found= [ k for k in keys if k  in dic.keys() ]
        if len(found)==0:
            found=False
            value=None
        elif len(found)==1:
            found=True
            value=dic[k]
        else:
            found=True
            value=[dic[k] for k in found]
            raise Exception(f"Found more than one key matching {keys}: {value}")
        return found, value

        
            

class fchk_data(BaseModel):
    mode            : str|None=None # mode e.g. SP or FORCE
    method          : str|None=None # electronic structure method used
    basis           : str           # basis set used
    atomic_numbers  : List[int] = []
    number_of_atoms : int = Field(default=-1, ge=1)
    coordinates     : List[List[float]]= []
    charge          : int 
    multiplicity    : int
    method          : str
    basis           : str

class fchk_reader(file_reader):
    comtag          : List[str]=[]
    fchk_data:fchk_data
    def __init__(self, inp_fi):
        # loads content of file and structures it into list
        super().__init__(inp_fi)
        self.fchk_data=self.extract_data()
    def extract_data(self):
        # Read the file

        # Numbered lines with truncated comments as process by partent
        lines_no=self.lines_nocom

        # The first line is a comment the second is the header
        def get_header(line):
            ar=line.split()
            assert len(ar)==3, f"Expected 3 entries on header line:\n   {line}"
            return ar
        header_line=lines_no[1] 
        mode, method, basis = get_header(header_line)

        # Get all lines with text
        text_line_search=r'^[A-Z]'
        line_idx_text= self.get_line_no(match=text_line_search)
        
        charge=self.get_line_item(self.get_line_no(match='Charge',single=True), tot_leng=3, data_type=int, get_item=2)
        multiplicity=self.get_line_item(self.get_line_no(match='Multiplicity',single=True), tot_leng=3, data_type=int, get_item=2)
        number_of_atoms=self.get_line_item(self.get_line_no(match='Number of atoms', single=True), tot_leng=5, get_item=4, data_type=int)
        
        block_ranges=[ [x,None] for x in line_idx_text ]
        for i,t in enumerate(line_idx_text):
            if i<len(line_idx_text)-1:
                next_num= line_idx_text[i+1]
            else:
                next_num=-1
            block_ranges[i][1]=next_num

        # Get the nuclear charges
        nuc_srch={'srch_str':'Atomic numbers', 'format':'int', 'target':"atomic_numbers"}
        coords_srch={'srch_str':'Current cartesian coordinates', 'format':'float_-1x3', 'target':"coordinates"}
        atomic_numbers= self.search_block(nuc_srch, block_ranges)
        coordinates=self.search_block(coords_srch, block_ranges)
        
        assert len(coordinates)>0
        assert len(atomic_numbers)>0
        return fchk_data(
            mode=mode, method=method, basis=basis,
            charge=charge, multiplicity=multiplicity, number_of_atoms=number_of_atoms, coordinates=coordinates, atomic_numbers=atomic_numbers
        )

    @validate_call
    def search_block(self,srch_di:dict, blocks:List[List[int]]):
        """ For Gaussian FCHK Files, find line that contain search string in provided list of lines
        then tries to convert to the right format (type and shape eg -1x3 for [:,3]) and saves in the class under traget
        string"""
        srch_str=srch_di['srch_str']
        format=srch_di['format']
        target=srch_di['target']

        the_block=[ (start,end)  for start,end in blocks if re.match(srch_str, self.lines[start])]
        assert len(the_block)==1, f"No block found for search: {srch_str} in file {self.input_file}"
        start,end=the_block[0]
        assert end>start+1, f"End and start are adjacent: {start}:{end} (srch_str={srch_str})"

        content=self.lines[start+1:end]

        # Read in in right format
        list_tmp=[]
        for l in content:
            tmp=l.strip().split()#['text'].strip().split()
            if len(format.split('_'))==2:
                format_type, dimension=format.split('_')
                dimension=tuple([int(x) for x in dimension.split('x')])
            else:
                format_type=format
                dimension=None
            conversion_di={'float':float, 'int':int,'str':str}
            if format_type not in conversion_di.keys():
                raise Exception(f"Conversion methods for requested format \'{format}\' not defined")
            try:
                tmp=[conversion_di[format_type](x) for x in tmp]
                list_tmp+=tmp
            except Exception as e:
                raise Exception(f"Could not read {format_type} line number {l}:\n{e}")
        
        if dimension!=None:
            list_tmp=np.array(list_tmp).reshape(dimension)
        return list_tmp


    @validate_call
    def make_geometry(self):
        dat=self.fchk_data
        return {
            'coordinates':dat.coordinates,
            'atom_types': [ nuclear_charge_to_element_symbol(x) for x in dat.atomic_numbers ],
            'multiplicity': dat.multiplicity,
            'charge': dat.charge,
        }



class geometry_reader(file_reader):
    coordinates=None
    atom_types=None
    num_sites=None
    charge=None
    multiplicity=None
    input_file=None
    def  __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.extract_data()
class read_xyz_file(geometry_reader):
    comtag=[]
    comments=None
    def  __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.extract_data()
    def extract_data(self):
        empty_lines=self.get_line_no(match=r'[ ]*', fullmatch=True)
        first_line=self.lines[0]
        try:
            site_number=int(first_line.strip())
        except Exception as e:
            raise Exception(f"Cannot convert to integer : \'{first_line.strip()}\'")

        # Comment
        comment_line = self.lines[1].strip()
        if comment_line != '':
            comments=comment_line
        
        at_ty_search=r'^\s*[A-Z][a-z,A-Z]? '
        at_lines=self.get_line_no(match=at_ty_search, fullmatch=False)
        if len(at_lines)==0:
            raise Exception(f"Could not read atom lines in {self.input_file}")
        if len(at_lines)!=site_number:
            raise Exception(f"XYZ file \'{self.input_file}\'should contain {site_number} of sites according to info, "+
                    f"but found {len(at_lines)}")
        
        atom_types=[]
        coordinates=[]
        for line in [ self.lines_nocom[i] for i in at_lines]:
            ar=line.split()
            if len(ar)!=4:
                raise Exception(f"Expected 4 element on line but found\n{line}")
            
            match=re.match(r'[A-Z][a-z]?', ar[0])
            if not bool(match):
                raise Exception(f"Could not finnd atom:\n{ar[0]}")
            at_ty=ar[0]
            if len(at_ty)==2:
                at_ty=at_ty[0].upper()+at_ty[1].lower()
            atom_types.append(at_ty)

            try:
                coor=ar[1:]
                coor=[ float(x) for x in coor ]
            except:
                raise Exception(f"{coor}")
            coordinates.append(coor)
        
        coordinates=np.array(coordinates)
        self.atom_types=atom_types
        self.coordinates=coordinates
    @validate_call
    def make_geometry(self):
        return {
            'coordinates':self.coordinates,
            'atom_types': self.atom_types,#[ nuclear_charge_to_element_symbol(x) for x in dat.atomic_numbers ],
            'multiplicity': None,# dat.multiplicity,
            'charge': None,#dat.charge,
        }
    def to_dict(self):
        keys=['coordinates','atom_types','num_sites','charge','multiplicity']
        return dict([ (k,getattr(self,k)) for k in keys])
class read_sdf_file(geometry_reader):
    comtag=[]
    comments=None
    conformation_id=None
    compound_id=None
    def  __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.extract_data()
    def extract_data(self):
        # Get number of sites
        no_num_sites=3 # where to find the site number
        self.num_sites=self.get_line_item(no_num_sites,get_item=0,data_type=int)
        # Read the next lines
        coordinates=[] ; atom_types=[]
        for i in range(no_num_sites+1, self.num_sites+no_num_sites+1):
            coor=self.get_line_item(i,get_item=[0,1,2],data_type=float)
            at_ty=self.get_line_item(i,get_item=3,data_type=str)
            coordinates.append(coor)
            atom_types.append(at_ty)
        self.coordinates=np.array(coordinates)
        self.atom_types=atom_types

        def get_pubchem_entry(string, number=None, type_conv=str, strict=False):
            empty_line_no=self.get_line_no(match='', fullmatch=True)
            target_line_no=self.get_line_no(match=string, fullmatch=True)
            if len(target_line_no)!=1:
                if strict:
                    raise Exception(f"Did not find only one line {string} for {self.input_file}:\n{target_line_no}")
                else:
                    return None
            else:
                target_line_no=target_line_no[0]
            limit=[ x for x in empty_line_no if x>target_line_no+1 ]
            if len(limit)==0:
                raise Exception
            else:
                limit=min(limit)
            target_lines=self.lines_nocom[target_line_no+1:limit]

            if type(number)==type(None):
                pass
            elif type(number)!=int:
                raise Exception()
            else:
                if len(target_lines)!=number:
                    raise Exception()

            try:
                conv_di={'int':int,
                        'str':str
                }
                target_values=[ conv_di[type_conv](x) for x in target_lines ]
                return target_values
            except Exception as e:
                raise Exception({e})

        self.molecule_name=self.lines_nocom[0].strip()

        ### Get the CID
        cid_str='> <PUBCHEM_COMPOUND_CID>'
        self.compound_id=get_pubchem_entry(cid_str, number=1, type_conv='int')[0]
        ### Get the conformation_ID
        conf_id_str='> <PUBCHEM_CONFORMER_ID>'
        self.conformation_id=get_pubchem_entry(conf_id_str, type_conv='str', strict=False)
        if type(self.conformation_id)==type(None):
            self.conformation_id=f"{self.compound_id}_1"
        elif len(self.conformation_id)==1:
            self.conformation_id=self.conformation_id[0]
            if not re.fullmatch(r'[0-9,A-Z]*', self.conformation_id):
                raise Exception(f"INVALID format of conformation id, expected numbers and capital characters:\n{self.conformation_id}")
        else:
            raise Exception()

        ### Check the charges
        # On a line of M CHG {num entrie} [ {at_index} {charge} ]
        charge_str=r'M[\ ]*CHG .*'
        charge_line_no=self.get_line_no(match=charge_str, fullmatch=True)
        if len(charge_line_no)==0:
            charge=0
        elif len(charge_line_no)==1:
            charge_line=self.lines_nocom[charge_line_no[0]]
            charge_element=[ int(x) for x  in charge_line.split()[2:] ]
            if len(charge_element)!=charge_element[0]*2+1:
                raise Exception(f"Wrong format {charge_line}")
            charge=0
            for i in charge_element[2::2]:
                charge+=i
            if charge!=0:
                raise Exception(f"Total Charge is not neutral")
        else:
            raise Exception()
        self.charge=charge
    def to_record(self):
        dic={
            'mol_name': self.molecule_name,
            'cid': self.compound_id,
            'conformation_id': self.conformation_id,
            'charge':self.charge,
            'coordinates':self.coordinates.tolist(),
            'atom_types':self.atom_types,
        }
        return dic
class read_psi4_output(file_reader):
    comtag=None
    def  __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.extract_data()
class read_psi4_meta(file_reader):
    comtag='#'
    def  __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.dic=self.read_meta_file()
    def read_meta_file(self):
        lines=[ l for l in self.lines_nocom if l.strip()!='' ]
        entries=[]
        for l in lines:
            ar = l.split(':')
            if len(ar)!=2:
                raise Exception(f"Not two elements separated by \':\'\n{line}")
            key, val = [x.strip() for x in ar ]
            entries+= [ (key,val) ]
        unique_keys=set( [ e[0] for e in entries ] )
        for key in list(unique_keys):
            duplicate= [e for e in entries if e[0].lower()==key.lower() ]
            if len(duplicate)!=1:
                raise Exception(f"Found duplicate keys for {key}: {duplicate}")
        return dict(entries)
class read_camcasp_solution(file_reader):
    comtag='!'
    def __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.dic=self.read_camcasp_solution_file()
    def read_camcasp_solution_file(self):
        ''' Read CAMCASP ISA solution FILE'''
        
        def get_block(tag, block_range=None, occurence=None):
            begin_str=rf"^[ ]*BEGIN {tag}"
            end_str=f"^[ ]*END {tag}"
            begin=self.get_line_no(begin_str, fullmatch=False)
            stop=self.get_line_no(end_str, fullmatch=False)
            if type(block_range)!=type(None):
                begin=[ x for x in begin if x>block_range[0] and x<block_range[1]]
                stop=[ x for x in stop if x>block_range[0] and x<block_range[1]]
            # Check integrity and create  start to stop pairs
            if len(begin)!=len(stop):
                raise Exception(f"Atom Begin and End statement unequal Frequency")
            whole_list, start_stop_pairs=[[] , [] ]
            for sta,end in zip(begin, stop):
                whole_list+=[sta,end]
                start_stop_pairs+= [ (sta,end) ]
            # Check if odering is good
            if any([ whole_list[i+1]-whole_list[i]<0 for i in range(len(whole_list)-1)  ]):
                raise Exception(f"The pairs are not well ordered {whole_list}")
            
            if type(occurence)==int:
                if not len(start_stop_pairs)==occurence:
                    raise Exception()
            elif type(occurence)==type(None):
                pass
            else:
                raise Exception()
            return start_stop_pairs

        def read_basis_functions(start, end):
            def get_lines(star,end):
                """ List of index lines (meta data) and list of list of values lines """
                """ Value lines can feature multiple pGTOs """
                valid_lines=[ x for x in self.lines_nocom[start+1:end] if x.strip()!='' ]
                info_lines=[]
                value_lines=[]
                for i,line in enumerate(valid_lines):
                    ar=line.split()
                    if len(ar)==6:
                        info_lines.append(line)
                        value_lines.append([])
                    elif len(ar)==3:
                        value_lines[-1].append( line )
                    else:
                        raise Exception(f"Excpected line to have 3 or 6 entries:{line}")
                assert len(info_lines)==len(value_lines), 'info and value line length not identical'
                return info_lines, value_lines
            def read_value_line(line):
                """ Reads a value line and return the values """
                v_ar=line.split()
                assert len(v_ar)==3, f"Value line should be length 6 but found {val}"
                index, normalization, exponent_scale =v_ar
                try:
                    normalization, exponent_scale = [ float(x) for x in [normalization, exponent_scale] ]
                except Exception as e:
                    raise Exception(f"Cannot convert line to to 3 number:\n{line}\n{e}")
                try:
                    index=int(index)
                except Exception as e:
                    raise Exception(f"Index (first item) cannot be converted to intger\n{info}\n{e}")
                return index, normalization, exponent_scale
            
            # Get a list of valid lines
            info_lines, value_lines= get_lines(start, end)

            # Parse the values of these lines
            indices, normalizations, exponent_scales, exponent_orders=[ [] for i in range(4)]
            for val_lines, info in zip(value_lines, info_lines):
                # Check the info line
                i_ar=info.split()
                assert len(i_ar)==6, f"Info line should be length 6 but found {val}\n INPUT FILE: {self.input_file}"
                index, site_index, symmetry,  l,  m, num_of_pGTO=i_ar
                assert symmetry=='S', f"Symmetry not of group S!\n{info}"
                num_of_pGTO=int(num_of_pGTO)
                
                # Get the values out of the exponent lines
                if len(val_lines)!=num_of_pGTO:
                    raise Exception(f"Excpected {num_of_pGTO} pGTOs")
                GTO_indices, GTO_exponent_scales, GTO_normalizations, GTO_exponent_orders=[ [] for i in range(4) ]
                for val_line in val_lines:
                    index, normalization, exponent_scale=read_value_line(val_line)
                    exponent_order=2 # since camcasp uses only gaussians atm
                    GTO_indices.append(index)
                    GTO_normalizations.append(normalization)
                    GTO_exponent_scales.append(exponent_scale)
                    GTO_exponent_orders.append(exponent_order)
                
                indices.append(GTO_indices)
                normalizations.append(GTO_normalizations)
                exponent_scales.append(GTO_exponent_scales)
                exponent_orders.append(GTO_exponent_orders)
            return {'normalizations':normalizations, 'exponent_scales':exponent_scales,
                    'exponent_orders':exponent_orders}

        def read_coefficients(start, end):
            valid_lines=[ x for x in self.lines_nocom[start+1:end] if x.strip()!='' ]
            coefficients,indices=[ [] for i in range(2) ]
            for line in valid_lines:
                ar=line.split()
                assert len(ar)==2, f"Not two elements on {line}"
                index, coefficient=ar
                try:
                    index=int(index)
                    coefficient=float(coefficient)
                except Exception as e:
                    raise Exception(f"Cannot read integer and float on line {line}\n{e}")
                coefficients.append(coefficient)
                indices.append(index)
            assert all([ i+1 in indices for i in range(len(valid_lines))])
            return {'coefficients':coefficients}
        
        # Get atom blocks and extract values
        atom_start_stop_pairs=get_block('ATOM')
        atoms=[]
        coordinates=[]
        atom_types=[]
        for start, stop in atom_start_stop_pairs:
            # Get Atom type -> search for TYPE, but this key also exist in ISA-SHAPEFUNCTION Block, hence dropt soltion that occur the
            tag='ISA-SHAPE-FUNCTION'
            shape_block=get_block(tag, block_range=[start,stop], occurence=1)[0]
            type_lino=self.get_line_no(r'^[ ]*TYPE', fullmatch=False, in_range=[start,stop])
            type_lino=[x for x in type_lino if x<shape_block[0] or x>shape_block[1] ]
            assert len(type_lino)==1, f"Did nof found one results for TYPE\n{type_lino}"
            type_lino=type_lino[0]
            # Read the line
            ar=self.lines[type_lino].split()
            assert len(ar)==2, f"Expected TYPE line with 2 items, found {self.lines[type_lino]}"
            atom_type=ar[1] 
            # Capitalize first letter, make rest small, assert that there are only 2 keys
            atom_type=atom_type.lower().capitalize()
            assert len(atom_type)<3, f"Raise Exception, atom type longer than 2 charcters found {self.lines[type_lino]}"
            atom_types.append(atom_type)

            # Get Coordinates
            tag='SITES'
            shape_block=get_block(tag, block_range=[start,stop], occurence=1)[0]
            assert shape_block[1]-shape_block[0]-1==1, f"Block {tag} should only feature one line, but does not"
            lines=self.lines_nocom[shape_block[0]+1:shape_block[1]]
            ar=lines[0].split()
            assert len(ar)==8
            index, tag, Z, x, y, z, type_tag, atom_type =ar
            try:
                x,y,z=[ float(a) for a  in [x,y,z]]
            except Exception as e:
                raise Exception(f"Could not convert xyz to floats:{[x,y,z]}")
            coordinates.append([x,y,z])

            tag='BASIS-FUNCTIONS'
            shape_block=get_block(tag, block_range=[start,stop], occurence=1)[0]
            dic=read_basis_functions(*shape_block)
            tag='COEFFICIENTS'
            shape_block=get_block(tag, block_range=[start,stop], occurence=1)[0]
            dic_coeff=read_coefficients(*shape_block)
            dic.update(dic_coeff)
            atoms.append(dic)
        return {'atoms':atoms, 'atom_types':atom_types, 'coordinates': coordinates}
    
class read_camcasp_basis(file_reader):
    comtag='!'
    def __init__(self, inp_fi):
        super().__init__(inp_fi)
        self.dic=self.read_camcasp_basis_file()
    def read_camcasp_basis_file(self):
        type_lines=self.get_line_no(match=r'^[ ]*[0-9][ ]*[spdfghi][ ]*', fullmatch=True)
        value_lines=self.get_line_no(match=r'^[ ]*[0-9]*\.[0-9]*([Ee]?[-]?[0-9]*)[ ]*[0-9]*\.[0-9]*[ ]*', fullmatch=True)

        # Verify that all valid lines (non-empty) have been read in!
        valid_lines=[ l for l in self.lines_num_nocom  if l[1]!='' ]
        for i,l in valid_lines:
            if not i in type_lines+value_lines:
                raise Exception(f"Did not read line {l}")

        link_type_to_value=[]
        for i in type_lines:
            local_value_lines=[]
            for j in range(i+1,len(self.lines)):
                if j in value_lines:
                    local_value_lines.append(j)
                else:
                    break
            if len(local_value_lines)!=1:
                raise Exception(f"Exepected only one line after {l}, found {self.lines[i]}")
            link_type_to_value+=[(i,local_value_lines)]
        link_type_to_value=dict(link_type_to_value)

        exponent_scales, normalizations, exponent_orders, = [ [] for i in range(3) ]
        for i, value_no in link_type_to_value.items():
            type_line=self.lines_nocom[i]
            if not re.match(r'^[ ]*[1][ ]*[s][ ]*$', type_line):
                print('check')
                continue
            for j in value_no:
                value_line=self.lines[j]
                try:
                    exponent_scale,normalization=value_line.split()
                    exponent_scale=float(exponent_scale)
                    normalization=float(normalization)
                    exponent_order=2
                except Exception as e:
                    raise Exception(f"Cannot read two float on the following line: {l} {e}")
                exponent_scales.append(exponent_scale)
                normalizations.append(normalization)
                exponent_orders.append(exponent_order)

        return {
            'exponent_scales':exponent_scales,
            'normalizations':normalizations,
            'exponent_orders': exponent_orders,
        }

from . import *
from .file_readers import file_reader

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
                #print(f"Line {type_line} not 1s function")
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
        }

class read_mom_file(file_reader):
    tag_to_rank={ 0:'mon', 1:'dip', 2:'qad', 3:'oct', 4:'hxd'}
    comtag=['!']
    dic={}
    atom_types=None
    atom_coordinates=None
    site_dmp_di=None
    def __init__(self,inp_fi):
        file_reader.__init__(self,inp_fi)
        self.extract_data()
    def extract_data(self):
        """ The moment file is structure by a commented header and blocks describing sites """
        # Orientation
        empty_lines=self.get_line_no(match=r'[ -]*', fullmatch=True)
        empty_lines+=self.get_line_no(match='', fullmatch=True)
        empty_lines=sorted( list(set(empty_lines)) )
        atom_lines=self.get_line_no(match=r'[A-Z][a-z]?[_]?[0-9]* ', fullmatch=False)
        if len(atom_lines)<1:
            raise Exception(f"Atom lines not found\n{self.input_file}")
        # Check if atom lines were recognized correctly
        for i,e in enumerate(empty_lines):
            if i<len(empty_lines)-1:
                # Skip if consecutive empty line
                if empty_lines[i+1] == e+1:
                    pass
                # In case the document does not end after this line
                elif e<len(self.lines)-1:
                    if not e+1 in atom_lines:
                        raise Exception(f"Line {e} is empty, Total lines {len(self.lines)}:\n{atom_lines}")
            if i==len(empty_lines)-1:
                if e!=len(self.lines)-1:
                    raise Exception

        # Get the atom blocks
        atom_entries=[]
        for a in atom_lines:
            check=False
            for e in empty_lines:
                if e>a:
                    atom_entries.append( self.lines_num_nocom[a:e])
                    check=True
                    break
            if not check:
                raise Exception()

        # Process the intervals
        atoms=[]
        for atom in atom_entries:
            # Read the first line containing meta data
            num,atom_data_line=atom[0]
            ar=atom_data_line.split()

            # Get atom type
            try:
                
                key='type'
                index=[ i for i,x in enumerate(ar) if x.lower()==key]
                assert len(index)==1, f"Expected on key for \'{key}\' (case-insensitive) on line:\n{self.lines_nocom[num]}"
                index=index[0]
                assert len(ar)>index+1, f"Found key \'{key}\' on terminal position but requires following argument"

                at_type=ar[index+1]
            except Exception as ex: 
                at_type='Not_Specified'
        
            # Get coordiantes
            try:
                assert len(ar)>4, f"Expected line of at least 4 (to accomodate coordinates), got: {atom_data_line}"
                xyz=np.array([ float(x) for x in ar[1:4]])
            except Exception as ex:
                raise Exception(f"Error in getting coordiantes for {atom_data_line}: {ex}")

            # Get the rank
            try:
                
                rank_index=[ i for i,x in enumerate(ar) if x.lower()=='rank']
                assert len(rank_index)==1, f"Expected on key for rank on line:\n{self.lines_nocom[num]}"
                rank_index=rank_index[0]
                assert len(ar)>rank_index+1, f"Found rank key on position {rank_index}"

                rank=ar[rank_index+1]
                rank=int(rank)

            except Exception as ex:
                raise Exception(f"Cannot read rank on line {atom_data_line}: {ex}")
            
            # Get the moments
            if len(atom)<2:
                raise Exception(f"No line for moments (atom block ends after one line)\n{atom}")
            # Get the moments as list of components
            try:
                moment_line=[ x[1].strip() for x in atom[1:] ]
                moment_line=[ x.split('!')[0] for x in moment_line ]
                moment_line=' '.join(moment_line).split()
                moment_comps=[ float(x) for x in moment_line ]
                moment_comps=np.array(moment_comps)
            except Exception as e:
                raise Exception(f"Could not read line:\n{atom[1:]}\n{e}")
            # Split the moments
            supposed_length_for_sph_rep=(rank+1)**2
            if len(moment_comps)!=supposed_length_for_sph_rep:
                print(f"Expected length of moment components as {supposed_length_for_sph_rep}, got: {len(moment_comps)}"+
                        "\n{}".format('\n'.join([x[1] for x in atom])))
            start=0
            mom_di={}
            if rank>4:
                print(f"Warning: found rank {rank} while only multipoles through rank 4 can be processed at the moment. Limiting read-in to rank 4.")
                rank=4
            for n in range(rank+1):
                if not n in self.tag_to_rank.keys():
                    raise Exception(f"Did not find key {n} in {self.tag_to_rank}")
                k=self.tag_to_rank[n]
                end=start+(n+1)*2-1
                mom_di.update({k: np.array(moment_comps[start:end]) })
                start=end

            atoms.append({'atom_type':at_type, 'xyz':xyz, 'rank':rank, 'mom_di':mom_di})

        self.atom_types=[ atom['atom_type'] for atom in atoms ]
        self.atom_coordinates=[ atom['xyz'] for atom in atoms ]
        self.site_ranks=[ atom['rank'] for atom in atoms ]
        self.site_dmp_di= [atom['mom_di'] for atom in atoms ]
from . import *

def atomic_charge_to_atom_type(Z):
    """ Return atom type when given an atomic Number """
    # Replacers for Thomas script
    PERIODIC_TABLE_STR = """
    H                                                                                                                           He
    Li  Be                                                                                                  B   C   N   O   F   Ne
    Na  Mg                                                                                                  Al  Si  P   S   Cl  Ar
    K   Ca  Sc                                                          Ti  V   Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr
    Rb  Sr  Y                                                           Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te  I   Xe
    Cs  Ba  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu  Hf  Ta  W   Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn
    Fr  Ra  Ac  Th  Pa  U   Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr  Rf  Db  Sg  Bh  Hs  Mt  Ds  Rg  Cn  Nh  Fl  Mc  Lv  Ts  Og
    """
    atom_types=PERIODIC_TABLE_STR.split()
    Z_to_atty_di=dict([ (i+1, at_ty) for i, at_ty in enumerate(atom_types)])
    assert Z in Z_to_atty_di.keys()
    return Z_to_atty_di[Z]

@validate_call
def load_json(file: str):
    assert os.path.isfile(file), f"Not a file: {file}"
    assert file.lower().endswith('.json'), f"File does end with json: {file}"
    with open(file, 'r') as rd:
        content=json.load(rd)
    return content
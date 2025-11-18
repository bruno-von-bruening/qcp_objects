from . import *

allowed_representation=Literal['S','C']

@val_call
def num_elements_per_rank(rank:int,rep:allowed_representation,symm:bool=False) -> int:
    """ return number of elements per moment"""
    if rep=='S' and not symm:
        raise Exception(f"Confusion in input spherical representation and not symmetric. Does not make sence")

    if rep=='S':
        return 2*rank+1
    else:
        if rank==0:
            return 1
        elif rank==1:
            return 3
        elif rank==2:
            if symm:
                return int( (rank+1)*(rank+2)/2 ) # half the number of elements in matrix and add the trace
            else:
                return 3**rank
        else: raise NotImplementedError(f"Rank {rank} not implememented yet")

@val_call
def triangular_index(i,j):
    """ This is upper (natural for c/python)"""
    if j>i:
        raise Exception(f"Expected upper triangular form, got {(i,j)}")
    return j * (j + 1) // 2 + i
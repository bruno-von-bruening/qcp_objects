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
        else:
            if symm:
                return int( (rank+1)*(rank+2)/2 ) # half the number of elements in matrix and add the trace
            else:
                return 3**rank
        #else: raise NotImplementedError(f"Rank {rank} not implememented yet")

@val_call
def triangular_index(*ints: int) -> int:
    """ Lower triangular matrix!!!"""
    if len(ints) == 0: raise ValueError("triangular_index requires at least one index")

    # Enforce Decreasing order: essential for correct indexing (lower triangular)
    if not all(ints[k] >= ints[k + 1] for k in range(len(ints) - 1)):
        raise Exception(f"Expected non-increasing (symmetric) ordering, got {ints}")

    # Rank-1: just return the index itself
    if len(ints) == 1:
        return ints[0]

    # Rank-2: standard upper-triangular indexing with i >= j
    if len(ints) == 2:
        i, j = ints
        return (i * (i + 1)) // 2 + j
    elif len(ints) == 3:
        i, j, k = ints
        return (i * (i + 1) * (i + 2)) // 6 + (j * (j + 1)) // 2 + k
    elif len(ints) == 4:
        i, j, k, l = ints
        return (i * (i + 1) * (i + 2) * (i + 3)) // 24 + (j * (j + 1) * (j + 2)) // 6 + (k * (k + 1)) // 2 + l
    else:
        raise NotImplementedError("triangular_index only implemented up to rank 4")

    # General rank r >= 3: combinatorial numbering of non-increasing tuples
    from math import comb

    r = len(ints)
    flat = 0
    for k, ik in enumerate(ints):
        flat += comb(ik + r - k, r - k)
    return flat
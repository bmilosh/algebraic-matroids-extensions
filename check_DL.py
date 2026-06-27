import itertools as it
from time import time

from sage.all import *

from lru_cache import LRUCache
from special_matroids import WeakMatroid
from timing import endlog as elog
from timing import log as slog
from timing import secondsToStr


def isModular(M, F1: set, F2: set) -> bool:
    '''
    Takes as input the matroid and the flats F1 and F2.
    Returns boolean indicating if they form a modular pair.
    '''

    rF1 = M.rank(F1)
    rF2 = M.rank(F2)
    rF12 = M.rank(F1.union(F2))
    rF1n2 = M.rank(F1.intersection(F2))
    return rF1 + rF2 == rF12 + rF1n2


def isDLDouble(M, X: set, Y: set):
    """
    Pair is a DL double if there is no T in X such that
    for every W in X, 
        r(T|W) = 0 if and only if r(Y|W) = r(Y|X), 
    where T and W are flats.
    """
    RHS = M.rank(X.union(Y)) - M.rank(X)  #rX
    intr_: frozenset = M.groundset()
    pos_part = []
    for W in sum([[frozenset(itm) for itm in it.combinations(X, i)] for i in range(1, len(X))], []): 
        if W == M.closure(W) and M.rank(Y.union(W)) - M.rank(W) == RHS:
            # W is a flat and r(Y|W) = r(Y|X)
            # we add it to the possible partition and update the intersection (T)
            pos_part.append(W)
            intr_ = intr_.intersection(W)

    # If the intersection of all the flats in the partition also has that r(Y|intr_) == RHS 
    # then we don't have a DL-double, since this intersection is also the
    # quasi-intersection of (X,Y). 
    isDLDb = not (intr_ == M.groundset() or M.rank(Y.union(intr_)) - M.rank(intr_) == RHS)

    if isDLDb and ADD_Y_INTERSECTION_TO_PARTITION:
        pos_part.append(Y.union(intr_))
    if not isDLDb:
        pos_part = []

    return isDLDb, pos_part, intr_



def getDLDouble(M):
    """
    Takes as input a matroid and returns the DL-doubles
    in the matroid along with the corresponding partition and intersection.
    """
    r = M.full_rank()

    if USE_ALL_FLATS:
        flats = [flat for i in range(2, r) for flat in M.flats(i)] 
        for X, Y in it.permutations(flats, 2):
            rX = M.rank(X); rY = M.rank(Y); rXY = M.rank(X.union(Y))
            if rX == 2 or rY == r - 1: 
                # Lemma 5.10
                continue 
            # if M.is_circuit(X.union(Y)) and X.isdisjoint(Y):
            #     # Proposition 5.9
            #     continue
            # if rX + rY == r + 1 and X not in circ_hyps: 
            #     # Proposition 5.11 (but only works for sparse-paving matroids, so needs more conditions)
            #     continue 
            if not isModular(M, X, Y):
                res, flats_, intr_ = isDLDouble(M, X, Y)
                if res:
                    yield X, Y, flats_, intr_

    else:
        for X in M.circuit_closures()[r - 1]:
            for Y in M.flats(2):
                if not isModular(M, X, Y):
                    res, flats_, intr_ = isDLDouble(M, X, Y)
                    if res:
                        yield X, Y, flats_, intr_


def CheckDLAlg(M) -> bool:
    '''
    Matroid staisfies DL if for every DL-double and a partition P,
    there is an extension of M that increases
    the rank of the intersection of all elements in P.
    '''
    for X, Y, partition_, intr_ in getDLDouble(M):  
        if intr_ in M.modular_cut(partition_):
            # No extension exists that increases the rank of the
            # intersection of all elements in the partition.
            # Therefore, the matroid does not satisfy the DL property.
            print(X, Y)
            return False
    return True


def GetExtensions(M, sbsets: list):
    '''
    Given a matroid and a nonmodular pair of flats,
    it returns all the extensions of the matroid in which 
    the corresponding modular cuts contain this pair.
    Output:   All extensions corresponding to sbsets.
              new_e -- The new element.
    '''

    new_e = f"{M.size()}"

    for ext in M.extensions(new_e, subsets=sbsets):
        yield ext, new_e


def RecursiveCheckDLAlg(M, depth=1):
    '''
    The recursive function.
    Takes as input the matroid to be checked
    and the depth at which DL is to be checked.
    Returns boolean indicating if the matroid is 
    DL at given depth.
    '''
    w = WeakMatroid(M)
    if depth != ORIGINAL_DEPTH:
        res = cache[depth].get(w)
        if res != -1:
            return res

    if depth == 1:
        ### We always check DL at depth 1. ###
        res = CheckDLAlg(M) 
        update_cache(w, res, 1)
        return res
    
    
    for X, Y, partition_, intr_ in getDLDouble(M): 
        for_cache = tuple([M] + sorted(
            [tuple(sorted(s)) for s in partition_]
        ))
        res = p_cache[depth].get(for_cache)
        if res == False:
            print(f'# flats are {X} and {Y} at depth {depth}; partition already checked')
            return False
        elif res == True:
            continue

        if intr_ in M.modular_cut(partition_):
            print(partition_)
            print(f'# flats are {X} and {Y} at depth {depth}; not checking these')
            p_cache[depth].put(for_cache, False)
            update_cache(w, False, depth)
            return False
        
        DLmatroid = False  ### The indicator we use when we find an appropriate DL extension. ###
        for N, i in GetExtensions(M, partition_):  
            if N.rank(intr_) == N.rank(intr_.union({i})):  
                # We should only check extensions that increase the rank of the intersection
                # of all elements in the partition.
                continue
            if RecursiveCheckDLAlg(N, depth-1): 
                ### Recursion takes place here. ###
                
                ### The moment we find a proper point extension that is DL we break for that depth. ###
                DLmatroid = True
                p_cache[depth].put(for_cache, True)
                break
        
        if not DLmatroid:
            print(f'# flats are {X} and {Y} at depth {depth}')
            p_cache[depth].put(for_cache, False)
            update_cache(w, False, depth)
            return False
    update_cache(w, True, depth)
    return True


def update_cache(matroid, result, depth):
    if depth != ORIGINAL_DEPTH:
        cache[depth].put(matroid, result)


if __name__ == "__main__":
    ADD_Y_INTERSECTION_TO_PARTITION = False
    # We can choose to check all pairs of flats or only pairs of circuit-hyperplanes and lines.
    USE_ALL_FLATS= True

    bol129075 = {
        ' Bol_129075_9_4_108 ' : ['0135', '0235', '1235', '0145', '0245', '1245', '0345', '1345', '2345', '0136', '0236', '1236', '0146', '0246', '1246', '0346', '1346', '2346', '0156', '0256', '1256', '0356', '1356', '0456', '2456', '3456', '0137', '0237', '1237', '0147', '0247', '1247', '0347', '1347', '2347', '0157', '0257', '1257', '0357', '1357', '2357', '0457', '1457', '3457', '0167', '0267', '1267', '0367', '2367', '0467', '1467', '2467', '3467', '1567', '2567', '3567', '4567', '0138', '0238', '1238', '0148', '0248', '1248', '0348', '1348', '2348', '0158', '0258', '1258', '0358', '1358', '2358', '0458', '1458', '2458', '3458', '0168', '0268', '1268', '0368', '1368', '2368', '0468', '1468', '2468', '0568', '1568', '2568', '3568', '4568', '0178', '0278', '1278', '0378', '1378', '2378', '1478', '2478', '3478', '0578', '1578', '2578', '4578', '0678', '1678', '3678', '4678', '5678']
    }

    mat_bol075 = Matroid(bases=bol129075[' Bol_129075_9_4_108 '], groundset='012345678')

    ORIGINAL_DEPTH = 6
    mat = "mat_bol075"

    cache = {i: LRUCache(2048) for i in range(1, ORIGINAL_DEPTH)}
    p_cache = {i: LRUCache(2048) for i in range(2, ORIGINAL_DEPTH + 1)}

    begin = time()
    slog(f'Start {ORIGINAL_DEPTH}-DL check for {mat}')

    with open(f"DL_run_{mat}_test3.log", "a+") as sfile:
        msg = f"{secondsToStr()} Solution file for {ORIGINAL_DEPTH}-DL check for matroid {mat}."
        wrapper = "%" * len(msg)
        sfile.write(wrapper + "\n")
        sfile.write(msg + "\n")
        sfile.write(wrapper + "\n")

        cgp = RecursiveCheckDLAlg(mat_bol075, ORIGINAL_DEPTH)
        print(cgp)
        if cgp:
            sfile.write(f"\n{secondsToStr()} Matroid {mat} is {ORIGINAL_DEPTH}-DL.\n\n")
        else:
            sfile.write(f"\n{secondsToStr()} Matroid {mat} is not {ORIGINAL_DEPTH}-DL.\n\n")
    elog(begin)
    for i in range(1, ORIGINAL_DEPTH):
        print(f"depth {i}: ", end=" ")
        print(cache[i].cache_info())
    for i in range(2, ORIGINAL_DEPTH + 1):
        print(f"partition cache depth {i}: ", end=" ")
        print(p_cache[i].cache_info())


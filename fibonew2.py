import itertools as it

from gurobipy import *

import config


def bs(v):
    s=set([])
    for j in range(0,config.vrbls):
        if v[j]==1:
            s.add(j)
    return s

# bi: passa un vector binari a enter
def bi(v):
    i=0
    for j in range(0,config.vrbls):
        i+=v[j]*2**j
    return i

# ib: passa un enter a vector binari
def ib(n):
    v=[]
    i=1
    for i in range(1,config.vrbls+1):
        v=v+[n%2]
        n=(n-n%2)//2
        i+=1
    return v

# sb: conjunt a vector binari
def sb(s):
    s = [int(it) for it in s]
    u=[0]*config.vrbls
    for i in range(0,config.vrbls): # [0..vars-1]:
    	if i in s:
            u[i]=1
    return u

def si(s):
    # i = 0
    return sum(2**int(itm) for itm in s)

# union: calcula l'enter corresponent a la unio de dos conjunts
def union(i,j):
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls):
        if v[k]==1:
            u[k]=1
    return(bi(u))

def contained(i,j):
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls):
        if u[k]==1 and v[k]==0:
            return False
    return True

def disjoint(i,j):
    u=ib(i)
    v=ib(j)
    for k in range(0,config.vrbls): # [0..vars-1]:
        if u[k]==v[k] and u[k]==1:
            return False
    return True

#####################################################
# LP function with some new and specific conditions #
#####################################################

def Init1m():
    config.p.setObjective(config.w[2**config.vrbls],GRB.MINIMIZE)
    #config.p.setParam("Presolve", 2)
    #currentlimit = p.Params.presolve
    #print("current presolve value is ",currentlimit)
    for i in range(1,config.Part): ####### [1,...,Part-1]:
        config.p.addConstr(config.w[2**config.vrbls]>=config.w[2**i])
    config.p.addConstr(config.w[0]==0)
    config.p.addConstr(config.w[1]==1)
    #config.p.update()


def Resol2m(output_flag=0, should_write=False, file_path=None):
    config.p.setParam('OutputFlag', output_flag)
    #config.p.setParam('MIPFocus',0)
    # config.p.write("test_polymat.mps")
    if should_write:
        config.p.write(file_path)
    res = config.p.optimize()
    # res = config.p.computeIIS()
    return res


def add_symmetry_conditions(perms: list[list[set[int]]], part_or_vrbls=config.Part):
    # Perms should have same rank
    for p in perms:
        config.p.addConstr(config.w[sum(2**j for j in p[0])] == config.w[sum(2**k for k in p[1])])
    
    # Supersets of perms should have same rank
    for s in sum([list(it.combinations(range(1, part_or_vrbls), i)) for i in range(2, part_or_vrbls)], []):
        X = set(s)
        for p in perms:
            if len(X) > max(len(p[0]), len(p[1])) and X.issuperset(p[0]) and X.isdisjoint(p[1]):
                config.p.addConstr(config.w[sum(2**j for j in X)] == config.w[sum(2**k for k in (X - p[0] | p[1]))])


def add_symmetry_conditions_new(perms: list[list[dict[int]]], part_or_vrbls=config.Part):
    # # Not sure about this
    # # Perms should have same rank
    # for p in perms:
    #     config.p.addConstr(config.w[sum(2**j for j in p[0])] == config.w[sum(2**k for k in p[1])])
    
    # Supersets of perms should have same rank
    for s in sum([list(it.combinations(range(1, part_or_vrbls), i)) for i in range(1, part_or_vrbls)], []):
        X = set(s)
        for p in perms:
            if not X.isdisjoint(set(p[0].keys())) or not X.isdisjoint(set(p[1].keys())):
                config.p.addConstr(config.w[sum(2**j for j in X)] == config.w[sum(2**k for k in [p[0].get(itm, itm) for itm in X])])


# Shannon0: H(S_i)>=0
def Shannon0():
    for i in range(0,config.vrbls):
        config.p.addConstr(config.w[2**i]>=0)


# Shannon1: H(S_Q)\geq H(S_{Q-i}) per tot i
def Shannon1():
    for i in range(0,config.vrbls):
        config.p.addConstr(config.w[2**(config.vrbls)-1]>=config.w[2**(config.vrbls)-1-2**i])


# Shannon2: Condicio de Submodularitat
def Shannon2():
    for i in range(0,2**config.vrbls):
        v=ib(i)
        for j in range(0,config.vrbls):
            if v[j]==0:
                for k in range(j+1,config.vrbls):
                    if v[k]==0:
                        config.p.addConstr(config.w[i+2**j]+config.w[i+2**k]>=config.w[i]+config.w[i+2**j+2**k])


def Shannonm():
	Shannon0()
	Shannon1()
	Shannon2()


def CE(v,i1,i2):
	i12=union(i1,i2)
	s=v[i12]-v[i2]
	return s


def MIC(v,i1,i2,i3):
	i13=union(i1,i3)
	i23=union(i2,i3)
	i123=union(i1,i23)
	s=v[i13]+v[i23]-v[i123]-v[i3]
	return s


def AKNew_def(X: set, Y: set, l: int):
    if isinstance(X, set):
        x_int = bi(sb(X))
        y_int = bi(sb(Y))
    else:
        x_int = X
        y_int = Y
        X = bs(ib(X))
        Y = bs(ib(Y))
    config.p.addConstr(CE(config.w, l, x_int) == 0) # g(Z|X) = 0
    # g(X'|Z) = g(X'|Y) for every X' subset of X
    subsets_of_X = sum([list(it.combinations(X, i)) for i in range(1, len(X) + 1)], [])
    for subset in subsets_of_X:
        sub_int = bi(sb(subset))
        config.p.addConstr(CE(config.w, sub_int, l) == CE(config.w, sub_int, y_int))


def MI(v,i1,i2):
	i12=union(i1,i2)
	s=v[i1]+v[i2]-v[i12]
	return s


def CI(i,j,k):
    config.p.addConstr(CE(config.w,k,i)==0)
    config.p.addConstr(CE(config.w,k,j)==0)
    config.p.addConstr(MI(config.w,i,j)==config.w[k])


def setgenerator(structure: dict):
    '''
    structure is the dictionary containing the matroid ports
    '''
    res = dict()
    for key, min_AS in structure.items():
        conv_AS = set()
        for min_set in min_AS:
            h = 0
            for num in min_set:
                h += 2**int(num)
            conv_AS.add(h)
        # de = set(crowd)
        res[key] = conv_AS
    return res


def AccStrCompatiblemnew(Ad):
    for j in range(0,2**config.Part-1,2): ######## [0,2,...,2**Part-2]:
        t=1
        for k in Ad:  
            #k = sum([2**int(ite) for ite in list(kin) ])
            if contained(k,j):
                t=0
        config.p.addConstr(config.w[j+1]-config.w[j]==t)

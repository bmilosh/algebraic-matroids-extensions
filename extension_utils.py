# from gurobipy import *
import gurobipy as guro
import typing
import itertools as it
import os
import sys


from sage.all import *


def binary_to_int(vector: typing.List[int], variables: int) -> int:
    """
    Converts binary vector to integer.
    `vector` has to have the same length as variables
    """
    res = 0
    for j in range(variables):
        res += (vector[j] * (2**j))
    return res


def sequence_to_int(s: typing.Sequence) -> int:
    """
    Converts a Python sequence to an integer.
    That is, given, say, [1, 4] it should return
    2**1 + 2**4
    Goes without saying that elements of sequence
    have to be numeric-type.
    """
    return sum(2**(int(itm)) for itm in s)


def int_to_binary_vector(num: int, variables: int) -> typing.List[int]:
    """
    COnverts an integer to a binary vector.
    2**num should be less than 2**variables to get
    a correct result
    """
    vector = []
    for _ in range(variables):
        vector += [num % 2]
        num = (num - (num % 2)) // 2
    return vector


def int_union(num1: int, num2: int, variables: int) -> int:
    """
    Computes the "union" of two integers.
    Does this by first converting them to binary vectors
    and working on the vectors.
    """
    v1 = int_to_binary_vector(num1, variables)
    v2 = int_to_binary_vector(num2, variables)
    for k in range(variables):
        if v2[k] == 1:
            v1[k] = 1
    return binary_to_int(v1, variables)


def conditional_entropy(vector, num1: int, num2: int, variables: int) -> int:
    union = int_union(num1, num2, variables)
    return vector[union] - vector[num2]


def mutual_information(vector, num1: int, num2: int, variables: int) -> int:
    union = int_union(num1, num2, variables)
    return vector[num1] + vector[num2] - vector[union]



class Utils:
    def __init__(self, model, constraint_vector, participants: int, variables: int) -> None:
        self.participants = participants
        self.variables = variables
        self.model = model
        self.constraint_vector = constraint_vector

    # @staticmethod
    def add_initialization_constraint(self):
        """
        `InitMatNew` in the former setup.
        Basically just adds the Shannon constraints
        """
        # Set the objective
        self.model.setObjective(0, guro.GRB.MINIMIZE)
        # Add Shannon inequality constraints
        for i in range(self.variables):
            # entropy is nonnegative
            self.model.addConstr(self.constraint_vector[2**i] >= 0) 
            # entropy of groundset > groundset less an element
            self.model.addConstr(self.constraint_vector[2**self.variables - 1] >= self.constraint_vector[(2**self.variables - 1) - (2**i)])
        
        #############################
        ## Submodularity condition ##
        #############################
        for i in range(2**self.variables):
            v = int_to_binary_vector(i, self.variables)
            for j in range(self.variables):
                if v[j] == 0:
                    for k in range(j + 1, self.variables):
                        if v[k] == 0:
                            self.model.addConstr(
                                self.constraint_vector[i + 2**j] + self.constraint_vector[i + 2**k] >= 
                                self.constraint_vector[i] + self.constraint_vector[i + 2**j + 2**k]
                            )
                            # self.submod_constraints_list.append([i + 2**j, i + 2**k, i, i + 2**j + 2**k])


    # @staticmethod
    def add_matroid_compatibility_constraints(self, matroid, groundset):
        """
        `MatroidCompatible` in the former setup
        Need to import sage for this to work
        """
        for i in range(len(groundset) + 1):
            for s in it.combinations(groundset, i):
                j = sequence_to_int(s)
                self.model.addConstr(self.constraint_vector[j] == matroid.rank(s))

    
    def CI(self, seq1: typing.Sequence, seq2: typing.Sequence, aux: int):
        i1 = sequence_to_int(seq1)
        i2 = sequence_to_int(seq2)
        self.model.addConstr(conditional_entropy(self.constraint_vector, aux, i1, self.variables) == 0)
        self.model.addConstr(conditional_entropy(self.constraint_vector, aux, i2, self.variables) == 0)
        self.model.addConstr(mutual_information(self.constraint_vector, i1, i2, self.variables) == self.constraint_vector[aux])

    
    def AK2(self, seq1: typing.Sequence, seq2: typing.Sequence, aux: int, matroid):
        i1 = sequence_to_int(seq1)
        i2 = sequence_to_int(seq2)
        self.model.addConstr(conditional_entropy(self.constraint_vector, aux, i1, self.variables) == 0)
        # sum([list(combinations(groundset, i)) for i in range(len(groundset)+1)], [])
        subsets = sum([list(it.combinations(seq1, num)) for num in range(1, len(seq1) + 1)], [])
        for sub in subsets:
            # Check if it's a flat
            # print(sub, matroid.rank(sub))
            # # print(matroid.flats(matroid.rank(sub)))
            # print(sub in [ f for f in matroid.flats(matroid.rank(sub))])
            # print(sub in matroid.flats(matroid.rank(sub)))
            # print(frozenset(sub) in matroid.flats(matroid.rank(sub)))
            # input("Waiting")
            # if frozenset(sub) in matroid.flats(matroid.rank(sub)):
            #     # add the constraint
            i1_sub = sequence_to_int(sub)
            self.model.addConstr(
                conditional_entropy(self.constraint_vector, i1_sub, aux, self.variables) == 
                conditional_entropy(self.constraint_vector, i1_sub, i2, self.variables) 
            )
        # self.model.addConstr(mutual_information(self.constraint_vector, i1, i2, self.variables) == self.constraint_vector[aux])


    # @staticmethod
    def solve_model(self, output_flag=0, should_write=False, file_path=None):
        """
        `Resol2m` in the former setup
        """
        self.model.setParam("OutputFlag", output_flag)
        # path = 'C:\\Users\\mikky\\OneDrive\\Escritorio\\For Oriol\\v8_AK2.lp'
        # if sys.platform == 'linux':
        #     path = os.popen(f'wslpath "{path}"').read().strip()

        # self.model.write(path)
        if should_write:
            self.model.write(file_path)

        result = self.model.optimize()

        return result



# class Utils:
#     def __init__(self, participants: int, variables: int) -> None:
#         self.participants = participants
#         self.variables = variables
#         self.model = None
#         self.constraint_vector = None
#         # self.init_constr_added = False
#         # self.submod_constraints_list = []
#         # self.mat_constr_added = False
#         # self.mat_constr_list = []

#     # @staticmethod
#     def add_initialization_constraint(self):
#         """
#         `InitMatNew` in the former setup.
#         Basically just adds the Shannon constraints
#         """
#         # Set the objective
#         self.model.setObjective(0, GRB.MINIMIZE)
#         # Add Shannon inequality constraints
#         for i in range(self.variables):
#             # entropy is nonnegative
#             self.model.addConstr(self.constraint_vector[2**i] >= 0) 
#             # entropy of groundset > groundset less an element
#             self.model.addConstr(self.constraint_vector[2**self.variables - 1] >= self.constraint_vector[(2**self.variables - 1) - (2**i)])
        
#         #############################
#         ## Submodularity condition ##
#         #############################
#         # if self.init_constr_added:
#         #     for l in self.submod_constraints_list:
#         #         self.model.addConstr(
#         #             self.constraint_vector[l[0]] + self.constraint_vector[l[1]] >= 
#         #             self.constraint_vector[l[2]] + self.constraint_vector[l[3]]
#         #         )
#         #         # self.submod_constraints_list.append([i + 2**j, i + 2**k, i, i + 2**j + 2**k])

#         # else:
#         #     for i in range(2**self.variables):
#         #         v = int_to_binary_vector(i, self.variables)
#         #         for j in range(self.variables):
#         #             if v[j] == 0:
#         #                 for k in range(j + 1, self.variables):
#         #                     if v[k] == 0:
#         #                         self.model.addConstr(
#         #                             self.constraint_vector[i + 2**j] + self.constraint_vector[i + 2**k] >= 
#         #                             self.constraint_vector[i] + self.constraint_vector[i + 2**j + 2**k]
#         #                         )
#         #                         self.submod_constraints_list.append([i + 2**j, i + 2**k, i, i + 2**j + 2**k])
#         #     self.init_constr_added = True
#         for i in range(2**self.variables):
#             v = int_to_binary_vector(i, self.variables)
#             for j in range(self.variables):
#                 if v[j] == 0:
#                     for k in range(j + 1, self.variables):
#                         if v[k] == 0:
#                             self.model.addConstr(
#                                 self.constraint_vector[i + 2**j] + self.constraint_vector[i + 2**k] >= 
#                                 self.constraint_vector[i] + self.constraint_vector[i + 2**j + 2**k]
#                             )
#                             # self.submod_constraints_list.append([i + 2**j, i + 2**k, i, i + 2**j + 2**k])


#     # @staticmethod
#     def add_matroid_compatibility_constraints(self, matroid, groundset):
#         """
#         `MatroidCompatible` in the former setup
#         Need to import sage for this to work
#         """
#         # if self.mat_constr_added:
#         #     for l in self.mat_constr_list:
#         #         self.model.addConstr(self.constraint_vector[l[0]] == l[1])
#         # else:
#         #     for i in range(len(groundset) + 1):
#         #         for s in it.combinations(groundset, i):
#         #             j = sequence_to_int(s)
#         #             r = matroid.rank(s)
#         #             self.model.addConstr(self.constraint_vector[j] == r)
#         #             self.mat_constr_list.append([j, r])
#         #     self.mat_constr_added = True
#         for i in range(len(groundset) + 1):
#             for s in it.combinations(groundset, i):
#                 j = sequence_to_int(s)
#                 self.model.addConstr(self.constraint_vector[j] == matroid.rank(s))

    
#     def CI(self, seq1: typing.Sequence, seq2: typing.Sequence, aux: int):
#         i1 = sequence_to_int(seq1)
#         i2 = sequence_to_int(seq2)
#         self.model.addConstr(conditional_entropy(self.constraint_vector, aux, i1, self.variables) == 0)
#         self.model.addConstr(conditional_entropy(self.constraint_vector, aux, i2, self.variables) == 0)
#         self.model.addConstr(mutual_information(self.constraint_vector, i1, i2, self.variables) == self.constraint_vector[aux])


#     # @staticmethod
#     def solve_model(self):
#         """
#         `Resol2m` in the former setup
#         """
#         self.model.setParam("OutputFlag", 0)
#         # path = 'C:\\Users\\mikky\\OneDrive\\Escritorio\\For Oriol\\v8_1.lp'
#         # if sys.platform == 'linux':
#         #     path = os.popen(f'wslpath "{path}"').read().strip()

#         # self.model.write(path)

#         result = self.model.optimize()

#         return result
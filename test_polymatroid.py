import gurobipy as guro
from sage.all import *

import config
from extension_utils import Utils
from fibonew2 import ib


def find_extension(method, combos, bases_dict: dict, file_path, gset, should_write=False, write_res_to=None):
    with open(file_path, "w+") as sol_file:
        for key, _ in bases_dict.items():

            with guro.Model(env=env) as config.p:
                sol_file.write(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
                sol_file.write(f"Checking {method} for {key}.\n\n")
                sol_file.write(f"Variables: {config.Part = }, {config.aux = }, {config.vrbls = }\n\n")
                sol_file.write(f"Sets: {combos}\n\n")
                config.w = config.p.addVars(range(0, 2**config.vrbls + 1), name="w")
                # config.p.setParam("Outputflag", 0.0)
                config.p.setParam("Method", 2.0)
                config.p.setParam("Crossover", 0.0)
                # config.p.setParam("BarHomogeneous", 1.0)

                utils = Utils(config.p, config.w, config.Part, config.vrbls)

                # Add model constraints
                utils.add_initialization_constraint()
                utils.add_matroid_compatibility_constraints(original_matroid, gset)
                if method == "CI":
                    for idx, items in enumerate(combos):
                        utils.CI(set(items[0]), set(items[1]), 2**(config.Part + idx))
                else:
                    for idx, items in enumerate(combos):
                        utils.AK2(set(items[0]), set(items[1]), 2**(config.Part + idx), original_matroid)

                utils.solve_model(output_flag=1, should_write=should_write, file_path=write_res_to)

                try:
                    sol_file.write(f"#### Objective value: {config.p.objVal}\n")
                    sol_file.write(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")
                    for i in range(2**config.vrbls):
                        bin_rep = ib(i)
                        rank = config.w[i].X
                        sol_file.write(f"{bin_rep} --> {rank}\n")
                    return True
                except Exception as e:
                    print(f"We got error: {e}")
                    print(f"Model status: {config.p.status}")
                    sol_file.write(f"We got error: {e}\n")
                    sol_file.write(f"Model status: {config.p.status}\n")
                    return False


Gset9 = [f"{i}" for i in range(9)]
Gset10 = [f"{i}" for i in range(10)]

ttt = {'Tictactoe' : ['01234','01246','01245','01248','01235','01238','01237','01256','01268','01267','01257','01278','02346','02345','02348','02347',
'02456','02468','02467','02458','02457','02478','02356','02368','02367','02358','02357','02378','02568','02567','02678','02578',
'01346','01345','01348','01347','01456','01468','01467','01458','01457','01478','01356','01368','01367','01358','01357','01378',
'01568','01567','01678','01578','03468','03467','03458','03457','03478','04568','04567','04678','04578','03568','03567','03578',
'05678','12346','12345','12348','12347','12456','12468','12467','12458','12457','12478','12356','12368','12367','12358','12357',
'12378','12568','12567','12678','12578','23456','23468','23467','23457','23478','24568','24567','24678','24578','23568','23567',
'23678','23578','13456','13468','13467','13458','13457','13478','14568','14567','14578','13568','13567','13678','13578','15678',
'34568','34567','34678','34578','45678','35678']
}

bol129075 = {
' Bol_129075_9_4_108 ' : ['0135', '0235', '1235', '0145', '0245', '1245', '0345', '1345', '2345', '0136', '0236', '1236', '0146', '0246', '1246', '0346', '1346', '2346', '0156', '0256', '1256', '0356', '1356', '0456', '2456', '3456', '0137', '0237', '1237', '0147', '0247', '1247', '0347', '1347', '2347', '0157', '0257', '1257', '0357', '1357', '2357', '0457', '1457', '3457', '0167', '0267', '1267', '0367', '2367', '0467', '1467', '2467', '3467', '1567', '2567', '3567', '4567', '0138', '0238', '1238', '0148', '0248', '1248', '0348', '1348', '2348', '0158', '0258', '1258', '0358', '1358', '2358', '0458', '1458', '2458', '3458', '0168', '0268', '1268', '0368', '1368', '2368', '0468', '1468', '2468', '0568', '1568', '2568', '3568', '4568', '0178', '0278', '1278', '0378', '1378', '2378', '1478', '2478', '3478', '0578', '1578', '2578', '4578', '0678', '1678', '3678', '4678', '5678'] ,
}

# MTTT = Matroid(ttt['Tictactoe'])
# tttd = MTTT.dual()
# MTTT_ext = tttd.coextension('9', subsets=list(tttd.circuit_closures()[3]))
# ttt_ext_bases = {'TTT_and_dual_coext' : [''.join(b) for b in MTTT_ext.bases()]}

original_matroid = Matroid(bol129075[' Bol_129075_9_4_108 '])

config.Part = 9
config.aux = 6
config.vrbls = config.Part + config.aux

# # These work!!!
# X = {'5', '9', '7', '4', '8'}; Y = {3, 6}    # for TTT_coext
# X1 = {'5', '2', '9', '4', '1'}; Y1 = {0, 3}  # for TTT_coext
# X2 = {2, 5, 8, 9, 10, 11}; Y2 = {1, 7}  # for TTT_coext

# for 129075  # This works!
combos = [
    ({8, 2, 7, 6}, {3, 5}, 0),        
    ({3, 9, 2, 5, 6}, {1, 4}, 1),   
    ({3, 0, 10, 2, 1, 4}, {5, 7}, 2), 
    ({3, 0, 10, 11, 2, 1, 4}, {6, 7}, 3),
    ({3, 8, 6, 4}, {9, 11}, 4),
    ({7, 9, 8, 2, 12, 6}, {10, 13}, 5),
] 


path = "test_poly_solution_bol129075_6AK.log"


with guro.Env() as env: 
    find_extension("AK", combos, bol129075, path, Gset9) 

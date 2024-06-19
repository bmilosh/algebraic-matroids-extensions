import gurobipy as guro

import config
from fibonew2 import (CI, AccStrCompatiblemnew, AKNew_def, Init1m, Resol2m,
                      Shannonm, add_symmetry_conditions, bi, ib, sb,
                      setgenerator)


def check_bound(method, combos, ports_dict, file_path, should_write=False, write_res_to=None, 
                use_spec=False, part_or_vrbls=config.Part, perms: list = None):
    generated_sets = setgenerator(ports_dict)
    with open(file_path, "w+") as sol_file:
        for key in generated_sets:
            access_structure = generated_sets[key]
            print(access_structure)

            with guro.Model(env=env) as config.p:
                sol_file.write(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
                sol_file.write(f"Checking {method} for {key}.\n\n")
                sol_file.write(f"Variables: {config.Part = }, {config.aux = }, {config.vrbls = }\n\n")
                sol_file.write(f"Sets: {combos}\n\n")
                config.w = config.p.addVars(range(0, 2**config.vrbls + 1), name="w")
                Init1m()
                Shannonm()
                AccStrCompatiblemnew(access_structure)
                if use_spec:
                    add_symmetry_conditions(perms, part_or_vrbls=part_or_vrbls)
                if method == "CI":
                    for X, Y, aux in combos:
                        CI(bi(sb(X)), bi(sb(Y)), 2**(config.Part + aux))
                else:
                    for X, Y, aux in combos:
                        AKNew_def(set(X), set(Y), 2**(config.Part + aux))
                Resol2m(output_flag=1, should_write=should_write, file_path=write_res_to)
                try:
                    sol_file.write(f"#### Objective value: {config.p.objVal}\n")
                    sol_file.write(f"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n")
                    for i in range(2**config.vrbls):
                        sol_file.write(f"{ib(i)} --> {config.w[i].X}\n")
                    return True
                except Exception as e:
                    print(f"We got error: {e}")
                    print(f"Model status: {config.p.status}")
                    return False


bol_075_ports = {
'Bol_129075_9_4_108 port 0' : ['12', '134', '1356', '1357', '1358', '1368', '1378', '1457', '1458', '1467', '1468', '1568', '1578', '1678', '234', '2357', '2358', '2367', '2368', '2378', '2456', '2458', '2467', '2468', '2568', '2578', '3456', '3457', '3458', '3467', '3568', '3678', '4568', '478', '567'] ,
}

Gset9 = [0, 1, 2, 3, 4, 5, 6, 7, 8]


# for 129075  # This works!
combos = [
    ({8, 2, 7, 6}, {3, 5}, 0),        
    ({3, 9, 2, 5, 6}, {1, 4}, 1),   
    ({3, 0, 10, 2, 1, 4}, {5, 7}, 2), 
    ({3, 0, 10, 11, 2, 1, 4}, {6, 7}, 3),
    ({3, 8, 6, 4}, {9, 11}, 4),
    ({7, 9, 8, 2, 12, 6}, {10, 13}, 5),
] 


config.Part = 9
config.aux = 6
config.vrbls = config.Part + config.aux


path = "test_ss_solution_bol075p0_6aux_AK.log"


with guro.Env() as env: 
    check_bound("AK", combos, bol_075_ports, path)

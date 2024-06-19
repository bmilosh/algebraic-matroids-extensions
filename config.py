import gurobipy as guro

Part = 10
aux = 3
vrbls = Part + aux

with guro.Env() as env, guro.Model(env=env) as p:
    w = p.addVars(range(0, 2**vrbls + 1), vtype=guro.GRB.INTEGER, name="w") # When we want only matroid extensions
    
    # # Needed for checking Dress-Lovasz extensions
    # b = p.addVar(vtype=guro.GRB.BINARY, name="b")

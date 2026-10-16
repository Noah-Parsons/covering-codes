# Direct k-rainbow domination number of Q_n (graph-theory definition), by SCIP.
import sys
from pyscipopt import Model, quicksum
def rainbow(n, k, tl=300):
    V = 1 << n
    m = Model(); m.hideOutput(); m.setParam("limits/time", tl)
    x = [[m.addVar(vtype="B") for a in range(k)] for v in range(V)]
    e = [m.addVar(vtype="B") for v in range(V)]           # e=1: f(v) empty
    for v in range(V):
        for a in range(k):
            m.addCons(x[v][a] <= 1 - e[v])
            m.addCons(quicksum(x[v ^ (1 << i)][a] for i in range(n)) >= e[v])
        m.addCons(quicksum(x[v]) >= 1 - e[v])
    m.setObjective(quicksum(x[v][a] for v in range(V) for a in range(k)), "minimize")
    m.optimize()
    return m.getStatus(), round(m.getObjVal()), m.getDualbound()
for k, ns in ((2, range(1, 7)), (3, range(1, 6))):
    for n in ns:
        st, val, bd = rainbow(n, k)
        print(f"gamma_r{k}(Q_{n}) = {val}  ({st}, bound {bd:.2f})", flush=True)

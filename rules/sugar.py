from sympy.polys.monomials import monomial_lcm

def total_degree(f):
    return max(sum(m) for m in f.monoms())

def pair_sugar(i, j, G, sugars):
    lcm = monomial_lcm(G[i].LM, G[j].LM)
    return max(sugars[i] - sum(G[i].LM), sugars[j] - sum(G[j].LM)) + sum(lcm)

def sugar_selection(B, G, state):
    order = state["ring"].order
    sugars = state["sugar"]
    lcms = state.setdefault("lcm", {})
    pair_sugars = state.setdefault("pair_sugar", {})

    def key(k):
        i, j = B[k]
        if (i, j) not in lcms:
            lcms[(i, j)] = monomial_lcm(G[i].LM, G[j].LM)
        if (i, j) not in pair_sugars:
            pair_sugars[(i, j)] = pair_sugar(i, j, G, sugars)
        return (pair_sugars[(i, j)], order(lcms[(i, j)]), i, j)

    return min(range(len(B)), key=key)
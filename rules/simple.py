from sympy.polys.monomials import monomial_lcm

def first_selection(B, G, state):
    return min(range(len(B)), key=lambda k: (B[k][1], B[k][0]))   # smallest j, then smallest i

def random_selection(B, G, state):
    return state["rng"].randrange(len(B))

def normal_selection(B, G, state):
    order = state["ring"].order
    lcms = state.setdefault("lcm", {})

    def key(k):
        i, j = B[k]
        if (i, j) not in lcms:
            lcms[(i, j)] = monomial_lcm(G[i].LM, G[j].LM)
        return (order(lcms[(i, j)]), j, i)

    return min(range(len(B)), key=key)

def degree_selection(B, G, state):
    order = state["ring"].order
    lcms = state.setdefault("lcm", {})

    def key(k):
        i, j = B[k]
        if (i, j) not in lcms:
            lcms[(i, j)] = monomial_lcm(G[i].LM, G[j].LM)
        m = lcms[(i, j)]
        return (sum(m), j, i)

    return min(range(len(B)), key=key)
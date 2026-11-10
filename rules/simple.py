def first_selection(B, G, state):
    return min(range(len(B)), key=lambda k: B[k])   # smallest (i, j)

def random_selection(B, G, state):
    return state["rng"].randrange(len(B))

def normal_selection(B, G, state):
    return B[0]     # TODO

def degree_selection(B, G, state):
    return B[0]     # TODO
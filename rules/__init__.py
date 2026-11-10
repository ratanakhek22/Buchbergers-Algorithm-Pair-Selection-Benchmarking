from .simple import first_selection, random_selection, normal_selection, degree_selection
from .sugar import sugar_selection

RULES_MAP = {
    "first": first_selection,
    "random": random_selection,
    "normal": normal_selection,
    "degree": degree_selection,
    "sugar": sugar_selection,
}
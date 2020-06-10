import random

import pyrbn
import pyrbnachem


def test_rbnmol():
    rbn_a = pyrbn.RBN.from_random(random.Random(42))
    rbn_b = pyrbn.RBN.from_random(random.Random(42))

    mol_a = pyrbnachem.RBNMol([rbn_a], {})
    mol_b = pyrbnachem.RBNMol([rbn_b], {})
    mol_ab = pyrbnachem.RBNMol([rbn_a, rbn_b], {})
    mol_ab = pyrbnachem.RBNMol([rbn_a, rbn_b], {})

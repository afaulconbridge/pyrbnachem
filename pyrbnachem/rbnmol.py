import itertools
import logging
from typing import Sequence
from pyrbn import RBN


class RBNMol:
    rbn = None
    atoms = []
    bonds = {}

    def __init__(self, atoms: Sequence[RBN], bonds):
        assert 0 < len(atoms), "must include at least one atom"
        assert len(bonds) == 0 or 0 >= min(
            itertools.chain(*bonds.keys())
        ), "bond to atom position less than zero"
        assert len(bonds) == 0 or len(atoms) > max(
            itertools.chain(*bonds.keys())
        ), f"bond to atom position larger than count {len(atoms)} {bonds}"

        # order atoms
        # zip atoms with an original postion, then sort them, then separate back into two parallel lists
        atom_numbers, atoms = zip(*sorted(enumerate(atoms), key=lambda x: x[1]))

        # update bonds as appropriate
        bonds_fixed = {}
        for key, value in bonds.items():
            i_a, i_b = key

            # use sorted atom numbering
            i_a_fixed = atom_numbers[i_a]
            i_b_fixed = atom_numbers[i_b]

            # swap if necessary
            if i_b_fixed < i_a_fixed:
                key_fixed = (i_a_fixed, i_b_fixed)
                # also swap value
                value_fixed = (value[1], value[0])
            else:
                key_fixed = (i_a_fixed, i_b_fixed)
                value_fixed = value

            # save into the new bond dictionary
            bonds_fixed[key_fixed] = value_fixed

        # change to use the new bond dictionary
        bonds = bonds_fixed

        # now that atoms and bonds have been tidied, store them
        # note that this isn't perfect - imagine a molecule with two of the same atoms but differnt bonding sites
        self.atoms = tuple(atoms)
        self.bonds = dict(bonds)

        # calculate rbn based on given atoms and bonds
        rbn_states = sum((a.states for a in atoms), ())
        rbn_funcs = sum((a.funcs for a in atoms), ())
        rbn_inputs = []
        # inputs must have extra offsets
        offset = 0
        for a in atoms:
            for inputs in a.inputs:
                inputs_offset = []
                for inpt in inputs:
                    inputs_offset.append(inpt + offset)
                rbn_inputs.append(list(inputs_offset))
            offset += a.n

        # apply bond changes
        for key in bonds:
            n_a, n_b = key
            (i_a, j_a), (i_b, j_b) = bonds[key]
            rbn_input_a = sum((x.n for x in atoms[:n_a])) + i_a
            rbn_input_b = sum((x.n for x in atoms[:n_b])) + i_b
            rbn_inputs[rbn_input_a][j_a] = rbn_input_b
            rbn_inputs[rbn_input_b][j_b] = rbn_input_a

        # initial RBN
        self.rbn = RBN(rbn_states, rbn_inputs, rbn_funcs)

        # calculate cycle
        self.cycle = self.rbn.get_cycle()

        # go to lowest state of cycle
        self.rbn = RBN(min(self.cycle), self.rbn.inputs, self.rbn.funcs)

    @classmethod
    def from_random(clzz, rng, n=5, k=2):
        to_return = clzz([RBN.from_random(rng, n, k)], {})
        return to_return

    def __repr__(self):
        return f"{self.__class__.__name__}({self.atoms}, {self.bonds})"

    # TODO __hash__
    # TODO __eq__
    # TODO __lt__

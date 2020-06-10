import abc
import itertools
import collections
from pyachem import AChem, Reaction
from .rbnmol import RBNMol


class RBNAchem(AChem):
    def _score(self, mol, n):
        n_pre = sum((x.n for x in mol.atoms[:n]))
        n_size = mol.atoms[n].n

        total = 0
        for step in mol.cycle:
            for nodestate in step[n_pre + 1 : n_pre + n_size + 1]:
                if nodestate:
                    total += 1
                else:
                    total -= 1
        return total

    def _bond(self, mol_a, mol_b, n_a, n_b, bondsite_a, bondsite_b):
        newatoms = mol_a.atoms + mol_b.atoms
        offset = len(mol_a.atoms)
        # add offset to existing bonds in mol_b
        newbonds = dict(
            itertools.chain(
                mol_a.bonds.items(),
                (
                    ((n_a + offset, n_b + offset), ((i_a, j_a), (i_b, j_b)))
                    for ((n_a, n_b), ((i_a, j_a), (i_b, j_b))) in mol_b.bonds.items()
                ),
            )
        )
        # add the latest bond
        i_a, j_a = bondsite_a
        i_b, j_b = bondsite_b
        newbonds[(n_a, n_b + offset)] = ((i_a, j_a), (i_b, j_b))

        return RBNMol(newatoms, newbonds)

    def stabilize(self, reactant):
        # atoms are always stable
        if len(reactant.atoms) == 1:
            return (reactant,)

        to_break = set()
        for bond in reactant.bonds:
            n_a, n_b = bond
            a_score = self._score(reactant, n_a)
            b_score = self._score(reactant, n_b)
            if abs(a_score - b_score) < 1:
                # bond stable
                pass
            else:
                to_break.add(bond)

        # if no bonds broke, nothing changed
        if len(to_break) == 0:
            return (reactant,)

        # as bonds have broken, work out what connected components remain
        components = set(frozenset([i]) for i in range(len(reactant.atoms)))
        for bond in reactant.bonds:
            if bond not in to_break:
                n_a, n_b = bond
                component_a = next((x for x in components if n_a in x))
                component_b = next((x for x in components if n_b in x))
                components.remove(component_a)
                components.remove(component_b)
                components.add(component_a | component_b)  # union

        # now turn those components into molecules
        products = []
        for component in sorted(components):
            component_sorted = sorted(component)
            atoms = [reactant.atoms[i] for i in component_sorted]
            bonds = {}
            for key, value in reactant.bonds.items():
                if key in to_break:
                    continue
                atom_a, atom_b = key
                if atom_a in component and atom_b in component:
                    atom_a_new = component_sorted.index(atom_a)
                    atom_b_new = component_sorted.index(atom_b)
                    bonds[(atom_a_new, atom_b_new)] = value
            product = RBNMol(atoms, bonds)
            # recursively check for stability
            products.extend(self.stabilize(product))

        return tuple(products)

    def _get_possible_bondsites(self, mol, atom_n):
        """
        Returns tuple of tuples of possible bondsites on a specific atom
        Each bond site is i,j tuple where it is the jth
        input to the ith node of the nth atom
        """
        used_bondsites = set()
        for ((n_a, n_b), ((i_a, j_a), (i_b, j_b))) in mol.bonds.items():
            if n_a == atom_n:
                used_bondsites.add((i_a, j_a))
            if n_b == atom_n:
                used_bondsites.add((i_b, j_b))

        bondsites = []
        for i in range(mol.atoms[atom_n].n):
            for j in range(mol.atoms[atom_n].k):
                # only replace self-bonded inputs
                if mol.atoms[atom_n].inputs[i][j] == i:
                    bondsite = (i, j)
                    if bondsite not in used_bondsites:
                        bondsites.append(bondsite)
        return tuple(bondsites)

    def react(self, reactants, rng):
        """
        Enable simulation of reaction network.

        Elastic reactions return the reactants as the products.
        """

        assert len(reactants) == 2

        # swap reactants half the time
        # TODO this might not be necessary
        reactants = tuple(reactants)
        if rng.choice([0, 1]):
            mol_a = reactants[0]
            mol_b = reactants[1]
        else:
            mol_a = reactants[1]
            mol_b = reactants[0]

        # pick random atoms in each molecule to collide
        atom_a = rng.choice(range(len(mol_a.atoms)))
        atom_b = rng.choice(range(len(mol_b.atoms)))

        # work out what bond to add
        # TODO consider using this in score?
        bondsites_a = self._get_possible_bondsites(mol_a, atom_a)
        bondsites_b = self._get_possible_bondsites(mol_b, atom_b)

        if not len(bondsites_a) or not len(bondsites_b):
            # at least one reactant cant form new bonds, elastic collision
            return None

        bondsite_a = rng.choice(bondsites_a)
        bondsite_b = rng.choice(bondsites_b)

        a_score = self._score(mol_a, atom_a)
        b_score = self._score(mol_b, atom_b)

        # if scores are close enough, bond
        # TODO make this "temperature" based
        if abs(a_score - b_score) < 1:
            product = self._bond(mol_a, mol_b, atom_a, atom_b, bondsite_a, bondsite_b)
            if not product:
                # unable to bond, elastic collision
                return None
            else:
                # return product inside a tuple
                # TODO stabilization & decomposition
                return Reaction(reactants, (product,))
        else:
            # no bond, elastic collision
            return None

    def all_reactions(self, reactants):
        """
        May be implemented by subclasses to enable enumeration
        of reaction network.
        """
        assert len(reactants) == 2
        mol_a, mol_b = reactants

        reactions = collections.Counter()

        for atom_a in range(len(mol_a.atoms)):
            bondsites_a = self._get_possible_bondsites(mol_a, atom_a)
            for atom_b in range(len(mol_b.atoms)):
                bondsites_b = self._get_possible_bondsites(mol_b, atom_b)
                for bondsite_a in bondsites_a:
                    for bondsite_b in bondsites_b:
                        product = self._bond(
                            mol_a, mol_b, atom_a, atom_b, bondsite_a, bondsite_b
                        )
                        # TODO stabilization & decomposition
                        reaction = Reaction(reactants, (product,))
                        reactions[reaction] += 1

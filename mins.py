#!/usr/bin/env python3
import numpy as np
from ase import Atoms, Atom
from ase.calculators.calculator import Calculator, all_changes
from ase.neighborlist import NewPrimitiveNeighborList
from ase.build import make_supercell
import ase
from ase.neighborlist import NeighborList

class SmallBasisCorrection(Calculator):
    implemented_properties = ['energy', 'forces']
    implemented_properties += ['stress', 'stresses']

    def __init__(self):
        """SmallBasisCorrection
        """
        self.s = {'O': 0.108, 'V': 0.051, 'Nb': 0.008, 'Mo': 0.024}
        self.gamma = {'O': 2.489, 'V': 3.847, 'Nb': 4.1, 'Mo': 3.85}
        Calculator.__init__(self)

    def getS(self):
        return self.s

    def getGamma(self):
        return self.gamma

    def setS(self, s):
        self.s = s

    def setGamma(self, gamma):
        self.gamma = gamma

    def fswitch(self, RAB, RAB_uncorr_min):
        val = 1

        w = 0.5
        if RAB > RAB_uncorr_min:
            val = 0
        if (RAB_uncorr_min - w <= RAB and RAB <= RAB_uncorr_min):
            val = (np.cos(np.pi * (RAB - RAB_uncorr_min + w) / w) + 1) / 2
        return val

    def dfswitch(self, RAB, RAB_uncorr_min):
        val = 0
        w = 0.5
        if (RAB_uncorr_min - w <= RAB and RAB <= RAB_uncorr_min):
            val = (-np.pi / w * np.sin(np.pi * (RAB - RAB_uncorr_min + w) / w)) / 2
        return val

    def calculate(self, atoms=None,
                  properties=None,
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        energy = 0.0
        n_atom = atoms.get_global_number_of_atoms()
        forces = np.zeros((n_atom, 3))
        stresses = np.zeros((n_atom, 3, 3))

        p_ele = list((atoms.get_chemical_symbols()))
        RAB_uncorr_min = 0
        for ele in p_ele:
            if ele in ["O", "V", "Nb", "Mo"]:
                RAB_cut_tmp = self.getGamma()[ele]
                if RAB_cut_tmp > RAB_uncorr_min:
                    RAB_uncorr_min = RAB_cut_tmp
        p_cell = atoms
        p_ele = list((atoms.get_chemical_symbols()))
        symb_p = atoms.get_chemical_symbols()
        s_p = np.array(symb_p).copy()
        g_p = np.array(symb_p).copy()
        s_po = atoms.get_initial_magnetic_moments()
        g_po = atoms.get_initial_magnetic_moments()
        for ele in set(p_ele):
            mask = (s_p == ele)
            if ele in ["O", "V", "Nb", "Mo"]:
                s_po[mask] = (float(self.getS()[ele]))
                g_po[mask] = (float(1 / 2 * self.getGamma()[ele]))
        cutoffs = []
        for ele in p_ele:
            if ele in self.getGamma():
                cutoffs.append(self.getGamma()[ele])
            else:
                cutoffs.append(0)
        nl = NeighborList(cutoffs, skin=0, self_interaction=False, bothways=True, primitive=NewPrimitiveNeighborList)
        nl.update(atoms)

        for a in range(len(p_cell)):
            indices, offsets = nl.get_neighbors(a)  # neighbour indices and offsets
            symb = atoms[indices].get_chemical_symbols()
            Dj = np.add(atoms.positions[indices], np.dot(offsets, atoms.get_cell()))

            D = np.asarray(np.subtract(Dj, atoms.positions[a]))
            d = np.asarray(np.linalg.norm(D, axis=1))
            el_i = atoms[a].symbol
            sj = s_po[indices]
            gj = g_po[indices]
            s = np.add(sj, s_po[a])
            g = np.array(np.add(gj, g_po[a]))
            e = 0
            mask = (d <= g) & (d != 0.)
            fswitch = d.copy()
            dfswitch = d.copy()
            mask1 = (d <= np.subtract(g, 0.5))
            mask2 = (d >= (np.subtract(g, 0.5))) & (d <= g)
            mask3 = (d > g) | (d == 0.)
            fswitch[mask2] = (np.cos(np.pi * (np.subtract(d, g)[mask2] + 0.5) / 0.5) + 1) / 2
            dfswitch[mask2] = (-np.pi / 0.5 * np.sin(np.pi * (np.subtract(d, g)[mask2] + 0.5) / 0.5)) / 2
            fswitch[mask1] = 1.
            dfswitch[mask1] = 0.
            fswitch[mask3] = 0.
            dfswitch[mask3] = 0.
            e = np.multiply(np.multiply(s, np.subtract(d, g)), fswitch)
            pairwise_forces = (np.asarray(np.multiply(s, fswitch))) + (
                np.multiply(np.multiply(s, np.subtract(d, g)), dfswitch))
            pairwise_forces = pairwise_forces[:, np.newaxis] * np.divide(D, d[:, None])
            forces[a, :] += pairwise_forces.sum(axis=0)
            energy += e.sum() * 0.5
            stresses[a] += 0.5 * np.dot(
                pairwise_forces.T, D
            )  # equivalent to outer product

        # no lattice, no stress
        # if self.atoms.cell.rank == 3:
            # stresses = full_3x3_to_voigt_6_stress(stresses)
            # self.results['stress'] = stresses.sum(axis=0) / self.atoms.get_volume()
            # self.results['stresses'] = stresses / self.atoms.get_volume()
        self.results['forces'] = forces
        self.results['energy'] = energy


class Element:

    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__nucleus_dic = {'H': '1', 'He': '2', 'Li': '3', 'Be': '4', 'B': '5', 'C': '6', 'N': '7', 'O': '8',
                              'F': '9', 'Ne': '10', 'Na': '11', 'Mg': '12', 'Al': '13', 'Si': '14', 'P': '15',
                              'S': '16', 'Cl': '17', 'Ar': '18', 'K': '19', 'Ca': '20', 'Sc': '21', 'Ti': '22',
                              'V': '23', 'Cr': '24', 'Mn': '25', 'Fe': '26', 'Co': '27', 'Ni': '28', 'Cu': '29',
                              'Zn': '30', 'Ga': '31', 'Ge': '32', 'As': '33', 'Se': '34', 'Br': '35', 'Kr': '36',
                              'Rb': '37', 'Sr': '38', 'Y': '39', 'Zr': '40', 'Nb': '41', 'Mo': '42', 'Tc': '43',
                              'Ru': '44', 'Rh': '45', 'Pd': '46', 'Ag': '47', 'Cd': '48', 'In': '49', 'Sn': '50',
                              'Sb': '51', 'Te': '52', 'I': '53', 'Xe': '54', 'Cs': '55', 'Ba': '56', 'La': '57',
                              'Ce': '58', 'Pr': '59', 'Nd': '60', 'Pm': '61', 'Sm': '62', 'Eu': '63', 'Gd': '64',
                              'Tb': '65', 'Dy': '66', 'Ho': '67', 'Er': '68', 'Tm': '69', 'Yb': '70', 'Lu': '71',
                              'Hf': '72', 'Ta': '73', 'W': '74', 'Re': '75', 'Os': '76', 'Ir': '77', 'Pt': '78',
                              'Au': '79', 'Hg': '80', 'Tl': '81', 'Pb': '82', 'Bi': '83', 'Po': '84', 'At': '85',
                              'Rn': '86', 'Fr': '87', 'Ra': '88', 'Ac': '89', 'Th': '90', 'Pa': '91', 'U': '92',
                              'Np': '93', 'Pu': '94', 'Am': '95', 'Cm': '96', 'Bk': '97', 'Cf': '98', 'Es': '99',
                              'Fm': '100', 'Md': '101', 'No': '102'}
        self.__nucleus = self.__nucleus_dic[name_species]
        self.__mass_dic = {'H': '1.00794', 'He': '4.002602', 'Li': '6.941', 'Be': '9.012182', 'B': '10.811',
                           'C': '12.0107', 'N': '14.0067', 'O': '15.9994', 'F': '18.9984032', 'Ne': '20.1797',
                           'Na': '22.9897692', 'Mg': '24.3050', 'Al': '26.9815386', 'Si': '28.0855', 'P': '30.973762',
                           'S': '32.065', 'Cl': '35.453', 'Ar': '39.948', 'K': '39.0983', 'Ca': '40.078',
                           'Sc': '44.955912', 'Ti': '47.867', 'V': '50.9415', 'Cr': '51.9961', 'Mn': '54.938045',
                           'Fe': '55.845', 'Co': '58.933', 'Ni': '58.6934', 'Cu': '63.546', 'Zn': '65.409',
                           'Ga': '69.723', 'Ge': '72.64', 'As': '74.92160', 'Se': '78.96', 'Br': '79.904',
                           'Kr': '83.798', 'Rb': '85.467', 'Sr': '87.62', 'Y': '88.905', 'Zr': '91.224',
                           'Nb': '92.90638', 'Mo': '95.94', 'Tc': '98', 'Ru': '101.07', 'Rh': '102.90550',
                           'Pd': '106.42', 'Ag': '107.8682', 'Cd': '112.411', 'In': '114.818', 'Sn': '118.710',
                           'Sb': '121.760', 'Te': '127.60', 'I': '126.90', 'Xe': '131.29', 'Cs': '132.90',
                           'Ba': '137.32', 'La': '138.90', 'Ce': '140.11', 'Pr': '140.90', 'Nd': '144.24', 'Pm': '145',
                           'Sm': '150.36', 'Eu': '151.96', 'Gd': '157.25', 'Tb': '158.92', 'Dy': '162.50',
                           'Ho': '164.93', 'Er': '167.25', 'Tm': '168.93', 'Yb': '173.04', 'Lu': '174.96',
                           'Hf': '178.49', 'Ta': '180.94', 'W': '183.84', 'Re': '186.20', 'Os': '190.23',
                           'Ir': '192.217', 'Pt': '195.084', 'Au': '196.966569', 'Hg': '200.59', 'Tl': '204.3833',
                           'Pb': '207.2', 'Bi': '208.98040', 'Po': '20', 'At': '21', 'Rn': '22', 'Fr': '22', 'Ra': '22',
                           'Ac': '22', 'Th': '232.03806', 'Pa': '231.03588', 'U': '238.02891', 'Np': '23', 'Pu': '24',
                           'Am': '24', 'Cm': '24', 'Bk': '24', 'Cf': '25', 'Es': '25', 'Fm': '25', 'Md': '25',
                           'No': '25'}
        self.__mass = self.__mass_dic[name_species]
        self.__atomic_number_dic = {'H': '1', 'He': '2', 'Li': '3', 'Be': '4', 'B': '5', 'C': '6', 'N': '7', 'O': '8',
                                    'F': '9', 'Ne': '10', 'Na': '11', 'Mg': '12', 'Al': '13', 'Si': '14', 'P': '15',
                                    'S': '16', 'Cl': '17', 'Ar': '18', 'K': '19', 'Ca': '20', 'Sc': '21', 'Ti': '22',
                                    'V': '23', 'Cr': '24', 'Mn': '25', 'Fe': '26', 'Co': '27', 'Ni': '28', 'Cu': '29',
                                    'Zn': '30', 'Ga': '31', 'Ge': '32', 'As': '33', 'Se': '34', 'Br': '35', 'Kr': '36',
                                    'Rb': '37', 'Sr': '38', 'Y': '39', 'Zr': '40', 'Nb': '41', 'Mo': '42', 'Tc': '43',
                                    'Ru': '44', 'Rh': '45', 'Pd': '46', 'Ag': '47', 'Cd': '48', 'In': '49', 'Sn': '50',
                                    'Sb': '51', 'Te': '52', 'I': '53', 'Xe': '54', 'Cs': '55', 'Ba': '56', 'La': '57',
                                    'Ce': '58', 'Pr': '59', 'Nd': '60', 'Pm': '61', 'Sm': '62', 'Eu': '63', 'Gd': '64',
                                    'Tb': '65', 'Dy': '66', 'Ho': '67', 'Er': '68', 'Tm': '69', 'Yb': '70', 'Lu': '71',
                                    'Hf': '72', 'Ta': '73', 'W': '74', 'Re': '75', 'Os': '76', 'Ir': '77', 'Pt': '78',
                                    'Au': '79', 'Hg': '80', 'Tl': '81', 'Pb': '82', 'Bi': '83', 'Po': '84', 'At': '85',
                                    'Rn': '86', 'Fr': '87', 'Ra': '88', 'Ac': '89', 'Th': '90', 'Pa': '91', 'U': '92',
                                    'Np': '93', 'Pu': '94', 'Am': '95', 'Cm': '96', 'Bk': '97', 'Cf': '98', 'Es': '99',
                                    'Fm': '100', 'Md': '101', 'No': '102'}
        self.__atomic_number = self.__atomic_number_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getNucleus(self):
        return self.__nucleus

    def getMass(self):
        return self.__mass

    def getAtomicNumber(self):
        return self.__atomic_number


class Lhartree:

    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__lhartree_dic = {'H': '4', 'He': '4', 'Li': '4', 'Be': '4', 'B': '4', 'C': '4', 'N': '4', 'O': '4',
                               'F': '4', 'Ne': '4', 'Na': '4', 'Mg': '4', 'Al': '4', 'Si': '4', 'P': '4', 'S': '4',
                               'Cl': '4', 'Ar': '4', 'K': '4', 'Ca': '4', 'Sc': '4', 'Ti': '4', 'V': '4', 'Cr': '4',
                               'Mn': '4', 'Fe': '4', 'Co': '4', 'Ni': '4', 'Cu': '4', 'Zn': '4', 'Ga': '4', 'Ge': '4',
                               'As': '4', 'Se': '4', 'Br': '4', 'Kr': '4', 'Rb': '4', 'Sr': '4', 'Y': '4', 'Zr': '4',
                               'Nb': '4', 'Mo': '4', 'Tc': '4', 'Ru': '4', 'Rh': '4', 'Pd': '4', 'Ag': '4', 'Cd': '4',
                               'In': '4', 'Sn': '4', 'Sb': '4', 'Te': '4', 'I': '4', 'Xe': '4', 'Cs': '4', 'Ba': '4',
                               'La': '4', 'Ce': '4', 'Pr': '4', 'Nd': '4', 'Pm': '4', 'Sm': '4', 'Eu': '4', 'Gd': '4',
                               'Tb': '4', 'Dy': '4', 'Ho': '4', 'Er': '4', 'Tm': '4', 'Yb': '4', 'Lu': '4', 'Hf': '4',
                               'Ta': '4', 'W': '4', 'Re': '4', 'Os': '4', 'Ir': '4', 'Pt': '4', 'Au': '4', 'Hg': '4',
                               'Tl': '4', 'Pb': '4', 'Bi': '4', 'Po': '4', 'At': '4', 'Rn': '4', 'Fr': '4', 'Ra': '4',
                               'Ac': '4', 'Th': '4', 'Pa': '4', 'U': '4', 'Np': '4', 'Pu': '4', 'Am': '4', 'Cm': '4',
                               'Bk': '4', 'Cf': '4', 'Es': '4', 'Fm': '4', 'Md': '4', 'No': '4'}
        self.__lhartree = self.__lhartree_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getLhartree(self):
        return self.__lhartree


class CutoffPotential:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__cut_pot_dic = {'H': '3.5  1.5  1.0', 'He': '3.5  1.5  1.0', 'Li': '3.5  1.5  1.0', 'Be': '3.5  1.5  1.0',
                              'B': '3.5  1.5  1.0', 'C': '3.5  1.5  1.0', 'N': '3.5  1.5  1.0', 'O': '3.5  1.5  1.0',
                              'F': '3.5  1.5  1.0', 'Ne': '3.5  1.5  1.0', 'Na': '4.0          1.5  1.0',
                              'Mg': '4.0          1.5  1.0', 'Al': '3.5          1.5  1.0',
                              'Si': '3.5          1.5  1.0', 'P': '3.5          1.5  1.0', 'S': '3.5          1.5  1.0',
                              'Cl': '3.5          1.5  1.0', 'Ar': '3.5          1.5  1.0',
                              'K': '4.0          1.5  1.0', 'Ca': '4.0          1.5  1.0',
                              'Sc': '3.5          1.5  1.0', 'Ti': '3.5          1.5  1.0',
                              'V': '3.5          1.5  1.0', 'Cr': '3.5          1.5  1.0',
                              'Mn': '3.5          1.5  1.0', 'Fe': '3.5          1.5  1.0', 'Co': '3.5  1.5  1.0',
                              'Ni': '3.5          1.5  1.0', 'Cu': '3.5  1.5  1.0', 'Zn': '3.5          1.5  1.0',
                              'Ga': '3.5          1.5  1.0', 'Ge': '3.5          1.5  1.0',
                              'As': '3.5          1.5  1.0', 'Se': '3.5          1.5  1.0',
                              'Br': '3.5          1.5  1.0', 'Kr': '3.5          1.5  1.0', 'Rb': '4.0  1.5  1.0',
                              'Sr': '4.0  1.5  1.0', 'Y': '4.0  1.5  1.0', 'Zr': '3.5  1.5  1.0', 'Nb': '3.5  1.5  1.0',
                              'Mo': '3.5  1.5  1.0', 'Tc': '3.5  1.5  1.0', 'Ru': '3.5  1.5  1.0',
                              'Rh': '3.5  1.5  1.0', 'Pd': '3.5  1.5  1.0', 'Ag': '3.5  1.5  1.0',
                              'Cd': '3.5  1.5  1.0', 'In': '3.5  1.5  1.0', 'Sn': '3.5  1.5  1.0',
                              'Sb': '3.5  1.5  1.0', 'Te': '3.5  1.5  1.0', 'I': '3.5  1.5  1.0', 'Xe': '3.5  1.5  1.0',
                              'Cs': '4.0  1.5  1.0', 'Ba': '4.0  1.5  1.0', 'La': '4.0  1.5  1.0',
                              'Ce': '4.0  1.5  1.0', 'Pr': '4.0  1.5  1.0', 'Nd': '4.0  1.5  1.0',
                              'Pm': '4.0  1.5  1.0', 'Sm': '4.0  1.5  1.0', 'Eu': '4.0  1.5  1.0',
                              'Gd': '4.0  1.5  1.0', 'Tb': '4.0  1.5  1.0', 'Dy': '4.0  1.5  1.0',
                              'Ho': '4.0  1.5  1.0', 'Er': '4.0  1.5  1.0', 'Tm': '4.0  1.5  1.0',
                              'Yb': '4.0  1.5  1.0', 'Lu': '4.0  1.5  1.0', 'Hf': '3.5  1.5  1.0',
                              'Ta': '3.5  1.5  1.0', 'W': '3.5  1.5  1.0', 'Re': '3.5  1.5  1.0', 'Os': '3.5  1.5  1.0',
                              'Ir': '3.5  1.5  1.0', 'Pt': '3.5  1.5  1.0', 'Au': '3.5  1.5  1.0',
                              'Hg': '3.5  1.5  1.0', 'Tl': '4.0  1.5  1.0', 'Pb': '4.0  1.5  1.0',
                              'Bi': '4.0  1.5  1.0', 'Po': '4.0  1.5  1.0', 'At': '4.0  1.5  1.0',
                              'Rn': '4.0  1.5  1.0', 'Fr': '4.0  1.5  1.0', 'Ra': '4.0  1.5  1.0',
                              'Ac': '4.0  1.5  1.0', 'Th': '4.0  1.5  1.0', 'Pa': '4.0  1.5  1.0', 'U': '4.0  1.5  1.0',
                              'Np': '4.0  1.5  1.0', 'Pu': '4.0  1.5  1.0', 'Am': '4.0  1.5  1.0',
                              'Cm': '4.0  1.5  1.0', 'Bk': '4.0  1.5  1.0', 'Cf': '4.0  1.5  1.0',
                              'Es': '4.0  1.5  1.0', 'Fm': '4.0  1.5  1.0', 'Md': '4.0  1.5  1.0',
                              'No': '4.0  1.5  1.0'}
        self.__cut_pot = self.__cut_pot_dic[name_species]
        self.__basis_dep_cutoff_dic = {'H': '1e-4', 'He': '1e-4', 'Li': '1e-4', 'Be': '1e-4', 'B': '1e-4', 'C': '1e-4',
                                       'N': '1e-4', 'O': '1e-4', 'F': '1e-4', 'Ne': '1e-4', 'Na': '1e-4', 'Mg': '1e-4',
                                       'Al': '1e-4', 'Si': '1e-4', 'P': '1e-4', 'S': '1e-4', 'Cl': '1e-4', 'Ar': '1e-4',
                                       'K': '1e-4', 'Ca': '1e-4', 'Sc': '1e-4', 'Ti': '1e-4', 'V': '1e-4', 'Cr': '1e-4',
                                       'Mn': '1e-4', 'Fe': '1e-4', 'Co': '1e-4', 'Ni': '1e-4', 'Cu': '1e-4',
                                       'Zn': '1e-4', 'Ga': '1e-4', 'Ge': '1e-4', 'As': '1e-4', 'Se': '1e-4',
                                       'Br': '1e-4', 'Kr': '1e-4', 'Rb': '1e-4', 'Sr': '1e-4', 'Y': '1e-4',
                                       'Zr': '1e-4', 'Nb': '1e-4', 'Mo': '1e-4', 'Tc': '1e-4', 'Ru': '1e-4',
                                       'Rh': '1e-4', 'Pd': '1e-4', 'Ag': '1e-4', 'Cd': '1e-4', 'In': '1e-4',
                                       'Sn': '1e-4', 'Sb': '1e-4', 'Te': '1e-4', 'I': '1e-4', 'Xe': '1e-4',
                                       'Cs': '1e-4', 'Ba': '1e-4', 'La': '1e-4', 'Ce': '1e-4', 'Pr': '1e-4',
                                       'Nd': '1e-4', 'Pm': '1e-4', 'Sm': '1e-4', 'Eu': '1e-4', 'Gd': '1e-4',
                                       'Tb': '1e-4', 'Dy': '1e-4', 'Ho': '1e-4', 'Er': '1e-4', 'Tm': '1e-4',
                                       'Yb': '1e-4', 'Lu': '1e-4', 'Hf': '1e-4', 'Ta': '1e-4', 'W': '1e-4',
                                       'Re': '1e-4', 'Os': '1e-4', 'Ir': '1e-4', 'Pt': '1e-4', 'Au': '1e-4',
                                       'Hg': '1e-4', 'Tl': '1e-4', 'Pb': '1e-4', 'Bi': '1e-4', 'Po': '1e-4',
                                       'At': '1e-4', 'Rn': '1e-4', 'Fr': '1e-4', 'Ra': '1e-4', 'Ac': '1e-4',
                                       'Th': '1e-4', 'Pa': '1e-4', 'U': '1e-4', 'Np': '1e-4', 'Pu': '1e-4',
                                       'Am': '1e-4', 'Cm': '1e-4', 'Bk': '1e-4', 'Cf': '1e-4', 'Es': '1e-4',
                                       'Fm': '1e-4', 'Md': '1e-4', 'No': '1e-4'}
        self.__basis_dep_cutoff = self.__basis_dep_cutoff_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getCutPot(self):
        return self.__cut_pot

    def getBasisDepCutoff(self):
        return self.__basis_dep_cutoff


class IntegrationGrid:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__radial_base_dic = {'H': '24 5.0', 'He': '27 5.0', 'Li': '29 5.0', 'Be': '31 5.0', 'B': '32 5.0',
                                  'C': '34 5.0', 'N': '35 5.0', 'O': '36 5.0', 'F': '37 5.0', 'Ne': '38 5.0',
                                  'Na': '40 5.5', 'Mg': '40 5.5', 'Al': '41 5.0', 'Si': '42 5.0', 'P': '43 5.0',
                                  'S': '44 5.0', 'Cl': '45 5.0', 'Ar': '46 5.0', 'K': '46 5.5', 'Ca': '47 5.5',
                                  'Sc': '47 5.0', 'Ti': '48 5.0', 'V': '49 5.0', 'Cr': '50 5.0', 'Mn': '50 5.0',
                                  'Fe': '51 5.0', 'Co': '52 5.0', 'Ni': '52 5.0', 'Cu': '53 5.0', 'Zn': '53 5.0',
                                  'Ga': '54 5.0', 'Ge': '54 5.0', 'As': '55 5.0', 'Se': '55 5.0', 'Br': '56 5.0',
                                  'Kr': '56 5.0', 'Rb': '57  5.5', 'Sr': '57  5.5', 'Y': '58  5.5', 'Zr': '58  5.0',
                                  'Nb': '59 5.0', 'Mo': '59 5.0', 'Tc': '60 5.0', 'Ru': '60 5.0', 'Rh': '61 5.0',
                                  'Pd': '62 5.0', 'Ag': '62 5.0', 'Cd': '62 5.0', 'In': '62 5.0', 'Sn': '63 5.0',
                                  'Sb': '63 5.0', 'Te': '64 5.0', 'I': '64  5.0', 'Xe': '64  5.0', 'Cs': '65  5.5',
                                  'Ba': '65  5.5', 'La': '65  5.5', 'Ce': '66  5.5', 'Pr': '66  5.5', 'Nd': '66  5.5',
                                  'Pm': '67  5.5', 'Sm': '67  5.5', 'Eu': '68  5.5', 'Gd': '68  5.5', 'Tb': '68  5.5',
                                  'Dy': '69  5.5', 'Ho': '69  5.5', 'Er': '69  5.5', 'Tm': '70  5.5', 'Yb': '70  5.5',
                                  'Lu': '70  5.5', 'Hf': '71  5.0', 'Ta': '71  5.0', 'W': '71  5.0', 'Re': '72  5.0',
                                  'Os': '72  5.0', 'Ir': '72 5.0', 'Pt': '72 5.0', 'Au': '73 5.0', 'Hg': '73 5.0',
                                  'Tl': '73 5.5', 'Pb': '74 5.5', 'Bi': '74 5.5', 'Po': '74 5.5', 'At': '74 5.5',
                                  'Rn': '75 5.5', 'Fr': '75 5.5', 'Ra': '75 5.5', 'Ac': '76 5.5', 'Th': '76 5.5',
                                  'Pa': '76 5.5', 'U': '76 5.5', 'Np': '77 5.5', 'Pu': '77 5.5', 'Am': '77 5.5',
                                  'Cm': '77 5.5', 'Bk': '78 5.5', 'Cf': '78 5.5', 'Es': '78 5.5', 'Fm': '78 5.5',
                                  'Md': '79 5.5', 'No': '79 5.5'}
        self.__radial_base = self.__radial_base_dic[name_species]
        self.__angular_grids_dic = {
            'H': 'specified\n      division   0.2421   50\n      division   0.3822  110\n      division   0.4799  194\n      division   0.5341  302\n#      division   0.5626  434\n#      division   0.5922  590\n#      division   0.6542  770\n#      division   0.6868 1202\n#      outer_grid  770\n      outer_grid  302\n',
            'He': 'specified\n      division   0.3349  110\n      division   0.4719  194\n      division   0.5352  302\n#      division   1.8809  770\n#      outer_grid    770\n      outer_grid    302\n',
            'Li': 'specified\n      division   0.4484  110\n      division   0.5659  194\n      division   0.6315  302\n#      division   0.6662  434\n#      division   0.8186  590\n#      division   0.9037  770\n#      division   6.2760  974\n#      outer_grid   974\n      outer_grid   302\n',
            'Be': 'specified\n      division   0.4283  110\n      division   0.4792  194\n      division   0.5061  302\n#      division   0.7227  434\n#      division   0.8724  590\n#      division   0.9555  770\n#      division   2.9770  974\n#      outer_grid   974\n      outer_grid 302\n',
            'B': 'specified\n      division   0.3742  110\n      division   0.5197  194\n      division   0.5753  302\n#      division   0.7664  434\n#      division   0.8392  770\n#      division   1.6522  974\n#      outer_grid   974\n      outer_grid   302\n',
            'C': 'specified\n      division   0.3326   50\n      division   0.5710  110\n      division   0.7727  194\n      division   0.8772  302\n#      division   0.9334  434\n#      division   0.9625  590\n#      division   0.9924  770\n#      division   1.0230  974\n#      division   1.4589 1202\n',
            'N': 'specified\n      division   0.2599   50\n      division   0.4601  110\n      division   0.5885  194\n      division   0.6503  302\n#      division   0.6939  434\n#      division   0.7396  590\n#      division   0.7632  770\n#      division   0.8122  974\n#      division   1.1604 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'O': 'specified\n      division   0.2659   50\n      division   0.4451  110\n      division   0.6052  194\n      division   0.7543  302\n#      division   0.8014  434\n#      division   0.8507  590\n#      division   0.8762  770\n#      division   0.9023  974\n#      division   1.2339 1202\n#      outer_grid 974\n      outer_grid 302\n',
            'F': 'specified\n      division   0.4014  110\n      division   0.5291  194\n      division   0.6019  302\n#      division   0.6814  434\n#      division   0.7989  590\n#      division   0.8965  770\n#      division   1.3427  974\n#      outer_grid   974\n      outer_grid   302\n',
            'Ne': 'specified\n      division   0.4737  110\n      division   0.6363  194\n      division   0.7448  302\n#      division   0.8348  590\n#      division   1.0034  770\n#      division   1.7611  974\n#      outer_grid   974\n      outer_grid   302\n',
            'Na': 'specified\n      division   0.5925  110\n      division   0.7843  194\n      division   1.0201  302\n#      division   1.1879  434\n#      division   1.3799  590\n#      division   1.4503  770\n#      division   7.0005  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Mg': 'specified\n      division   0.7029   50\n      division   0.9689  110\n      division   1.1879  194\n      division   1.3129  302\n#      division   1.4867  434\n#      division   1.6018  590\n#      division   1.8611  770\n#      division   1.9576  974\n#      division   2.2261 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Al': 'specified\n      division   0.6594  110\n      division   0.8170  194\n      division   0.9059  302\n#      division   1.0363  434\n#      division   1.1443  590\n#      division   1.2621  770\n#      division   2.8177  974\n#      outer_grid   974\n      outer_grid   302\n',
            'Si': 'specified\n      division   0.5866   50\n      division   0.9616  110\n      division   1.2249  194\n      division   1.3795  302\n#      division   1.4810  434\n#      division   1.5529  590\n#      division   1.6284  770\n#      division   1.7077  974\n#      division   2.4068 1202\n#      outer_grid   974\n      outer_grid 302\n',
            'P': 'specified\n      division   0.5527   50\n      division   0.8801  110\n      division   1.1447  194\n      division   1.3165  302\n#      division   1.4113  434\n#      division   1.4781  590\n#      division   1.5482  770\n#      division   1.5845  974\n#      division   2.2606 1202\n      outer_grid  302\n',
            'S': 'specified\n      division   0.4665  110\n      division   0.5810  194\n      division   0.7139  302\n#      division   0.8274  434\n#      division   0.9105  590\n#      division   1.0975  770\n#      division   1.2028  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Cl': 'specified\n      division   0.4412  110\n      division   0.5489  194\n      division   0.6734  302\n#      division   0.7794  434\n#      division   0.9402  590\n#      division   1.0779  770\n#      division   1.1792  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Ar': 'specified\n      division   0.5855  110\n      division   0.8590  194\n      division   0.9692  302\n#      division   1.1235  590\n#      division   1.1911  770\n#      division   1.2623  974\n#      outer_grid  974\n      outer_grid  302\n',
            'K': 'specified\n      division   0.5285  110\n      division   0.7831  194\n      division   0.9986  302\n#      division   1.0594  434\n#      division   1.1569  590\n#      division   1.2994  770\n#      division   1.4587  974\n#      outer_grid 974\n      outer_grid 302\n',
            'Ca': 'specified\n',
            'Sc': 'specified\n      division   0.6021   50\n      division   1.1116  110\n      division   1.4663  194\n      division   1.6660  302\n#      division   1.8551  434\n#      division   2.0245  590\n#      division   2.2132  770\n#      division   2.5421  974\n#      division   3.1021 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Ti': 'specified\n      division   0.5171   50\n      division   0.9824  110\n      division   1.2917  194\n      division   1.4940  302\n#      division   1.6934  434\n#      division   1.8425  590\n#      division   2.1901  770\n#      division   2.2896  974\n#      division   2.8244 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'V': 'specified\n      division   0.4553   50\n      division   0.8707  110\n      division   1.1666  194\n      division   1.3737  302\n#      division   1.5524  434\n#      division   1.8303  590\n#      division   2.0330  770\n#      division   2.0767  974\n#      division   2.5907 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Cr': 'specified\n      division   0.4331   50\n      division   0.8246  110\n      division   1.1008  194\n      division   1.3188  302\n#      division   1.4867  434\n#      division   1.7819  590\n#      division   1.9339  770\n#      division   1.9742  974\n#      division   2.4437 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Mn': 'specified\n      division   0.4222   50\n      division   0.8072  110\n      division   1.0787  194\n      division   1.2927  302\n#      division   1.4573  434\n#      division   1.8560  590\n#      division   1.8945  770\n#      division   1.9339  974\n#      division   2.3905 1202\n      outer_grid   302\n',
            'Fe': 'specified\n      division   0.4337   50\n      division   0.8154  110\n      division   1.1046  194\n      division   1.3713  302\n#      division   1.5424  434\n#      division   1.7365  590\n#      division   1.9587  770\n#      division   1.9990  974\n#      division   2.4139 1202\n#      outer_grid  1202\n      outer_grid  302\n',
            'Co': 'specified\n      division   0.4668   50\n      division   0.8401  110\n      division   1.1973  194\n      division   1.4237  302\n#      division   1.5981  434\n#      division   1.7961  590\n#      division   1.9829  770\n#      division   2.0231  974\n#      division   2.4367 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Ni': 'specified\n      division   0.4449   50\n      division   0.8232  110\n      division   1.1299  194\n      division   1.4513  302\n#      division   1.6613  434\n#      division   1.8317  590\n#      division   1.9829  770\n#      division   2.0231  974\n#      division   2.4367 1202\n      outer_grid  302\n',
            'Cu': 'specified\n      division   0.5231   50\n      division   0.8642  110\n      division   1.1767  194\n      division   1.5041  302\n#      division   1.9293  434\n#      division   2.0065  590\n#      division   2.0466  770\n#      division   2.0877  974\n#      division   2.4589 1202\n      outer_grid   302\n',
            'Zn': 'specified\n      division   0.5114   50\n      division   0.8989  110\n      division   1.2692  194\n      division   1.6226  302\n#      division   1.7854  434\n#      division   2.0877  590\n#      division   2.1298  770\n#      division   2.1730  974\n#      division   2.5659 1202\n      outer_grid  302\n',
            'Ga': 'specified\n      division   0.5103   50\n      division   0.8880  110\n      division   1.2009  194\n      division   1.5000  302\n#      division   1.7093  434\n#      division   1.8791  590\n#      division   1.9525  770\n#      division   2.3801 1202\n#      outer_grid  1202\n      outer_grid  302\n',
            'Ge': 'specified\n      division   0.0947  110\n      division   0.1314  194\n      division   0.7746  302\n#      division   0.8710  434\n#      division   0.9770  590\n#      division   1.1356  770\n#      division   2.6430  974\n#      outer_grid  974\n      outer_grid  302\n',
            'As': 'specified\n      division   0.4982   50\n      division   0.9113  110\n      division   1.1593  194\n      division   1.4959  302\n#      division   1.6697  434\n#      division   1.8319  590\n#      division   1.9752  770\n#      division   2.0131  974\n#      division   2.4015 1202\n#      outer_grid  1202\n      outer_grid  302\n',
            'Se': 'specified\n      division   0.0830  110\n      division   0.1357  194\n      division   0.7377  302\n#      division   0.8610  434\n#      division   0.9640  590\n#      division   1.0773  770\n#      division   2.5539  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Br': 'specified\n      division   0.0871  110\n      division   0.1400  194\n      division   0.7896  302\n#      division   0.8837  434\n#      division   0.9869  590\n#      division   1.0613  770\n#      division   2.6835  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Kr': 'specified\n      division   0.3980   50\n      division   0.7174  110\n      division   1.0235  194\n      division   1.1824  302\n#      division   1.3889  434\n#      division   1.8888  590\n#      division   1.9972  770\n#      division   2.1543  974\n#      division   2.4715 1202\n#      outer_grid  1202\n      outer_grid  302\n',
            'Rb': 'specified\n      division   0.1250  110\n      division   0.9394  194\n      division   1.1230  302\n#      division   1.2051  434\n#      division   1.2929  590\n#      division   1.3869  770\n#      division   7.0005  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Sr': 'specified\n      division   0.6981  110\n      division   0.9394  194\n      division   1.1230  302\n#      division   1.2482  434\n#      division   1.3391  590\n#      division   1.4365  770\n#      division   7.0005  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Y': 'specified\n      division   0.7193   50\n      division   1.2925  110\n      division   1.6473  194\n      division   1.8976  302\n#      division   2.1161  434\n#      division   2.4151  590\n#      division   2.7220  770\n#      division   2.7789  974\n#      division   3.4772 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Zr': 'specified\n      division   0.5825   50\n      division   1.1060  110\n      division   1.4586  194\n      division   1.7061  302\n#      division   1.9320  434\n#      division   2.2803  590\n#      division   2.4151  770\n#      division   2.4626  974\n#      division   3.1649 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Nb': 'specified\n      division   0.5255   50\n      division   0.9829  110\n      division   1.2922  194\n      division   1.6123  302\n#      division   1.9879  434\n#      division   2.0980  590\n#      division   2.1365  770\n#      division   2.1759  974\n#      division   2.8558 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Mo': 'specified\n      division   0.4847   50\n      division   0.9168  110\n      division   1.2280  194\n      division   1.6402  302\n#      division   1.8849  434\n#      division   2.0237  590\n#      division   2.0980  974\n#      division   2.7405 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Tc': 'specified\n      division   0.4840   50\n      division   0.9061  110\n      division   1.2286  194\n      division   1.6333  302\n#      division   1.9054  434\n#      division   1.9732  590\n#      division   2.0439  770\n#      division   2.0806  974\n#      division   2.7039 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Ru': 'specified\n      division   0.4743   50\n      division   0.8754  110\n      division   1.1882  194\n      division   1.6059  302\n#      division   1.8727  434\n#      division   1.9389  590\n#      division   2.0082  770\n#      division   2.0439  974\n#      division   2.6509 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Rh': 'specified\n      division   0.4929   50\n      division   0.8959  110\n      division   1.2091  194\n      division   1.6265  302\n#      division   1.8928  434\n#      division   1.9931  590\n#      division   2.0637  770\n#      division   2.6689 1202\n#      outer_grid   974\n      outer_grid   302\n',
            'Pd': 'specified\n      division   0.5211   50\n      division   0.9161  110\n      division   1.2296  194\n      division   1.5678  302\n#      division   1.9785  434\n#      division   2.0474  590\n#      division   2.1195  770\n#      division   2.1568  974\n#      division   2.7392 1202\n#       outer_grid  974\n       outer_grid  302\n',
            'Ag': 'specified\n      division   0.5617   50\n      division   0.9788  110\n      division   1.2700  194\n      division   1.5678  302\n#      division   2.1195  434\n#      division   2.1949  590\n#      division   2.2741  770\n#      division   2.8497 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Cd': 'specified\n      division   0.5937   50\n      division   1.0282  110\n      division   1.3769  194\n      division   1.7301  302\n#      division   2.2341  434\n#      division   2.2741  590\n#      division   2.3152  770\n#      division   2.3574  974\n#      division   2.9077 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'In': 'specified\n      division   0.1831  110\n      division   0.9161  194\n      division   1.0115  302\n#      division   1.1156  434\n#      division   1.1524  590\n#      division   1.2296  770\n#      division   7.0005  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Sn': 'specified\n      division   0.1666  110\n      division   0.9058  302\n#      division   0.9669  434\n#      division   1.0315  590\n#      division   1.0999  770\n#      division   3.0459  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Sb': 'specified\n      division   0.1144  110\n      division   0.1571  194\n      division   0.8765  302\n#      division   0.9669  434\n#      division   1.0315  590\n#      division   1.0999  770\n#      division   3.0459  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Te': 'specified\n      division   0.1259  110\n      division   0.8959  194\n      division   0.9864  302\n#      division   1.1196  434\n#      division   1.1922  590\n#      division   1.3098  770\n#      division   2.9404  974\n#      outer_grid  974\n      outer_grid  302\n',
            'I': 'specified\n      division   0.1103  110\n      division   0.1515  194\n      division   0.9554  302\n#      division   1.1196  590\n#      division   1.1922  770\n#      division   6.1948  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Xe': 'specified\n      division   0.7862  110\n      division   1.1196  194\n      division   1.4850  302\n#      division   1.6329  434\n#      division   1.6858  590\n#      division   1.7978  770\n#      division   1.9188  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Cs': 'specified\n      division   0.7542  110\n      division   1.0056  194\n      division   1.2887  302\n#      division   1.4138  434\n#      division   1.5042  590\n#      division   1.6519  770\n#      outer_grid  974\n      outer_grid  302\n',
            'Ba': 'specified\n      division   0.6752  110\n      division   0.9746  194\n      division   1.2241  302\n#      division   1.3850  434\n#      division   1.4734  590\n#      division   1.6010  770\n#      division   4.8366  974\n#      outer_grid  974\n      outer_grid  302\n',
            'La': 'specified\n      division   0.1164  110\n      division   0.8770  194\n      division   0.9952  302\n#      division   1.1042  434\n#      division   1.1747  590\n#      division   1.2496  770\n#      division   4.2775  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Ce': 'specified\n      division   0.1028  110\n      division   0.1495  194\n      division   0.8411  302\n#      division   0.9338  434\n#      division   0.9935  590\n#      division   1.0783  770\n#      division   3.5126  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Pr': 'specified\n      division   0.0809  110\n      division   0.1276  194\n      division   0.7726  302\n#      division   0.8590  434\n#      division   0.9338  590\n#      division   1.0351  770\n#      division   3.3134  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Nd': 'specified\n      division   0.0851  110\n      division   0.1329  194\n      division   0.6933  302\n#      division   0.8063  434\n#      division   0.9338  590\n#      division   1.0141  770\n#      division   3.0511  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Pm': 'specified\n      division   0.0906  110\n      division   0.6102  194\n      division   0.6960  302\n#      division   0.8074  434\n#      division   0.9141  590\n#      division   1.0120  770\n#      division   3.0660  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Sm': 'specified\n      division   0.0906  110\n      division   0.1284  194\n      division   0.6960  302\n#      division   0.8074  434\n#      division   0.9141  590\n#      division   1.0120  770\n#      division   3.0660  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Eu': 'specified\n      division   0.1049  110\n      division   0.1343  194\n      division   0.6986  302\n#      division   0.8085  434\n#      division   0.9136  590\n#      division   1.0099  770\n#      division   3.0016  974\n#      outer_grid   974\n      outer_grid  302\n',
            'Gd': 'specified\n      division   0.0917  110\n      division   0.1291  194\n      division   0.7135  302\n#      division   0.8085  434\n#      division   0.9322  590\n#      division   1.0099  770\n#      division   2.9262  974\n#      outer_grid   974\n      outer_grid  302\n',
            'Tb': 'specified\n      division   0.0876  110\n      division   0.1343  194\n      division   0.7135  302\n#      division   0.8085  434\n#      division   0.9322  590\n#      division   1.0509  770\n#      division   3.0016  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Dy': 'specified\n      division   0.0887  110\n      division   0.1200  194\n      division   0.7773  302\n#      division   0.8774  434\n#      division   0.9501  590\n#      division   1.0482  770\n#      division   3.1772  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Ho': 'specified\n      division   0.0971  110\n      division   0.1105  194\n      division   0.1402  302\n#      division   0.8774  434\n#      division   0.9690  590\n#      division   1.0482  770\n#      division   3.1772  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Er': 'specified\n      division   0.1015  110\n      division   0.1349  194\n      division   0.8600  302\n#      division   0.9314  434\n#      division   1.0079  590\n#      division   1.1114  770\n#      division   3.2637  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Tm': 'specified\n      division   0.1069  110\n      division   0.1797  194\n      division   1.0059  302\n#      division   1.0865  434\n#      division   1.1732  590\n#      division   1.2665  770\n#      division   7.7895  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Yb': 'specified\n      division   0.1305  110\n      division   0.8949  194\n      division   1.1509  302\n#      division   1.2665  434\n#      division   1.3413  590\n#      division   1.4207  770\n#      division   7.7895  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Lu': 'specified\n      division   0.0940  110\n      division   0.8603  194\n      division   0.9866  302\n#      division   1.1076  434\n#      division   1.1732  590\n#      division   1.2190  770\n#      division   7.7895  974\n#      outer_grid  974\n      outer_grid  302\n',
            'Hf': 'specified\n      division   0.4447   50\n      division   1.1303  110\n      division   1.4795  194\n      division   1.7333  302\n#      division   1.9508  434\n#      division   2.2755  590\n#      division   2.4268  770\n#      division   2.5081  974\n#      division   3.0443 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Ta': 'specified\n      division   0.3792   50\n      division   1.0232  110\n      division   1.3396  194\n      division   1.5892  302\n#      division   1.8380  434\n#      division   2.1374  590\n#      division   2.2049  770\n#      division   2.2755  974\n#      division   2.8291 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'W': 'specified\n      division   0.3522   50\n      division   0.9662  110\n      division   1.2839  194\n      division   1.5443  302\n#      division   1.7847  434\n#      division   2.0413  590\n#      division   2.1047  770\n#      division   2.1708  974\n#      division   2.7309 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Re': 'specified\n      division   0.3533   50\n      division   0.9557  110\n      division   1.3010  194\n      division   1.6061  302\n#      division   1.8277  434\n#      division   2.0267  590\n#      division   2.0887  770\n#      division   2.1534  974\n#      division   2.6985 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Os': 'specified\n      division   0.3468   50\n      division   0.9290  110\n      division   1.2652  194\n      division   1.5835  302\n#      division   1.8546  434\n#      division   1.9966  590\n#      division   2.0267  770\n#      division   2.0887  974\n#      division   2.6529 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Ir': 'specified\n      division   0.3664   50\n      division   0.9423  110\n      division   1.2652  194\n      division   1.6525  302\n#      division   1.8819  434\n#      division   2.0267  590\n#      division   2.0887  770\n#      division   2.1534  974\n#      division   2.6985 1202\n#      outer_grid  974\n      outer_grid  302\n',
            'Pt': 'specified\n      division   0.4222   50\n      division   0.9557  110\n      division   1.2477  194\n      division   1.5393  302\n#      division   1.9382  434\n#      division   2.0887  590\n#      division   2.1534  770\n#      division   2.2208  974\n#      division   2.6985 1202\n#      outer_grid    974\n      outer_grid    302\n',
            'Au': 'specified\n      division   0.5066   50\n      division   0.9861  110\n      division   1.2821  194\n      division   1.5344  302\n#      division   2.0427  434\n#      division   2.1690  590\n#      division   2.2710  770\n#      division   2.3066  974\n#      division   2.7597 1202\n#      outer_grid 974\n      outer_grid 302\n',
            'Hg': 'specified\n      division   0.5485   50\n      division   1.0425  110\n      division   1.3361  194\n      division   1.5779  302\n#      division   2.1690  434\n#      division   2.2710  590\n#      division   2.3801  770\n#      division   2.4181  974\n#      division   2.8573 1202\n#      outer_grid  974\n      outer_grid 302\n',
            'Tl': 'specified\n      division   0.1054  110\n      division   0.1577  194\n      division   0.2156  302\n#      division   1.0186  434\n#      division   1.1590  590\n#      division   1.2472  770\n#      division   7.7807  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Pb': 'specified\n      division   0.1064  110\n      division   0.1579  194\n      division   0.2150  302\n#      division   1.0164  434\n#      division   1.1544  590\n#      division   1.1970  770\n#      division   7.7779  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Bi': 'specified\n      division   0.1064  110\n      division   0.1579  194\n      division   0.2150  302\n#      division   1.0164  434\n#      division   1.1133  590\n#      division   1.1970  770\n#      division   7.7779  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Po': 'specified\n      division   0.1022  110\n      division   0.1528  194\n      division   0.2150  302\n#      division   1.0164  434\n#      division   1.1133  590\n#      division   1.1755  770\n#      division   7.7779  974\n#      outer_grid  974\n      outer_grid 302\n',
            'At': 'specified\n      division   0.1106  110\n      division   0.1579  194\n      division   1.0736  302\n#      division   1.1970  434\n#      division   1.2869  590\n#      division   1.4091  770\n#      division   7.7779  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Rn': 'specified\n      division   0.5994  110\n      division   0.8610  194\n      division   1.0898  302\n#      division   1.2801  434\n#      division   1.4253  590\n#      division   7.7751  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Fr': 'specified\n      division   0.5994  110\n      division   0.8769  194\n      division   1.1095  302\n#      division   1.2801  434\n#      division   1.4253  590\n#      division   7.7751  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Ra': 'specified\n      division   0.6236  110\n      division   0.9264  194\n      division   1.1500  302\n#      division   1.3507  434\n#      division   1.5599  590\n#      division   1.6475  770\n#      outer_grid  974\n      outer_grid 302\n',
            'Ac': 'specified\n      division   0.1040  110\n      division   0.1635  194\n      division   1.0865  302\n#      division   1.2735  434\n#      division   1.3667  590\n#      division   1.4413  770\n#      division   7.7724  974\n#      outer_grid 974\n      outer_grid 302\n',
            'Th': 'specified\n      division   0.0887  110\n      division   0.1635  194\n      division   0.9942  302\n#      division   1.0302  434\n#      division   1.1660  590\n#      division   1.2294  770\n#      division   7.7724  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Pa': 'specified\n      division   0.0746  110\n      division   0.1252  194\n      division   0.1687  302\n#      division   0.1905  434\n#      division   0.9942  590\n#      division   1.0865  770\n#      division   7.7724  974\n#      outer_grid  974\n      outer_grid 302\n',
            'U': 'specified\n      division   0.0850  110\n      division   0.1081  194\n      division   0.1389  302\n#      division   0.1794  434\n#      division   0.9255  590\n#      division   1.0302  770\n#      division   7.7724  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Np': 'specified\n      division   0.0861  110\n      division   0.1172  194\n      division   0.1637  302\n#      division   0.1740  434\n#      division   0.9579  590\n#      division   1.0832  770\n#      division   7.7698  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Pu': 'specified\n      division   0.0660  110\n      division   0.0897  194\n      division   0.1215  302\n#      division   0.8157  434\n#      division   0.9246  590\n#      division   1.0644  770\n#      division   7.7698  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Am': 'specified\n      division   0.0757  110\n      division   0.1049  194\n      division   0.1394  302\n#      division   0.8765  434\n#      division   0.9579  590\n#      division   1.1022  770\n#      division   7.7698  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Cm': 'specified\n      division   0.0861  110\n      division   0.1215  194\n      division   0.1394  302\n#      division   0.1793  434\n#      division   0.9579  590\n#      division   1.0832  770\n#      division   7.7698  974\n#      outer_grid 974\n      outer_grid 302\n',
            'Bk': 'specified\n      division   0.0801  110\n      division   0.1179  194\n      division   0.1540  302\n#      division   0.2184  434\n#      division   0.9565  590\n#      division   1.0799  770\n#      division   7.7672  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Cf': 'specified\n      division   0.0703  110\n      division   0.1058  194\n      division   0.1445  302\n#      division   0.2305  434\n#      division   0.9565  590\n#      division   1.0799  770\n#      division   7.7672  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Es': 'specified\n      division   0.0871  110\n      division   0.1221  194\n      division   0.1398  302\n#      division   0.1793  434\n#      division   0.9565  590\n#      division   1.0799  770\n#      division   7.7672  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Fm': 'specified\n      division   0.0735  110\n      division   0.1138  194\n      division   0.1445  302\n#      division   0.2429  434\n#      division   1.0077  590\n#      division   1.0799  770\n#      division   7.7672  974\n#      outer_grid  974\n      outer_grid 302\n',
            'Md': 'specified\n      division   0.0916  110\n      division   0.1314  194\n      division   0.1740  302\n#      division   1.0230  434\n#      division   1.0952  590\n#      division   1.1723  770\n#      division   7.7647  974\n#      outer_grid  974\n      outer_grid 302\n',
            'No': 'specified\n      division   0.1145  110\n      division   0.1640  194\n      division   1.0952  302\n#      division   1.2546  590\n#      division   1.4376  770\n#      division   7.7647  974\n#      outer_grid  974\n      outer_grid 302\n'}
        self.__angular_grids = self.__angular_grids_dic[name_species]
        self.__radial_multiplier_dic = {'H': '1', 'He': '1', 'Li': '1', 'Be': '1', 'B': '1', 'C': '1', 'N': '1',
                                        'O': '1', 'F': '1', 'Ne': '1', 'Na': '1', 'Mg': '1', 'Al': '1', 'Si': '1',
                                        'P': '1', 'S': '1', 'Cl': '1', 'Ar': '1', 'K': '1', 'Ca': '1', 'Sc': '1',
                                        'Ti': '1', 'V': '1', 'Cr': '1', 'Mn': '1', 'Fe': '1', 'Co': '1', 'Ni': '1',
                                        'Cu': '1', 'Zn': '1', 'Ga': '1', 'Ge': '1', 'As': '1', 'Se': '1', 'Br': '1',
                                        'Kr': '1', 'Rb': '1', 'Sr': '1', 'Y': '1', 'Zr': '1', 'Nb': '1', 'Mo': '1',
                                        'Tc': '1', 'Ru': '1', 'Rh': '1', 'Pd': '1', 'Ag': '1', 'Cd': '1', 'In': '1',
                                        'Sn': '1', 'Sb': '1', 'Te': '1', 'I': '1', 'Xe': '1', 'Cs': '1', 'Ba': '1',
                                        'La': '1', 'Ce': '1', 'Pr': '1', 'Nd': '1', 'Pm': '1', 'Sm': '1', 'Eu': '1',
                                        'Gd': '1', 'Tb': '1', 'Dy': '1', 'Ho': '1', 'Er': '1', 'Tm': '1', 'Yb': '1',
                                        'Lu': '1', 'Hf': '1', 'Ta': '1', 'W': '1', 'Re': '1', 'Os': '1', 'Ir': '1',
                                        'Pt': '1', 'Au': '1', 'Hg': '1', 'Tl': '1', 'Pb': '1', 'Bi': '1', 'Po': '1',
                                        'At': '1', 'Rn': '1', 'Fr': '1', 'Ra': '1', 'Ac': '1', 'Th': '1', 'Pa': '1',
                                        'U': '1', 'Np': '1', 'Pu': '1', 'Am': '1', 'Cm': '1', 'Bk': '1', 'Cf': '1',
                                        'Es': '1', 'Fm': '1', 'Md': '1', 'No': '1'}
        self.__radial_multiplier = self.__radial_multiplier_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getAngularGrids(self):
        return self.__angular_grids

    def getRadialBase(self):
        return self.__radial_base

    def getRadialMultiplier(self):
        return self.__radial_multiplier


class MinimalBasis:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__valence_dic = {'H': '#     valence basis states\n    valence      1  s   1.\n',
                              'He': '#     valence basis states\n    valence      1  s   2.\n',
                              'Li': '#     valence basis states\n    valence      2  s   1.\n',
                              'Be': '#     valence basis states\n    valence      2  s   2.\n',
                              'B': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   1.\n',
                              'C': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   2.\n',
                              'N': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   3.\n',
                              'O': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   4.\n',
                              'F': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   5.\n',
                              'Ne': '#     valence basis states\n    valence      2  s   2.\n    valence      2  p   6.\n',
                              'Na': '#     valence basis states\n    valence      3  s   1.\n    valence      2  p   6.\n',
                              'Mg': '#     valence basis states\n    valence      3  s   2.\n    valence      2  p   6.\n',
                              'Al': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   1.\n',
                              'Si': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   2.\n',
                              'P': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   3.\n',
                              'S': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   4.\n',
                              'Cl': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   5.\n',
                              'Ar': '#     valence basis states\n    valence      3  s   2.\n    valence      3  p   6.\n',
                              'K': '#     valence basis states\n    valence      4  s   1.\n    valence      3  p   6.\n',
                              'Ca': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n',
                              'Sc': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   1.\n',
                              'Ti': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   2.\n',
                              'V': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   3.\n',
                              'Cr': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   4.\n',
                              'Mn': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   5.\n',
                              'Fe': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   6.\n',
                              'Co': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   7.\n',
                              'Ni': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d   8.\n',
                              'Cu': '#     valence basis states\n    valence      4  s   1.\n    valence      3  p   6.\n    valence      3  d  10.\n',
                              'Zn': '#     valence basis states\n    valence      4  s   2.\n    valence      3  p   6.\n    valence      3  d  10.\n',
                              'Ga': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   1.\n    valence      3  d  10.\n',
                              'Ge': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   2.\n    valence      3  d  10.\n',
                              'As': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   3.\n    valence      3  d  10.\n',
                              'Se': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   4.\n    valence      3  d  10.\n',
                              'Br': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   5.\n    valence      3  d  10.\n',
                              'Kr': '#     valence basis states\n    valence      4  s   2.\n    valence      4  p   6.\n    valence      3  d  10.\n',
                              'Rb': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      3  d  10.\n',
                              'Sr': '#     valence basis states\n    valence      5  s   2.\n    valence      4  p   6.\n    valence      3  d  10.\n',
                              'Y': '#     valence basis states\n    valence      5  s   2.\n    valence      4  p   6.\n    valence      4  d   1.\n',
                              'Zr': '#     valence basis states\n    valence      5  s   2.\n    valence      4  p   6.\n    valence      4  d   2.\n',
                              'Nb': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   4.\n',
                              'Mo': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   5.\n',
                              'Tc': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   6.\n',
                              'Ru': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   7.\n',
                              'Rh': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   8.\n',
                              'Pd': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d   9.\n',
                              'Ag': '#     valence basis states\n    valence      5  s   1.\n    valence      4  p   6.\n    valence      4  d  10.\n',
                              'Cd': '#     valence basis states\n    valence      5  s   2.\n    valence      4  p   6.\n    valence      4  d  10.\n',
                              'In': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   1.\n    valence      4  d  10.\n',
                              'Sn': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   2.\n    valence      4  d  10.\n',
                              'Sb': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   3.\n    valence      4  d  10.\n',
                              'Te': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   4.\n    valence      4  d  10.\n',
                              'I': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   5.\n    valence      4  d  10.\n',
                              'Xe': '#     valence basis states\n    valence      5  s   2.\n    valence      5  p   6.\n    valence      4  d  10.\n',
                              'Cs': '#     valence basis states\n    valence      6  s   1.\n    valence      5  p   6.\n    valence      4  d  10.\n',
                              'Ba': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      4  d  10.\n',
                              'La': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   1.\n    valence      4  f   0.\n',
                              'Ce': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   1.\n    valence      4  f   1.\n',
                              'Pr': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   3.\n',
                              'Nd': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   4.\n',
                              'Pm': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   5.\n',
                              'Sm': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   6.\n',
                              'Eu': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   7.\n',
                              'Gd': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   1.\n    valence      4  f   7.\n',
                              'Tb': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f   9.\n',
                              'Dy': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f  10.\n',
                              'Ho': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f  11.\n',
                              'Er': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f  12.\n',
                              'Tm': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f  13.\n',
                              'Yb': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   0.\n    valence      4  f  14.\n',
                              'Lu': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   1.\n    valence      4  f  14.\n',
                              'Hf': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   2.\n    valence      4  f  14.\n',
                              'Ta': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   3.\n    valence      4  f  14.\n',
                              'W': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   4.\n    valence      4  f  14.\n',
                              'Re': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   5.\n    valence      4  f  14.\n',
                              'Os': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   6.\n    valence      4  f  14.\n',
                              'Ir': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d   7.\n    valence      4  f  14.\n',
                              'Pt': '#     valence basis states\n    valence      6  s   1.\n    valence      5  p   6.\n    valence      5  d   9.\n    valence      4  f  14.\n',
                              'Au': '#     valence basis states\n    valence      6  s   1.\n    valence      5  p   6.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Hg': '#     valence basis states\n    valence      6  s   2.\n    valence      5  p   6.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Tl': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   1.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Pb': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   2.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Bi': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   3.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Po': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   4.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'At': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   5.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Rn': '#     valence basis states\n    valence      6  s   2.\n    valence      6  p   6.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Fr': '#     valence basis states\n    valence      7  s   1.\n    valence      6  p   6.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Ra': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      5  d  10.\n    valence      4  f  14.\n',
                              'Ac': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   1.\n    valence      5  f   0.\n',
                              'Th': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   2.\n    valence      5  f   0.\n',
                              'Pa': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   1.\n    valence      5  f   2.\n',
                              'U': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   1.\n    valence      5  f   3.\n',
                              'Np': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   1.\n    valence      5  f   4.\n',
                              'Pu': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f   6.\n',
                              'Am': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f   7.\n',
                              'Cm': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   1.\n    valence      5  f   7.\n',
                              'Bk': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f   9.\n',
                              'Cf': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f  10.\n',
                              'Es': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f  11.\n',
                              'Fm': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f  12.\n',
                              'Md': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f  13.\n',
                              'No': '#     valence basis states\n    valence      7  s   2.\n    valence      6  p   6.\n    valence      6  d   0.\n    valence      5  f  14.\n'}
        self.__valence = self.__valence_dic[name_species]
        self.__ion_occ_dic = {'H': '    ion_occ      1  s   0.5\n', 'He': '    ion_occ      1  s   1.\n',
                              'Li': '    ion_occ      1  s   2.\n', 'Be': '    ion_occ      2  s   1.\n',
                              'B': '    ion_occ      2  s   1.\n',
                              'C': '    ion_occ      2  s   1.\n    ion_occ      2  p   1.\n',
                              'N': '    ion_occ      2  s   1.\n    ion_occ      2  p   2.\n',
                              'O': '    ion_occ      2  s   1.\n    ion_occ      2  p   3.\n',
                              'F': '    ion_occ      2  s   1.\n    ion_occ      2  p   4.\n',
                              'Ne': '    ion_occ      2  s   1.\n    ion_occ      2  p   5.\n',
                              'Na': '    ion_occ      2  s   2.\n    ion_occ      2  p   6.\n',
                              'Mg': '    ion_occ      2  s   2.\n    ion_occ      2  p   6.\n',
                              'Al': '    ion_occ      3  s   1.\n    ion_occ      2  p   6.\n',
                              'Si': '    ion_occ      3  s   1.\n    ion_occ      3  p   1.\n',
                              'P': '    ion_occ      3  s   1.\n    ion_occ      3  p   2.\n',
                              'S': '    ion_occ      3  s   1.\n    ion_occ      3  p   3.\n',
                              'Cl': '    ion_occ      3  s   1.\n    ion_occ      3  p   4.\n',
                              'Ar': '    ion_occ      3  s   1.\n    ion_occ      3  p   5.\n',
                              'K': '    ion_occ      3  s   2.\n    ion_occ      3  p   6.\n',
                              'Ca': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n',
                              'Sc': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n',
                              'Ti': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   1.\n',
                              'V': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   2.\n',
                              'Cr': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   3.\n',
                              'Mn': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   4.\n',
                              'Fe': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   5.\n',
                              'Co': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   6.\n',
                              'Ni': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   7.\n',
                              'Cu': '    ion_occ      4  s   0.\n    ion_occ      3  p   6.\n    ion_occ      3  d   9.\n',
                              'Zn': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d   9.\n',
                              'Ga': '    ion_occ      4  s   1.\n    ion_occ      3  p   6.\n    ion_occ      3  d  10.\n',
                              'Ge': '    ion_occ      4  s   1.\n    ion_occ      4  p   1.\n    ion_occ      3  d  10.\n',
                              'As': '    ion_occ      4  s   1.\n    ion_occ      4  p   2.\n    ion_occ      3  d  10.\n',
                              'Se': '    ion_occ      4  s   1.\n    ion_occ      4  p   3.\n    ion_occ      3  d  10.\n',
                              'Br': '    ion_occ      4  s   1.\n    ion_occ      4  p   4.\n    ion_occ      3  d  10.\n',
                              'Kr': '    ion_occ      4  s   1.\n    ion_occ      4  p   5.\n    ion_occ      3  d  10.\n',
                              'Rb': '    ion_occ      5  s   0.\n    ion_occ      4  p   6.\n    ion_occ      3  d  10.\n',
                              'Sr': '    ion_occ      5  s   1.\n    ion_occ      4  p   6.\n    ion_occ      3  d  10.\n',
                              'Y': '    ion_occ      5  s   1.\n    ion_occ      4  p   6.\n    ion_occ      3  d  10.\n',
                              'Zr': '    ion_occ      5  s   1.\n    ion_occ      4  p   6.\n    ion_occ      4  d   1.\n',
                              'Nb': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   3.\n',
                              'Mo': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   4.\n',
                              'Tc': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   5.\n',
                              'Ru': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   6.\n',
                              'Rh': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   7.\n',
                              'Pd': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   8.\n',
                              'Ag': '    ion_occ     5  s   0.\n    ion_occ     4  p   6.\n    ion_occ     4  d   9.\n',
                              'Cd': '    ion_occ     5  s   1.\n    ion_occ     4  p   6.\n    ion_occ     4  d   9.\n',
                              'In': '    ion_occ     5  s   1.\n    ion_occ     5  p   0.\n    ion_occ     4  d  10.\n',
                              'Sn': '    ion_occ     5  s   1.\n    ion_occ     5  p   1.\n    ion_occ     4  d  10.\n',
                              'Sb': '    ion_occ     5  s   1.\n    ion_occ     5  p   2.\n    ion_occ     4  d  10.\n',
                              'Te': '    ion_occ     5  s   1.\n    ion_occ     5  p   3.\n    ion_occ     4  d  10.\n',
                              'I': '    ion_occ      5  s   1.\n    ion_occ      5  p   4.\n    ion_occ      4  d  10.\n',
                              'Xe': '    ion_occ      5  s   1.\n    ion_occ      5  p   5.\n    ion_occ      4  d  10.\n',
                              'Cs': '    ion_occ      6  s   0.\n    ion_occ      5  p   6.\n    ion_occ      4  d  10.\n',
                              'Ba': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      4  d  10.\n',
                              'La': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   0.\n',
                              'Ce': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   0.\n',
                              'Pr': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   2.\n',
                              'Nd': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   3.\n',
                              'Pm': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   4.\n',
                              'Sm': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   5.\n',
                              'Eu': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   6.\n',
                              'Gd': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   6.\n',
                              'Tb': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   8.\n',
                              'Dy': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f   9.\n',
                              'Ho': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f  10.\n',
                              'Er': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f  11.\n',
                              'Tm': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f  12.\n',
                              'Yb': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f  13.\n',
                              'Lu': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   0.\n    ion_occ      4  f  14.\n',
                              'Hf': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   1.\n    ion_occ      4  f  14.\n',
                              'Ta': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   2.\n    ion_occ      4  f  14.\n',
                              'W': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   3.\n    ion_occ      4  f  14.\n',
                              'Re': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   4.\n    ion_occ      4  f  14.\n',
                              'Os': '    ion_occ      6  s   1.\n    ion_occ      5  p   6.\n    ion_occ      5  d   5.\n    ion_occ      4  f  14.\n',
                              'Ir': '    ion_occ     6  s   1.\n    ion_occ     5  p   6.\n    ion_occ     5  d   6.\n    ion_occ     4  f   14.\n',
                              'Pt': '    ion_occ     6  s   0.\n    ion_occ     5  p   6.\n    ion_occ     5  d   8.\n    ion_occ     4  f   14.\n',
                              'Au': '    ion_occ     6  s   0.\n    ion_occ     5  p   6.\n    ion_occ     5  d   9.\n    ion_occ     4  f   14.\n',
                              'Hg': '    ion_occ     6  s   1.\n    ion_occ     5  p   6.\n    ion_occ     5  d   9.\n    ion_occ     4  f   14.\n',
                              'Tl': '    ion_occ     6  s    1.\n    ion_occ     6  p    0.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Pb': '    ion_occ     6  s    1.\n    ion_occ     6  p    1.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Bi': '    ion_occ     6  s    1.\n    ion_occ     6  p    2.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Po': '    ion_occ     6  s    1.\n    ion_occ     6  p    3.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'At': '    ion_occ     6  s    1.\n    ion_occ     6  p    4.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Rn': '    ion_occ     6  s    1.\n    ion_occ     6  p    5.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Fr': '    ion_occ     7  s    0.\n    ion_occ     6  p    6.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Ra': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     5  d   10.\n    ion_occ     4  f   14.\n',
                              'Ac': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    0.\n',
                              'Th': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    0.\n',
                              'Pa': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    1.\n',
                              'U': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    2.\n',
                              'Np': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    3.\n',
                              'Pu': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    5.\n',
                              'Am': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    6.\n',
                              'Cm': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    6.\n',
                              'Bk': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    8.\n',
                              'Cf': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f    9.\n',
                              'Es': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f   10.\n',
                              'Fm': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f   11.\n',
                              'Md': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f   12.\n',
                              'No': '    ion_occ     7  s    1.\n    ion_occ     6  p    6.\n    ion_occ     6  d    0.\n    ion_occ     5  f   13.\n'}
        self.__ion_occ = self.__ion_occ_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getValence(self):
        return self.__valence

    def getIonOcc(self):
        return self.__ion_occ


class FirstTier:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__ft_first_dic = {'H': '     hydro 2 s 2.1\n', 'He': '     hydro 1 s 0.85\n', 'Li': '     hydro 2 p 1.6\n',
                               'Be': '     ionic 2 p auto\n', 'B': '     hydro 2 p 1.4\n', 'C': '     hydro 2 p 1.7\n',
                               'N': '     hydro 2 p 1.8\n', 'O': '     hydro 2 p 1.8\n', 'F': '     hydro 2 p 1.7\n',
                               'Ne': '     hydro 3 d 6\n', 'Na': '     hydro 2 p 1.2\n', 'Mg': '     hydro 2 p 1.5\n',
                               'Al': '     ionic 3 d auto\n', 'Si': '     hydro 3 d 4.2\n',
                               'P': '     ionic 3 d auto\n', 'S': '     ionic 3 d auto\n',
                               'Cl': '     ionic 3 d auto\n', 'Ar': '     ionic 3 d auto\n',
                               'K': '     ionic 3 d auto\n', 'Ca': '     ionic 3 d auto\n',
                               'Sc': '     hydro 4 f 6.8\n', 'Ti': '     hydro 4 f 8\n', 'V': '     hydro 4 f 9\n',
                               'Cr': '     hydro 4 f 9.6\n', 'Mn': '     hydro 4 f 9.6\n', 'Fe': '     hydro 4 f 9.4\n',
                               'Co': '     hydro 3 p 5.8\n', 'Ni': '     hydro 3 p 6\n', 'Cu': '     ionic 4 p auto\n',
                               'Zn': '     hydro 2 p 1.7\n', 'Ga': '     hydro 2 p 1.2\n', 'Ge': '     hydro 2 p 1.4\n',
                               'As': '     hydro 3 d 4\n', 'Se': '     hydro 3 d 4.3\n', 'Br': '     hydro 3 d 4.6\n',
                               'Kr': '     hydro 3 d 4.5\n', 'Rb': '     hydro 3 d 4.5\n',
                               'Sr': '     ionic 4 d auto\n', 'Y': '     hydro 4 f 5.4\n', 'Zr': '     hydro 4 f 7.2\n',
                               'Nb': '     hydro 4 f 7.8\n', 'Mo': '     hydro 4 f 8.4\n', 'Tc': '     hydro 4 f 8.6\n',
                               'Ru': '     hydro 4 f 8.8\n', 'Rh': '     hydro 4 f 8.6\n',
                               'Pd': '     ionic 5 p auto\n', 'Ag': '     ionic 5 p auto\n',
                               'Cd': '     hydro 2 p 1.6\n', 'In': '     hydro 3 p 2.1\n', 'Sn': '     hydro 2 p 1.3\n',
                               'Sb': '     hydro 3 d 3.5\n', 'Te': '     hydro 3 d 3.7\n', 'I': '     hydro 3 d 4\n',
                               'Xe': '     hydro 3 d 3.8\n', 'Cs': '     hydro 3 d 3.9\n',
                               'Ba': '     ionic 5 d auto\n', 'La': '     hydro 4 d 4.6     \n',
                               'Ce': '     hydro 4 d 4.8\n', 'Pr': '     hydro 3 d 4.9\n', 'Nd': '     hydro 3 d 5\n',
                               'Pm': '     hydro 3 d 5.2\n', 'Sm': '     hydro 3 d 5.2\n', 'Eu': '     hydro 3 d 5.4\n',
                               'Gd': '     hydro 3 d 5.8\n', 'Tb': '     hydro 3 d 5.4\n', 'Dy': '     hydro 3 d 2.2\n',
                               'Ho': '     ionic 5 d auto\n', 'Er': '     ionic 5 d auto\n',
                               'Tm': '     ionic 5 d auto\n', 'Yb': '     hydro 2 p 1\n', 'Lu': '     hydro 2 p 1.3\n',
                               'Hf': '     hydro 4 f 6\n', 'Ta': '     hydro 4 f 7\n', 'W': '     hydro 4 f 7.8\n',
                               'Re': '     hydro 4 f 8\n', 'Os': '     hydro 4 f 8.2\n', 'Ir': '     hydro 4 f 8.2\n',
                               'Pt': '     hydro 4 f 7.4\n', 'Au': '     ionic 6 p auto\n',
                               'Hg': '     hydro 2 p 1.7\n', 'Tl': '     hydro 3 p 2.1\n', 'Pb': '     hydro 3 p 2.3\n',
                               'Bi': '     hydro 2 p 1.4\n', 'Po': '     hydro 3 d 3.5\n', 'At': '     hydro 3 d 3.7\n',
                               'Rn': '     hydro 3 d 3.6\n', 'Fr': '     hydro 3 d 3.6\n',
                               'Ra': '     ionic 6 d auto\n', 'Ac': '     ionic 5 f auto\n',
                               'Th': '     ionic 5 f auto\n', 'Pa': '     hydro 3 d 2.5\n', 'U': '     hydro 3 d 5\n',
                               'Np': '     hydro 3 d 5.2\n', 'Pu': '     hydro 3 d 5\n', 'Am': '     hydro 3 d 5.2\n',
                               'Cm': '     hydro 3 d 2.7\n', 'Bk': '     hydro 3 d 5.2\n', 'Cf': '     hydro 3 d 5.2\n',
                               'Es': '     hydro 3 d 5.2\n', 'Fm': '     hydro 3 d 5.2\n', 'Md': '     hydro 3 d 5.2\n',
                               'No': '     hydro 2 p 1\n'}
        self.__ft_first = self.__ft_first_dic[name_species]
        self.__ft_second_dic = {'H': '     hydro 2 p 3.5\n', 'He': '     hydro 2 p 3.5\n', 'Li': '     hydro 2 s 2\n',
                                'Be': '     hydro 3 s 2.9\n', 'B': '     hydro 3 d 4.8\n', 'C': '     hydro 3 d 6\n',
                                'N': '     hydro 3 d 6.8\n', 'O': '     hydro 3 d 7.6\n', 'F': '     hydro 3 d 7.4\n',
                                'Ne': '     hydro 2 p 1.5\n', 'Na': '     hydro 3 s 1.8\n',
                                'Mg': '     ionic 3 d auto\n', 'Al': '     ionic 3 p auto\n',
                                'Si': '     hydro 2 p 1.4\n', 'P': '     ionic 3 p auto\n', 'S': '     hydro 2 p 1.8\n',
                                'Cl': '     hydro 2 p 1.9\n', 'Ar': '     ionic 4 p auto\n',
                                'K': '     ionic 4 p auto\n', 'Ca': '     ionic 4 p auto\n',
                                'Sc': '     ionic 4 p auto\n', 'Ti': '     hydro 3 d 2.7\n', 'V': '     hydro 3 d 3\n',
                                'Cr': '     hydro 3 d 3.1\n', 'Mn': '     hydro 3 d 3.2\n',
                                'Fe': '     hydro 2 p 2.2\n', 'Co': '     hydro 4 f 8.2\n', 'Ni': '     hydro 4 f 9\n',
                                'Cu': '     hydro 4 f 7.4\n', 'Zn': '     hydro 3 s 2.9\n',
                                'Ga': '     hydro 3 d 3.8\n', 'Ge': '     hydro 3 d 4.3\n',
                                'As': '     hydro 2 p 1.5\n', 'Se': '     hydro 2 p 1.6\n',
                                'Br': '     hydro 2 p 1.7\n', 'Kr': '     hydro 3 p 3.1\n',
                                'Rb': '     hydro 3 p 2.5\n', 'Sr': '     ionic 5 p auto\n',
                                'Y': '     hydro 2 p 1.3\n', 'Zr': '     ionic 4 d auto\n',
                                'Nb': '     hydro 3 d 2.6\n', 'Mo': '     hydro 3 d 2.8\n',
                                'Tc': '     ionic 4 d auto\n', 'Ru': '     ionic 4 d auto\n',
                                'Rh': '     ionic 5 p auto\n', 'Pd': '     hydro 4 f 8\n', 'Ag': '     hydro 4 f 7.6\n',
                                'Cd': '     hydro 4 f 7\n', 'In': '#     hydro 4 f 7.6\n', 'Sn': '     hydro 3 d 3.7\n',
                                'Sb': '     ionic 5 p auto\n', 'Te': '#     hydro 4 f 6\n',
                                'I': '#     hydro 4 f 6.4\n', 'Xe': '#     hydro 4 f 6.2\n',
                                'Cs': '     hydro 4 f 6.4\n', 'Ba': '     ionic 4 f auto\n',
                                'La': '     hydro 4 f 6.2     \n', 'Ce': '     hydro 5 g 11.2\n',
                                'Pr': '     hydro 2 p 1.3\n', 'Nd': '     hydro 5 g 11.2\n',
                                'Pm': '     hydro 5 g 11.6\n', 'Sm': '     hydro 5 g 11.6\n',
                                'Eu': '     hydro 4 f 8.2\n', 'Gd': '     hydro 4 f 9\n', 'Tb': '     hydro 4 f 8.2\n',
                                'Dy': '     hydro 4 f 8\n', 'Ho': '     hydro 4 f 7.8\n', 'Er': '     hydro 4 f 7.4\n',
                                'Tm': '     hydro 3 p 3.2\n', 'Yb': '     hydro 3 d 1.6\n',
                                'Lu': '     ionic 6 d auto\n', 'Hf': '     hydro 3 d 6\n', 'Ta': '     hydro 4 d 5.6\n',
                                'W': '     hydro 4 d 5.8\n', 'Re': '     hydro 3 d 7\n', 'Os': '     ionic 6 p auto\n',
                                'Ir': '     ionic 6 p auto\n', 'Pt': '     ionic 6 p auto\n',
                                'Au': '     hydro 4 f 7.4\n', 'Hg': '     hydro 4 f 7\n', 'Tl': '     hydro 4 f 7.6\n',
                                'Pb': '     hydro 4 f 7.6\n', 'Bi': '     ionic 5 d auto\n', 'Po': '     hydro 4 f 6\n',
                                'At': '     hydro 4 f 6.4\n', 'Rn': '     ionic 5 f auto\n',
                                'Fr': '     hydro 4 f 6.4\n', 'Ra': '     ionic 5 f auto\n', 'Ac': '     hydro 4 d 4\n',
                                'Th': '     hydro 4 f 5.2\n', 'Pa': '     hydro 5 g 10.8\n',
                                'U': '     hydro 5 g 11.6\n', 'Np': '     hydro 5 g 12.4\n',
                                'Pu': '     hydro 5 g 12\n', 'Am': '     hydro 5 g 12.4\n',
                                'Cm': '     hydro 5 g 13.2\n', 'Bk': '     hydro 5 g 12.4\n',
                                'Cf': '     hydro 5 g 12.4\n', 'Es': '     hydro 5 g 12.4\n',
                                'Fm': '     hydro 5 g 12\n', 'Md': '     hydro 5 g 12\n', 'No': '     hydro 3 d 4.7\n'}
        if name_species in self.__ft_second_dic:
            self.__ft_second = self.__ft_second_dic[name_species]
        else:
            self.__ft_second = ''
        self.__ft_third_dic = {'Li': '     hydro 3 d 2.6\n', 'Be': '     hydro 3 d 3.5\n', 'B': '     hydro 2 s 4\n',
                               'C': '     hydro 2 s 4.9\n', 'N': '     hydro 3 s 5.8\n', 'O': '     hydro 3 s 6.4\n',
                               'F': '     hydro 3 s 6.8\n', 'Ne': '     hydro 3 s 7.6\n', 'Na': '     hydro 3 d 3.8\n',
                               'Mg': '     hydro 3 s 2.4\n', 'Al': '#     hydro 4 f 4.7\n',
                               'Si': '#     hydro 4 f 6.2\n', 'P': '#     hydro 4 f 6.2\n', 'S': '#     hydro 4 f 7\n',
                               'Cl': '#     hydro 4 f 7.4\n', 'Ar': '#     hydro 4 f 7.4\n',
                               'K': '     hydro 4 s 3.1\n', 'Ca': '     hydro 3 d 2.3\n', 'Sc': '     ionic 4 d auto\n',
                               'Ti': '     ionic 4 p auto\n', 'V': '     ionic 4 p auto\n',
                               'Cr': '     ionic 4 p auto\n', 'Mn': '     hydro 2 p 2\n',
                               'Fe': '#     hydro 5 g 12.4\n', 'Co': '     hydro 3 d 5.4\n',
                               'Ni': '#     hydro 5 g 12.4\n', 'Cu': '     hydro 3 s 2.6\n',
                               'Zn': '     hydro 4 p 5.4\n', 'Ga': '#     hydro 4 f 6.8\n',
                               'Ge': '#     hydro 4 f 7.4\n', 'As': '#     hydro 4 f 6.8\n',
                               'Se': '#     hydro 4 f 7.2\n', 'Br': '#     hydro 4 f 7.6\n',
                               'Kr': '#     hydro 4 f 7.4\n', 'Rb': '#     hydro 4 f 6.6\n',
                               'Sr': '#     hydro 4 f 5.6\n', 'Y': '     ionic 4 d auto\n',
                               'Zr': '     ionic 5 p auto\n', 'Nb': '     ionic 5 p auto\n',
                               'Mo': '     ionic 5 p auto\n', 'Tc': '     ionic 5 p auto\n',
                               'Ru': '     ionic 5 p auto\n', 'Rh': '     ionic 4 d auto\n',
                               'Pd': '#     hydro 5 g 10\n', 'Ag': '     hydro 3 s 2.6\n', 'Cd': '     hydro 3 s 2.8\n',
                               'In': '     hydro 3 d 3.3\n', 'Sn': '#     hydro 4 f 7.4\n',
                               'Sb': '#     hydro 4 f 6.8\n', 'Te': '     hydro 3 p 2.7\n', 'I': '     hydro 2 p 1.6\n',
                               'Xe': '     ionic 6 p auto\n', 'Cs': '     hydro 3 p 2.3\n',
                               'Ba': '     hydro 3 p 2.7\n', 'La': '     hydro 5 g 10      \n',
                               'Ce': '     hydro 4 f 7.6\n', 'Pr': '     hydro 4 f 8\n', 'Nd': '     hydro 4 f 7.6\n',
                               'Pm': '     hydro 4 f 7.8\n', 'Sm': '     hydro 4 f 7.8\n',
                               'Eu': '     hydro 5 g 11.6\n', 'Gd': '     hydro 5 g 13.2\n',
                               'Tb': '     hydro 2 p 1.4\n', 'Dy': '     hydro 2 p 1.3\n', 'Ho': '     hydro 2 p 1.2\n',
                               'Er': '     hydro 2 p 1.2\n', 'Tm': '     hydro 4 f 7\n', 'Yb': '     hydro 4 f 5.6\n',
                               'Lu': '     hydro 4 f 6.6\n', 'Hf': '     ionic 6 p auto\n',
                               'Ta': '     ionic 6 p auto\n', 'W': '     ionic 6 p auto\n',
                               'Re': '     ionic 6 p auto\n', 'Os': '     ionic 5 d auto\n',
                               'Ir': '#     hydro 5 g 10.8\n', 'Pt': '#     hydro 5 g 9.8\n',
                               'Au': '     ionic 6 s auto\n', 'Hg': '     ionic 6 s auto\n',
                               'Tl': '     hydro 3 d 3.4\n', 'Pb': '     hydro 3 d 3.5\n', 'Bi': '     hydro 4 f 7.6\n',
                               'Po': '     hydro 3 p 2.6\n', 'At': '     hydro 2 p 1.6\n', 'Rn': '     hydro 2 p 1.5\n',
                               'Fr': '     ionic 7 p auto\n', 'Ra': '     hydro 3 p 2.4\n',
                               'Ac': '     ionic 7 p auto\n', 'Th': '     hydro 4 d 4.6\n',
                               'Pa': '     hydro 2 p 1.5\n', 'U': '     hydro 2 p 1.9\n', 'Np': '     hydro 2 p 2\n',
                               'Pu': '     hydro 2 p 1.8\n', 'Am': '     hydro 2 p 1.8\n', 'Cm': '     hydro 4 f 8.8\n',
                               'Bk': '     hydro 2 p 1.7\n', 'Cf': '     hydro 2 p 1.6\n', 'Es': '     hydro 4 f 8\n',
                               'Fm': '     hydro 4 f 8\n', 'Md': '     hydro 4 f 7.6\n', 'No': '     hydro 4 f 5.8\n'}
        if name_species in self.__ft_third_dic:
            self.__ft_third = self.__ft_third_dic[name_species]
        else:
            self.__ft_third = ''
        self.__ft_fourth_dic = {'Al': '     ionic 3 s auto\n', 'Si': '     ionic 3 s auto\n',
                                'P': '#     hydro 5 g 8.6\n', 'S': '     ionic 3 s auto\n',
                                'Cl': '     ionic 3 s auto\n', 'Ar': '     hydro 3 s 4.5\n',
                                'K': '#     hydro 4 f 5.6\n', 'Ca': '#     hydro 4 f 4.8\n',
                                'Sc': '#     hydro 5 g 10.4\n', 'Ti': '#     hydro 5 g 11.6\n',
                                'V': '#     hydro 5 g 12.8\n', 'Cr': '#     hydro 5 g 13.6\n',
                                'Mn': '#     hydro 5 g 13.6\n', 'Fe': '     hydro 3 d 3.1\n',
                                'Co': '#     hydro 5 g 12\n', 'Ni': '     hydro 3 d 5.2\n', 'Cu': '     hydro 3 d 5\n',
                                'Zn': '     hydro 4 f 7.8\n', 'Ga': '     ionic 4 s auto\n',
                                'Ge': '     hydro 3 s 3.4\n', 'As': '     ionic 4 s auto\n',
                                'Se': '     ionic 4 s auto\n', 'Br': '     ionic 4 s auto\n',
                                'Kr': '     hydro 3 s 4.2\n', 'Rb': '     hydro 4 s 2.9\n',
                                'Sr': '     ionic 5 s auto\n', 'Y': '#     hydro 5 g 8.4\n',
                                'Zr': '#     hydro 5 g 10.4\n', 'Nb': '#     hydro 5 g 11.2\n',
                                'Mo': '#     hydro 5 g 12\n', 'Tc': '#     hydro 5 g 12.4\n',
                                'Ru': '#     hydro 5 g 12.4\n', 'Rh': '#     hydro 5 g 11.6\n',
                                'Pd': '     hydro 3 s 2.6\n', 'Ag': '#     hydro 5 g 9.8\n',
                                'Cd': '     hydro 3 p 5.2\n', 'In': '     hydro 3 s 2.9\n',
                                'Sn': '     ionic 5 s auto\n', 'Sb': '     ionic 5 s auto\n',
                                'Te': '     ionic 5 s auto\n', 'I': '     ionic 5 s auto\n',
                                'Xe': '     ionic 6 s auto\n', 'Cs': '     hydro 4 s 2.7\n',
                                'Ba': '     hydro 4 s 3.3\n', 'La': '     hydro 2 p 1.5     \n',
                                'Ce': '     hydro 2 p 1.8\n', 'Pr': '     hydro 5 g 11.2\n',
                                'Nd': '     hydro 2 p 1.4\n', 'Pm': '     hydro 2 p 1.4\n',
                                'Sm': '     hydro 2 p 1.4\n', 'Eu': '     hydro 2 p 1.4\n',
                                'Gd': '     hydro 2 p 1.9\n', 'Tb': '     hydro 5 g 12.4\n',
                                'Dy': '     hydro 5 g 12\n', 'Ho': '     hydro 5 g 11.6\n',
                                'Er': '     hydro 5 g 11.2\n', 'Tm': '     hydro 5 g 10.4\n',
                                'Yb': '     hydro 5 g 8.4\n', 'Lu': '#     hydro 5 g 10.4\n',
                                'Hf': '#     hydro 5 g 10.8\n', 'Ta': '#     hydro 5 g 11.6\n',
                                'W': '#     hydro 5 g 12.4\n', 'Re': '#     hydro 5 g 12\n',
                                'Os': '#     hydro 5 g 12\n', 'Ir': '     hydro 3 d 2.9\n',
                                'Pt': '     ionic 6 s auto\n', 'Au': '#     hydro 5 g 10\n',
                                'Hg': '#     hydro 5 g 9.6\n', 'Tl': '     hydro 3 s 3\n',
                                'Pb': '#     hydro 5 g 9.8\n', 'Bi': '     hydro 3 s 3.3\n',
                                'Po': '     ionic 6 s auto\n', 'At': '#     hydro 5 g 9\n', 'Rn': '#     hydro 5 g 8\n',
                                'Fr': '     hydro 4 s 2.8\n', 'Ra': '#     hydro 5 g 6.8\n',
                                'Ac': '     hydro 5 g 9.8\n', 'Th': '     hydro 5 g 10.4\n', 'Pa': '     hydro 4 f 8\n',
                                'U': '#     hydro 6 h 14.8\n', 'Np': '#     hydro 6 h 15.6\n',
                                'Pu': '     hydro 5 f 7.2\n', 'Am': '     hydro 4 f 8.8\n',
                                'Cm': '     hydro 2 p 2.1\n', 'Bk': '     hydro 4 f 8.6\n',
                                'Cf': '     hydro 4 f 8.4\n', 'Es': '     hydro 2 p 1.5\n',
                                'Fm': '     hydro 2 p 1.5\n', 'Md': '     hydro 2 p 1.2\n',
                                'No': '     hydro 4 s 4.4\n'}
        if name_species in self.__ft_fourth_dic:
            self.__ft_fourth = self.__ft_fourth_dic[name_species]
        else:
            self.__ft_fourth = ''
        self.__ft_fifth_dic = {'P': '     ionic 3 s auto\n', 'Cl': '#     hydro 5 g 10.4\n',
                               'Ca': '     ionic 4 s auto\n', 'Sc': '     ionic 4 s auto\n',
                               'Ti': '     ionic 4 s auto\n', 'V': '     ionic 4 s auto\n',
                               'Cr': '     ionic 4 s auto\n', 'Mn': '     hydro 3 s 3.3\n',
                               'Fe': '     ionic 4 s auto\n', 'Co': '     ionic 4 s auto\n',
                               'Ni': '     ionic 4 s auto\n', 'Cu': '#     hydro 5 g 10.4\n',
                               'Zn': '     hydro 3 d 4.5\n', 'Y': '     ionic 5 s auto\n',
                               'Zr': '     ionic 5 s auto\n', 'Nb': '     ionic 5 s auto\n',
                               'Mo': '     ionic 5 s auto\n', 'Tc': '     ionic 5 s auto    \n',
                               'Ru': '     hydro 3 s 2.4\n', 'Rh': '     hydro 3 s 2.5\n', 'Pd': '     hydro 3 d 2.5\n',
                               'Ag': '     hydro 4 d 8.4\n', 'Cd': '#     hydro 5 g 10.0\n',
                               'La': '     hydro 4 s 4.1     \n', 'Ce': '     hydro 3 s 2.7  \n',
                               'Pr': '     ionic 6 s auto \n', 'Nd': '     hydro 3 s 2.6 \n',
                               'Pm': '     ionic 6 s auto \n', 'Sm': '     ionic 6 s auto \n',
                               'Eu': '     hydro 4 s 4.0  \n', 'Gd': '     hydro 3 s 3.0 \n',
                               'Tb': '     hydro 4 s 4.0   \n', 'Dy': '     hydro 4 s 4.0   \n',
                               'Ho': '     hydro 4 s 4.2  \n', 'Er': '     hydro 4 s 4.0   \n',
                               'Tm': '     hydro 4 s 4\n', 'Yb': '     hydro 4 s 4.2 \n',
                               'Lu': '     hydro 4 s 4.4  \n', 'Hf': '     hydro 4 s 4.7  \n',
                               'Ta': '     ionic 6 s auto\n', 'W': '     ionic 6 s auto\n',
                               'Re': '     ionic 6 s auto\n', 'Os': '     ionic 6 s auto\n',
                               'Ir': '     ionic 6 s auto\n', 'Pt': '     hydro 3 d 2.6\n',
                               'Au': '#     hydro 6 h 12.8\n', 'Hg': '     hydro 4 p 4.7\n',
                               'Tl': '#     hydro 5 g 10\n', 'Pb': '     hydro 3 s 3.2\n',
                               'Bi': '#     hydro 5 g 10.4\n', 'Po': '#     hydro 5 g 9\n',
                               'At': '     ionic 6 s auto\n', 'Rn': '     hydro 3 s 3.6\n',
                               'Fr': '#     hydro 5 g 8.2\n', 'Ra': '     hydro 4 s 3.3\n',
                               'Ac': '     hydro 4 f 5.4\n', 'Th': '     hydro 3 p 3.2\n',
                               'Pa': '     hydro 4 s 4.6  \n', 'U': '     hydro 4 f 8.2\n',
                               'Np': '     hydro 4 f 8.6\n', 'Pu': '     hydro 4 s 4.0 \n',
                               'Am': '     hydro 3 s 2.8  \n', 'Cm': '     hydro 1 s 0.7  \n',
                               'Bk': '     ionic 7 s auto  \n', 'Cf': '     hydro 4 s 4.3  \n',
                               'Es': '     hydro 4 s 4.4  \n', 'Fm': '     hydro 4 s 4.6  \n',
                               'Md': '     hydro 4 s 4.2\n', 'No': '     hydro 5 g 9.8\n'}
        if name_species in self.__ft_fifth_dic:
            self.__ft_fifth = self.__ft_fifth_dic[name_species]
        else:
            self.__ft_fifth = ''
        self.__ft_sixth_dic = {'Cd': '     hydro 3 d 3.8\n', 'Au': '     hydro 3 d 2.5\n',
                               'Hg': '     hydro 4 d 7.8   \n', 'Ac': '     hydro 4 s 3.8  \n',
                               'Th': '     ionic 7 s auto  \n', 'U': '     hydro 3 s 2.8  \n',
                               'Np': '     hydro 3 s 2.9  \n'}
        if name_species in self.__ft_sixth_dic:
            self.__ft_sixth = self.__ft_sixth_dic[name_species]
        else:
            self.__ft_sixth = ''
        self.__ft_s_dic = {'H': '     hydro 2 s 2.1\n', 'He': '     hydro 1 s 0.85\n', 'Li': '     hydro 2 s 2\n',
                           'Be': '     hydro 3 s 2.9\n', 'B': '     hydro 2 s 4\n', 'C': '     hydro 2 s 4.9\n',
                           'N': '     hydro 3 s 5.8\n', 'O': '     hydro 3 s 6.4\n', 'F': '     hydro 3 s 6.8\n',
                           'Ne': '     hydro 3 s 7.6\n', 'Na': '     hydro 3 s 1.8\n', 'Mg': '     hydro 3 s 2.4\n',
                           'Al': '     ionic 3 s auto\n', 'Si': '     ionic 3 s auto\n', 'P': '     ionic 3 s auto\n',
                           'S': '     ionic 3 s auto\n', 'Cl': '     ionic 3 s auto\n', 'Ar': '     hydro 3 s 4.5\n',
                           'K': '     hydro 4 s 3.1\n', 'Ca': '     ionic 4 s auto\n', 'Sc': '     ionic 4 s auto\n',
                           'Ti': '     ionic 4 s auto\n', 'V': '     ionic 4 s auto\n', 'Cr': '     ionic 4 s auto\n',
                           'Mn': '     hydro 3 s 3.3\n', 'Fe': '     ionic 4 s auto\n', 'Co': '     ionic 4 s auto\n',
                           'Ni': '     ionic 4 s auto\n', 'Cu': '     hydro 3 s 2.6\n', 'Zn': '     hydro 3 s 2.9\n',
                           'Ga': '     ionic 4 s auto\n', 'Ge': '     hydro 3 s 3.4\n', 'As': '     ionic 4 s auto\n',
                           'Se': '     ionic 4 s auto\n', 'Br': '     ionic 4 s auto\n', 'Kr': '     hydro 3 s 4.2\n',
                           'Rb': '     hydro 4 s 2.9\n', 'Sr': '     ionic 5 s auto\n', 'Y': '     ionic 5 s auto\n',
                           'Zr': '     ionic 5 s auto\n', 'Nb': '     ionic 5 s auto\n', 'Mo': '     ionic 5 s auto\n',
                           'Tc': '     ionic 5 s auto    \n', 'Ru': '     hydro 3 s 2.4\n',
                           'Rh': '     hydro 3 s 2.5\n', 'Pd': '     hydro 3 s 2.6\n', 'Ag': '     hydro 3 s 2.6\n',
                           'Cd': '     hydro 3 s 2.8\n', 'In': '     hydro 3 s 2.9\n', 'Sn': '     ionic 5 s auto\n',
                           'Sb': '     ionic 5 s auto\n', 'Te': '     ionic 5 s auto\n', 'I': '     ionic 5 s auto\n',
                           'Xe': '     ionic 6 s auto\n', 'Cs': '     hydro 4 s 2.7\n', 'Ba': '     hydro 4 s 3.3\n',
                           'La': '     hydro 4 s 4.1     \n', 'Ce': '     hydro 3 s 2.7  \n',
                           'Pr': '     ionic 6 s auto \n', 'Nd': '     hydro 3 s 2.6 \n',
                           'Pm': '     ionic 6 s auto \n', 'Sm': '     ionic 6 s auto \n',
                           'Eu': '     hydro 4 s 4.0  \n', 'Gd': '     hydro 3 s 3.0 \n',
                           'Tb': '     hydro 4 s 4.0   \n', 'Dy': '     hydro 4 s 4.0   \n',
                           'Ho': '     hydro 4 s 4.2  \n', 'Er': '     hydro 4 s 4.0   \n', 'Tm': '     hydro 4 s 4\n',
                           'Yb': '     hydro 4 s 4.2 \n', 'Lu': '     hydro 4 s 4.4  \n',
                           'Hf': '     hydro 4 s 4.7  \n', 'Ta': '     ionic 6 s auto\n', 'W': '     ionic 6 s auto\n',
                           'Re': '     ionic 6 s auto\n', 'Os': '     ionic 6 s auto\n', 'Ir': '     ionic 6 s auto\n',
                           'Pt': '     ionic 6 s auto\n', 'Au': '     ionic 6 s auto\n', 'Hg': '     ionic 6 s auto\n',
                           'Tl': '     hydro 3 s 3\n', 'Pb': '     hydro 3 s 3.2\n', 'Bi': '     hydro 3 s 3.3\n',
                           'Po': '     ionic 6 s auto\n', 'At': '     ionic 6 s auto\n', 'Rn': '     hydro 3 s 3.6\n',
                           'Fr': '     hydro 4 s 2.8\n', 'Ra': '     hydro 4 s 3.3\n', 'Ac': '     hydro 4 s 3.8  \n',
                           'Th': '     ionic 7 s auto  \n', 'Pa': '     hydro 4 s 4.6  \n',
                           'U': '     hydro 3 s 2.8  \n', 'Np': '     hydro 3 s 2.9  \n', 'Pu': '     hydro 4 s 4.0 \n',
                           'Am': '     hydro 3 s 2.8  \n', 'Cm': '     hydro 1 s 0.7  \n',
                           'Bk': '     ionic 7 s auto  \n', 'Cf': '     hydro 4 s 4.3  \n',
                           'Es': '     hydro 4 s 4.4  \n', 'Fm': '     hydro 4 s 4.6  \n', 'Md': '     hydro 4 s 4.2\n',
                           'No': '     hydro 4 s 4.4\n'}
        if name_species in self.__ft_s_dic:
            self.__ft_s = self.__ft_s_dic[name_species]
        else:
            self.__ft_s = ''
        self.__ft_p_dic = {'H': '     hydro 2 p 3.5\n', 'He': '     hydro 2 p 3.5\n', 'Li': '     hydro 2 p 1.6\n',
                           'Be': '     ionic 2 p auto\n', 'B': '     hydro 2 p 1.4\n', 'C': '     hydro 2 p 1.7\n',
                           'N': '     hydro 2 p 1.8\n', 'O': '     hydro 2 p 1.8\n', 'F': '     hydro 2 p 1.7\n',
                           'Ne': '     hydro 2 p 1.5\n', 'Na': '     hydro 2 p 1.2\n', 'Mg': '     hydro 2 p 1.5\n',
                           'Al': '     ionic 3 p auto\n', 'Si': '     hydro 2 p 1.4\n', 'P': '     ionic 3 p auto\n',
                           'S': '     hydro 2 p 1.8\n', 'Cl': '     hydro 2 p 1.9\n', 'Ar': '     ionic 4 p auto\n',
                           'K': '     ionic 4 p auto\n', 'Ca': '     ionic 4 p auto\n', 'Sc': '     ionic 4 p auto\n',
                           'Ti': '     ionic 4 p auto\n', 'V': '     ionic 4 p auto\n', 'Cr': '     ionic 4 p auto\n',
                           'Mn': '     hydro 2 p 2\n', 'Fe': '     hydro 2 p 2.2\n', 'Co': '     hydro 3 p 5.8\n',
                           'Ni': '     hydro 3 p 6\n', 'Cu': '     ionic 4 p auto\n', 'Zn': '     hydro 4 p 5.4\n',
                           'Ga': '     hydro 2 p 1.2\n', 'Ge': '     hydro 2 p 1.4\n', 'As': '     hydro 2 p 1.5\n',
                           'Se': '     hydro 2 p 1.6\n', 'Br': '     hydro 2 p 1.7\n', 'Kr': '     hydro 3 p 3.1\n',
                           'Rb': '     hydro 3 p 2.5\n', 'Sr': '     ionic 5 p auto\n', 'Y': '     hydro 2 p 1.3\n',
                           'Zr': '     ionic 5 p auto\n', 'Nb': '     ionic 5 p auto\n', 'Mo': '     ionic 5 p auto\n',
                           'Tc': '     ionic 5 p auto\n', 'Ru': '     ionic 5 p auto\n', 'Rh': '     ionic 5 p auto\n',
                           'Pd': '     ionic 5 p auto\n', 'Ag': '     ionic 5 p auto\n', 'Cd': '     hydro 3 p 5.2\n',
                           'In': '     hydro 3 p 2.1\n', 'Sn': '     hydro 2 p 1.3\n', 'Sb': '     ionic 5 p auto\n',
                           'Te': '     hydro 3 p 2.7\n', 'I': '     hydro 2 p 1.6\n', 'Xe': '     ionic 6 p auto\n',
                           'Cs': '     hydro 3 p 2.3\n', 'Ba': '     hydro 3 p 2.7\n',
                           'La': '     hydro 2 p 1.5     \n', 'Ce': '     hydro 2 p 1.8\n',
                           'Pr': '     hydro 2 p 1.3\n', 'Nd': '     hydro 2 p 1.4\n', 'Pm': '     hydro 2 p 1.4\n',
                           'Sm': '     hydro 2 p 1.4\n', 'Eu': '     hydro 2 p 1.4\n', 'Gd': '     hydro 2 p 1.9\n',
                           'Tb': '     hydro 2 p 1.4\n', 'Dy': '     hydro 2 p 1.3\n', 'Ho': '     hydro 2 p 1.2\n',
                           'Er': '     hydro 2 p 1.2\n', 'Tm': '     hydro 3 p 3.2\n', 'Yb': '     hydro 2 p 1\n',
                           'Lu': '     hydro 2 p 1.3\n', 'Hf': '     ionic 6 p auto\n', 'Ta': '     ionic 6 p auto\n',
                           'W': '     ionic 6 p auto\n', 'Re': '     ionic 6 p auto\n', 'Os': '     ionic 6 p auto\n',
                           'Ir': '     ionic 6 p auto\n', 'Pt': '     ionic 6 p auto\n', 'Au': '     ionic 6 p auto\n',
                           'Hg': '     hydro 4 p 4.7\n', 'Tl': '     hydro 3 p 2.1\n', 'Pb': '     hydro 3 p 2.3\n',
                           'Bi': '     hydro 2 p 1.4\n', 'Po': '     hydro 3 p 2.6\n', 'At': '     hydro 2 p 1.6\n',
                           'Rn': '     hydro 2 p 1.5\n', 'Fr': '     ionic 7 p auto\n', 'Ra': '     hydro 3 p 2.4\n',
                           'Ac': '     ionic 7 p auto\n', 'Th': '     hydro 3 p 3.2\n', 'Pa': '     hydro 2 p 1.5\n',
                           'U': '     hydro 2 p 1.9\n', 'Np': '     hydro 2 p 2\n', 'Pu': '     hydro 2 p 1.8\n',
                           'Am': '     hydro 2 p 1.8\n', 'Cm': '     hydro 2 p 2.1\n', 'Bk': '     hydro 2 p 1.7\n',
                           'Cf': '     hydro 2 p 1.6\n', 'Es': '     hydro 2 p 1.5\n', 'Fm': '     hydro 2 p 1.5\n',
                           'Md': '     hydro 2 p 1.2\n', 'No': '     hydro 2 p 1\n'}
        if name_species in self.__ft_p_dic:
            self.__ft_p = self.__ft_p_dic[name_species]
        else:
            self.__ft_p = ''
        self.__ft_d_dic = {'Li': '     hydro 3 d 2.6\n', 'Be': '     hydro 3 d 3.5\n', 'B': '     hydro 3 d 4.8\n',
                           'C': '     hydro 3 d 6\n', 'N': '     hydro 3 d 6.8\n', 'O': '     hydro 3 d 7.6\n',
                           'F': '     hydro 3 d 7.4\n', 'Ne': '     hydro 3 d 6\n', 'Na': '     hydro 3 d 3.8\n',
                           'Mg': '     ionic 3 d auto\n', 'Al': '     ionic 3 d auto\n', 'Si': '     hydro 3 d 4.2\n',
                           'P': '     ionic 3 d auto\n', 'S': '     ionic 3 d auto\n', 'Cl': '     ionic 3 d auto\n',
                           'Ar': '     ionic 3 d auto\n', 'K': '     ionic 3 d auto\n', 'Ca': '     hydro 3 d 2.3\n',
                           'Sc': '     ionic 4 d auto\n', 'Ti': '     hydro 3 d 2.7\n', 'V': '     hydro 3 d 3\n',
                           'Cr': '     hydro 3 d 3.1\n', 'Mn': '     hydro 3 d 3.2\n', 'Fe': '     hydro 3 d 3.1\n',
                           'Co': '     hydro 3 d 5.4\n', 'Ni': '     hydro 3 d 5.2\n', 'Cu': '     hydro 3 d 5\n',
                           'Zn': '     hydro 3 d 4.5\n', 'Ga': '     hydro 3 d 3.8\n', 'Ge': '     hydro 3 d 4.3\n',
                           'As': '     hydro 3 d 4\n', 'Se': '     hydro 3 d 4.3\n', 'Br': '     hydro 3 d 4.6\n',
                           'Kr': '     hydro 3 d 4.5\n', 'Rb': '     hydro 3 d 4.5\n', 'Sr': '     ionic 4 d auto\n',
                           'Y': '     ionic 4 d auto\n', 'Zr': '     ionic 4 d auto\n', 'Nb': '     hydro 3 d 2.6\n',
                           'Mo': '     hydro 3 d 2.8\n', 'Tc': '     ionic 4 d auto\n', 'Ru': '     ionic 4 d auto\n',
                           'Rh': '     ionic 4 d auto\n', 'Pd': '     hydro 3 d 2.5\n', 'Ag': '     hydro 4 d 8.4\n',
                           'Cd': '     hydro 3 d 3.8\n', 'In': '     hydro 3 d 3.3\n', 'Sn': '     hydro 3 d 3.7\n',
                           'Sb': '     hydro 3 d 3.5\n', 'Te': '     hydro 3 d 3.7\n', 'I': '     hydro 3 d 4\n',
                           'Xe': '     hydro 3 d 3.8\n', 'Cs': '     hydro 3 d 3.9\n', 'Ba': '     ionic 5 d auto\n',
                           'La': '     hydro 4 d 4.6     \n', 'Ce': '     hydro 4 d 4.8\n',
                           'Pr': '     hydro 3 d 4.9\n', 'Nd': '     hydro 3 d 5\n', 'Pm': '     hydro 3 d 5.2\n',
                           'Sm': '     hydro 3 d 5.2\n', 'Eu': '     hydro 3 d 5.4\n', 'Gd': '     hydro 3 d 5.8\n',
                           'Tb': '     hydro 3 d 5.4\n', 'Dy': '     hydro 3 d 2.2\n', 'Ho': '     ionic 5 d auto\n',
                           'Er': '     ionic 5 d auto\n', 'Tm': '     ionic 5 d auto\n', 'Yb': '     hydro 3 d 1.6\n',
                           'Lu': '     ionic 6 d auto\n', 'Hf': '     hydro 3 d 6\n', 'Ta': '     hydro 4 d 5.6\n',
                           'W': '     hydro 4 d 5.8\n', 'Re': '     hydro 3 d 7\n', 'Os': '     ionic 5 d auto\n',
                           'Ir': '     hydro 3 d 2.9\n', 'Pt': '     hydro 3 d 2.6\n', 'Au': '     hydro 3 d 2.5\n',
                           'Hg': '     hydro 4 d 7.8   \n', 'Tl': '     hydro 3 d 3.4\n', 'Pb': '     hydro 3 d 3.5\n',
                           'Bi': '     ionic 5 d auto\n', 'Po': '     hydro 3 d 3.5\n', 'At': '     hydro 3 d 3.7\n',
                           'Rn': '     hydro 3 d 3.6\n', 'Fr': '     hydro 3 d 3.6\n', 'Ra': '     ionic 6 d auto\n',
                           'Ac': '     hydro 4 d 4\n', 'Th': '     hydro 4 d 4.6\n', 'Pa': '     hydro 3 d 2.5\n',
                           'U': '     hydro 3 d 5\n', 'Np': '     hydro 3 d 5.2\n', 'Pu': '     hydro 3 d 5\n',
                           'Am': '     hydro 3 d 5.2\n', 'Cm': '     hydro 3 d 2.7\n', 'Bk': '     hydro 3 d 5.2\n',
                           'Cf': '     hydro 3 d 5.2\n', 'Es': '     hydro 3 d 5.2\n', 'Fm': '     hydro 3 d 5.2\n',
                           'Md': '     hydro 3 d 5.2\n', 'No': '     hydro 3 d 4.7\n'}
        if name_species in self.__ft_d_dic:
            self.__ft_d = self.__ft_d_dic[name_species]
        else:
            self.__ft_d = ''
        self.__ft_f_dic = {'Al': '     hydro 4 f 4.7\n', 'Si': '     hydro 4 f 6.2\n', 'P': '     hydro 4 f 6.2\n',
                           'S': '     hydro 4 f 7\n', 'Cl': '     hydro 4 f 7.4\n', 'Ar': '     hydro 4 f 7.4\n',
                           'K': '     hydro 4 f 5.6\n', 'Ca': '     hydro 4 f 4.8\n', 'Sc': '     hydro 4 f 6.8\n',
                           'Ti': '     hydro 4 f 8\n', 'V': '     hydro 4 f 9\n', 'Cr': '     hydro 4 f 9.6\n',
                           'Mn': '     hydro 4 f 9.6\n', 'Fe': '     hydro 4 f 9.4\n', 'Co': '     hydro 4 f 8.2\n',
                           'Ni': '     hydro 4 f 9\n', 'Cu': '     hydro 4 f 7.4\n', 'Zn': '     hydro 4 f 7.8\n',
                           'Ga': '     hydro 4 f 6.8\n', 'Ge': '     hydro 4 f 7.4\n', 'As': '     hydro 4 f 6.8\n',
                           'Se': '     hydro 4 f 7.2\n', 'Br': '     hydro 4 f 7.6\n', 'Kr': '     hydro 4 f 7.4\n',
                           'Rb': '     hydro 4 f 6.6\n', 'Sr': '     hydro 4 f 5.6\n', 'Y': '     hydro 4 f 5.4\n',
                           'Zr': '     hydro 4 f 7.2\n', 'Nb': '     hydro 4 f 7.8\n', 'Mo': '     hydro 4 f 8.4\n',
                           'Tc': '     hydro 4 f 8.6\n', 'Ru': '     hydro 4 f 8.8\n', 'Rh': '     hydro 4 f 8.6\n',
                           'Pd': '     hydro 4 f 8\n', 'Ag': '     hydro 4 f 7.6\n', 'Cd': '     hydro 4 f 7\n',
                           'In': '     hydro 4 f 7.6\n', 'Sn': '     hydro 4 f 7.4\n', 'Sb': '     hydro 4 f 6.8\n',
                           'Te': '     hydro 4 f 6\n', 'I': '     hydro 4 f 6.4\n', 'Xe': '     hydro 4 f 6.2\n',
                           'Cs': '     hydro 4 f 6.4\n', 'Ba': '     ionic 4 f auto\n',
                           'La': '     hydro 4 f 6.2     \n', 'Ce': '     hydro 4 f 7.6\n', 'Pr': '     hydro 4 f 8\n',
                           'Nd': '     hydro 4 f 7.6\n', 'Pm': '     hydro 4 f 7.8\n', 'Sm': '     hydro 4 f 7.8\n',
                           'Eu': '     hydro 4 f 8.2\n', 'Gd': '     hydro 4 f 9\n', 'Tb': '     hydro 4 f 8.2\n',
                           'Dy': '     hydro 4 f 8\n', 'Ho': '     hydro 4 f 7.8\n', 'Er': '     hydro 4 f 7.4\n',
                           'Tm': '     hydro 4 f 7\n', 'Yb': '     hydro 4 f 5.6\n', 'Lu': '     hydro 4 f 6.6\n',
                           'Hf': '     hydro 4 f 6\n', 'Ta': '     hydro 4 f 7\n', 'W': '     hydro 4 f 7.8\n',
                           'Re': '     hydro 4 f 8\n', 'Os': '     hydro 4 f 8.2\n', 'Ir': '     hydro 4 f 8.2\n',
                           'Pt': '     hydro 4 f 7.4\n', 'Au': '     hydro 4 f 7.4\n', 'Hg': '     hydro 4 f 7\n',
                           'Tl': '     hydro 4 f 7.6\n', 'Pb': '     hydro 4 f 7.6\n', 'Bi': '     hydro 4 f 7.6\n',
                           'Po': '     hydro 4 f 6\n', 'At': '     hydro 4 f 6.4\n', 'Rn': '     ionic 5 f auto\n',
                           'Fr': '     hydro 4 f 6.4\n', 'Ra': '     ionic 5 f auto\n', 'Ac': '     hydro 4 f 5.4\n',
                           'Th': '     hydro 4 f 5.2\n', 'Pa': '     hydro 4 f 8\n', 'U': '     hydro 4 f 8.2\n',
                           'Np': '     hydro 4 f 8.6\n', 'Pu': '     hydro 5 f 7.2\n', 'Am': '     hydro 4 f 8.8\n',
                           'Cm': '     hydro 4 f 8.8\n', 'Bk': '     hydro 4 f 8.6\n', 'Cf': '     hydro 4 f 8.4\n',
                           'Es': '     hydro 4 f 8\n', 'Fm': '     hydro 4 f 8\n', 'Md': '     hydro 4 f 7.6\n',
                           'No': '     hydro 4 f 5.8\n'}
        if name_species in self.__ft_f_dic:
            self.__ft_f = self.__ft_f_dic[name_species]
        else:
            self.__ft_f = ''
        self.__ft_g_dic = {'P': '     hydro 5 g 8.6\n', 'Cl': '     hydro 5 g 10.4\n', 'Sc': '     hydro 5 g 10.4\n',
                           'Ti': '     hydro 5 g 11.6\n', 'V': '     hydro 5 g 12.8\n', 'Cr': '     hydro 5 g 13.6\n',
                           'Mn': '     hydro 5 g 13.6\n', 'Fe': '     hydro 5 g 12.4\n', 'Co': '     hydro 5 g 12\n',
                           'Ni': '     hydro 5 g 12.4\n', 'Cu': '     hydro 5 g 10.4\n', 'Y': '     hydro 5 g 8.4\n',
                           'Zr': '     hydro 5 g 10.4\n', 'Nb': '     hydro 5 g 11.2\n', 'Mo': '     hydro 5 g 12\n',
                           'Tc': '     hydro 5 g 12.4\n', 'Ru': '     hydro 5 g 12.4\n', 'Rh': '     hydro 5 g 11.6\n',
                           'Pd': '     hydro 5 g 10\n', 'Ag': '     hydro 5 g 9.8\n', 'Cd': '     hydro 5 g 10.0\n',
                           'La': '     hydro 5 g 10      \n', 'Ce': '     hydro 5 g 11.2\n',
                           'Pr': '     hydro 5 g 11.2\n', 'Nd': '     hydro 5 g 11.2\n', 'Pm': '     hydro 5 g 11.6\n',
                           'Sm': '     hydro 5 g 11.6\n', 'Eu': '     hydro 5 g 11.6\n', 'Gd': '     hydro 5 g 13.2\n',
                           'Tb': '     hydro 5 g 12.4\n', 'Dy': '     hydro 5 g 12\n', 'Ho': '     hydro 5 g 11.6\n',
                           'Er': '     hydro 5 g 11.2\n', 'Tm': '     hydro 5 g 10.4\n', 'Yb': '     hydro 5 g 8.4\n',
                           'Lu': '     hydro 5 g 10.4\n', 'Hf': '     hydro 5 g 10.8\n', 'Ta': '     hydro 5 g 11.6\n',
                           'W': '     hydro 5 g 12.4\n', 'Re': '     hydro 5 g 12\n', 'Os': '     hydro 5 g 12\n',
                           'Ir': '     hydro 5 g 10.8\n', 'Pt': '     hydro 5 g 9.8\n', 'Au': '     hydro 5 g 10\n',
                           'Hg': '     hydro 5 g 9.6\n', 'Tl': '     hydro 5 g 10\n', 'Pb': '     hydro 5 g 9.8\n',
                           'Bi': '     hydro 5 g 10.4\n', 'Po': '     hydro 5 g 9\n', 'At': '     hydro 5 g 9\n',
                           'Rn': '     hydro 5 g 8\n', 'Fr': '     hydro 5 g 8.2\n', 'Ra': '     hydro 5 g 6.8\n',
                           'Ac': '     hydro 5 g 9.8\n', 'Th': '     hydro 5 g 10.4\n', 'Pa': '     hydro 5 g 10.8\n',
                           'U': '     hydro 5 g 11.6\n', 'Np': '     hydro 5 g 12.4\n', 'Pu': '     hydro 5 g 12\n',
                           'Am': '     hydro 5 g 12.4\n', 'Cm': '     hydro 5 g 13.2\n', 'Bk': '     hydro 5 g 12.4\n',
                           'Cf': '     hydro 5 g 12.4\n', 'Es': '     hydro 5 g 12.4\n', 'Fm': '     hydro 5 g 12\n',
                           'Md': '     hydro 5 g 12\n', 'No': '     hydro 5 g 9.8\n'}
        if name_species in self.__ft_g_dic:
            self.__ft_g = self.__ft_g_dic[name_species]
        else:
            self.__ft_g = ''
        self.__ft_h_dic = {'Au': '     hydro 6 h 12.8\n', 'U': '     hydro 6 h 14.8\n', 'Np': '     hydro 6 h 15.6\n'}
        if name_species in self.__ft_h_dic:
            self.__ft_h = self.__ft_h_dic[name_species]
        else:
            self.__ft_h = ''

    def getNameSpecies(self):
        return self.__name_species

    def getFirst(self):
        return self.__ft_first

    def getSecond(self):
        return self.__ft_second

    def getThird(self):
        return self.__ft_third

    def getFourth(self):
        return self.__ft_fourth

    def getFifth(self):
        return self.__ft_fifth

    def getSixth(self):
        return self.__ft_sixth

    def getS(self):
        return self.__ft_s

    def getP(self):
        return self.__ft_p

    def getD(self):
        return self.__ft_d

    def getF(self):
        return self.__ft_f

    def getG(self):
        return self.__ft_g

    def getH(self):
        return self.__ft_h

    def getFirstTotalTier(self):
        return '\n#  "First tier" \n     {}     {} {} {} {} {}'.format(self.getFirst(), self.getSecond(),
                                                                       self.getThird(), self.getFourth(),
                                                                       self.getFifth(), self.getSixth())


class SecondTier:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__first_dic = {}
        self.__first = self.__first_dic[name_species]
        self.__second_dic = {}
        self.__second = self.__second_dic[name_species]
        self.__third_dic = {}
        self.__third = self.__second_dic[name_species]
        self.__fourth_dic = {}
        self.__fourth = self.__second_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getFirst(self):
        return self.__first

    def getSecond(self):
        return self.__second

    def getThird(self):
        return self.__third

    def getFourth(self):
        return self.__fourth

    def getSecondTotalTier(self):
        return '\n#  "Second tier" - improvements: -12.89 meV to -1.83 meV'


class ThirdTier:
    def __init__(self, name_species):  # Konstruktor
        self.__name_species = name_species
        self.__first_dic = {}
        self.__first = self.__first_dic[name_species]
        self.__second_dic = {}
        self.__second = self.__second_dic[name_species]
        self.__third_dic = {}
        self.__third = self.__second_dic[name_species]
        self.__fourth_dic = {}
        self.__fourth = self.__second_dic[name_species]

    def getNameSpecies(self):
        return self.__name_species

    def getFirst(self):
        return self.__first

    def getSecond(self):
        return self.__second

    def getThird(self):
        return self.__third

    def getFourth(self):
        return self.__fourth

    def getThirdTotalTier(self):
        return '\n#  "Third tier" - improvements: -0.25 meV to -0.12 meV'


class Species:
    def __init__(self, name_species):
        self.__species = Element(name_species)
        self.__name = name_species
        self.__lhartree = Lhartree(name_species)
        self.__cutoff_pot = CutoffPotential(name_species)
        self.__integration_grid = IntegrationGrid(name_species)
        self.__minimal_basis = MinimalBasis(name_species)
        self.__first_tier = FirstTier(name_species)
        #	self.__second_tier = SecondTier(name_species)
        #		self.__third_tier = ThirdTier(name_species)
        self.__include_fcts = 's'

    def getName(self):
        return self.__name

    def getSpecies(self):
        return self.__species

    def getLhartree(self):
        return self.__lhartree

    def getCutoffPot(self):
        return self.__cutoff_pot

    def getIntegrationGrid(self):
        return self.__integration_grid

    def getMinimalBasis(self):
        return self.__minimal_basis

    def getFirstTier(self):
        return self.__first_tier

    def getSecondTier(self):
        return self.__second_tier

    def getThirdTier(self):
        return self.__third_tier

    def getIncludeFcts(self):
        return self.__include_fcts

    def write_file(self):
        first_tier_string = '\n#  "First tier" \n'
        if self.getIncludeFcts().find('s') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getS())
        if self.getIncludeFcts().find('p') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getP())
        if self.getIncludeFcts().find('d') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getD())
        if self.getIncludeFcts().find('f') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getF())
        if self.getIncludeFcts().find('g') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getG())
        if self.getIncludeFcts().find('h') > -1:
            first_tier_string += '{}'.format(self.getFirstTier().getH())
        if self.getIncludeFcts().find('6') > -1:
            first_tier_string = '{}'.format(self.getFirstTier().getFirstTotalTier())
        if self.getIncludeFcts().find('5') > -1:
            first_tier_string = '     {}     {} {} {} {}'.format(self.getFirstTier().getFirst(),
                                                                 self.getFirstTier().getSecond(),
                                                                 self.getFirstTier().getThird(),
                                                                 self.getFirstTier().getFourth(),
                                                                 self.getFirstTier().getFifth())
        if self.getIncludeFcts().find('4') > -1:
            first_tier_string = '     {}     {} {} {}'.format(self.getFirstTier().getFirst(),
                                                              self.getFirstTier().getSecond(),
                                                              self.getFirstTier().getThird(),
                                                              self.getFirstTier().getFourth())
        if self.getIncludeFcts().find('3') > -1:
            first_tier_string = '     {}     {} {} '.format(self.getFirstTier().getFirst(),
                                                            self.getFirstTier().getSecond(),
                                                            self.getFirstTier().getThird())
        if self.getIncludeFcts().find('2') > -1:
            first_tier_string = '     {}     {}  '.format(self.getFirstTier().getFirst(),
                                                          self.getFirstTier().getSecond())
        if self.getIncludeFcts().find('1') > -1:
            first_tier_string = '     {}      '.format(self.getFirstTier().getFirst())
        if self.getIncludeFcts().find('0') > -1:
            first_tier_string = ' '
        if len(str(self.__species.getAtomicNumber())) == 1:
            f = open('0{}_{}_default'.format(self.__species.getAtomicNumber(), self.getName()), 'w')
        else:
            f = open('{}_{}_default'.format(self.__species.getAtomicNumber(), self.getName()), 'w')
        f.write('################################################################################')
        f.write('\n#')
        f.write('\n#  FHI-aims code project')
        f.write('\n#  VB, Fritz-Haber Institut, 2009')
        f.write('\n#')
        f.write('\n#  Suggested "light" defaults for H atom (to be pasted into control.in file)')
        f.write('\n#  Be sure to double-check any results obtained with these settings for post-processing,')
        f.write('\n#  e.g., with the "tight" defaults and larger basis sets.')
        f.write('\n#')
        f.write('################################################################################')
        f.write('\n  species        {}'.format(self.getSpecies().getNameSpecies()))
        f.write('\n#     global species definitions')
        # TS surf parameters
        if self.getSpecies().getNameSpecies() == 'Ag':
            f.write('\n    hirshfeld_param     {} {} {}'.format(122, 15.4, 2.57))
        f.write('\n    nucleus             {}'.format(self.getSpecies().getNucleus()))
        f.write('\n    mass                {}'.format(self.getSpecies().getMass()))
        f.write('\n#')
        f.write('\n    l_hartree           {}'.format(self.getLhartree().getLhartree()))
        f.write('\n#')
        f.write('\n    cut_pot             {}'.format(self.getCutoffPot().getCutPot()))
        f.write('\n    basis_dep_cutoff    {}'.format(self.getCutoffPot().getBasisDepCutoff()))
        f.write('\n#')
        f.write('\n    radial_base         {}'.format(self.getIntegrationGrid().getRadialBase()))
        f.write('\n    radial_multiplier   {}'.format(self.getIntegrationGrid().getRadialMultiplier()))
        f.write('\n    angular_grids       {}'.format(self.getIntegrationGrid().getAngularGrids()))
        f.write('\n################################################################################')
        f.write('\n#')
        f.write('\n#  Definition of "minimal" basis')
        f.write('\n#')
        f.write('\n#')
        f.write('\n################################################################################')
        f.write('\n{}'.format(self.getMinimalBasis().getValence()))
        f.write('\n#     ion occupancy')
        f.write('\n    {}'.format(self.getMinimalBasis().getIonOcc()))
        f.write('\n################################################################################')
        f.write('\n#')
        f.write('\n#  Suggested additional basis functions. For production calculations, ')
        f.write('\n#  uncomment them one after another (the most important basis functions are')
        f.write('\n#  listed first).')
        f.write('\n#')
        f.write('#  Basis constructed for dimers: 0.5 A, 0.7 A, 1.0 A, 1.5 A, 2.5 A')
        f.write('\n#')
        f.write('\n################################################################################')
        f.write('\n{}'.format(first_tier_string))
        # f.write('\n{}'.format(self.getSecondTier().getSecondTotalTier()))
        # f.write('\n{}'.format(self.getThirdTier().getThirdTotalTier()))
        f.close()
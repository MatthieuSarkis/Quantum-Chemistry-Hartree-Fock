# -*- coding: utf-8 -*-
#
# Written by Matthieu Sarkis, https://github.com/MatthieuSarkis
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

import numpy as np
from src.info import *
from src.molecule import Molecule
from src.primitive_gaussian import PrimitiveGaussian
from src.utilities import frobenius_norm

class HartreeFock():
    
    def __init__(self,
                 molecule: Molecule,
                 number_electrons: int,
                 ) -> None:
        
        self.molecule = molecule
        self.number_electrons = number_electrons
        
        self.total_number_orbitals = 0
        for atom in molecule.atoms:
            self.total_number_orbitals += atom.number_atomic_orbitals
            
        self.S = None
        self.T = None
        self.V = None
        self.coulomb_exchange = None
        self.H_core = None
        self.C = None
        self.density_matrix = None

        self._initialize_matrices()
        
        self.hamiltonian_computed = False
        
    def _initialize_matrices(self) -> None:
        
        self.S = np.zeros((self.total_number_orbitals, self.total_number_orbitals))
        self.T = np.zeros((self.total_number_orbitals, self.total_number_orbitals))
        self.V = np.zeros((self.total_number_orbitals, self.total_number_orbitals))
        self.coulomb_exchange = np.zeros((self.total_number_orbitals, self.total_number_orbitals, self.total_number_orbitals, self.total_number_orbitals))
        self.H_core = np.zeros((self.total_number_orbitals, self.total_number_orbitals))
        self.C = np.zeros((self.total_number_orbitals, int(self.number_electrons / 2)))
        self.density_matrix = np.zeros((self.total_number_orbitals, self.total_number_orbitals))

    def compute_Hamiltonian_contributions(self) -> None:
        
        self._initialize_matrices()
        
        for idx_a, atom_a in enumerate(self.molecule.atoms):
            for idxo_a, orbital_a in enumerate(atom_a.atomic_orbitals):
                for gaussian_a in orbital_a.basis_primitive_gaussians:
                    for idx_b, atom_b in enumerate(self.molecule.atoms):
                        for idxo_b, orbital_b in enumerate(atom_b.atomic_orbitals):
                            for gaussian_b in orbital_b.basis_primitive_gaussians:
                                a = (idx_a + 1) * (idxo_a + 1) - 1
                                b = (idx_b + 1) * (idxo_b + 1) - 1
        
                                self.S[a, b] += gaussian_a.weight * gaussian_b.weight * PrimitiveGaussian.overlap(gaussian_a, gaussian_b)                                
                                self.T[a, b] += gaussian_a.weight * gaussian_b.weight * PrimitiveGaussian.kinetic(gaussian_a, gaussian_b)
                                
                                for atom in self.molecule.atoms:
                                    self.V[a, b] += gaussian_a.weight * gaussian_b.weight * PrimitiveGaussian.potential(gaussian_a, gaussian_b, atom)
                                    
                                for idx_c, atom_c in enumerate(self.molecule.atoms):
                                    for idxo_c, orbital_c in enumerate(atom_c.atomic_orbitals):
                                        for gaussian_c in orbital_c.basis_primitive_gaussians:
                                            for idx_d, atom_d in enumerate(self.molecule.atoms):
                                                for idxo_d, orbital_d in enumerate(atom_d.atomic_orbitals):
                                                    for gaussian_d in orbital_d.basis_primitive_gaussians:
                                                        c = (idx_c + 1) * (idxo_c + 1) - 1
                                                        d = (idx_d + 1) * (idxo_d + 1) - 1
                                                        
                                                        self.coulomb_exchange[a, b, c, d] += gaussian_a.weight * gaussian_b.weight * gaussian_c.weight * gaussian_d.weight *  PrimitiveGaussian.multi(gaussian_a, gaussian_b, gaussian_c, gaussian_d)
                                                        
        self.H_core = self.T + self.V
        self.hamiltonian_computed = True
        
    def _diagonalize_S(self) -> np.array:
        
        _, U = np.linalg.eig(self.S)
        s = U.T.dot(self.S.dot(U))
        s_inverse_root = np.diag(1 / np.sqrt(np.diagonal(s)))
        X = U.dot(s_inverse_root.dot(U.T))
        X = np.dot(U, np.dot(s_inverse_root, U.T))
        
        return X
    
    def self_consistent_field(self,
                              epsilon: float = 10e-5,
                              ) -> None:
            
        if not self.hamiltonian_computed:
            self.compute_Hamiltonian_contributions()
                
        diff = 1.0
        M = self.total_number_orbitals
        self.density_matrix = np.zeros((M, M))
        density_matrix_old = self.density_matrix.copy()
        X = self._diagonalize_S()
        
        while diff > epsilon:
            
            G = np.zeros((M, M))
            
            for i in range(M):
                for j in range(M):
                    for k in range(M):
                        for l in range(M):
                            G[i, j] += self.density_matrix[k, l] * (self.coulomb_exchange[i, j, l, k] - self.coulomb_exchange[i, k, l, j] / 2)
                            
            F = self.H_core + G
            F_Sbasis = np.dot(X.T, np.dot(F, X))
            eigenvalues_F_Sbasis, C_Sbasis = np.linalg.eig(F_Sbasis)
            
            ordered_idx = eigenvalues_F_Sbasis.argsort()
            eigenvalues_F_Sbasis = eigenvalues_F_Sbasis[ordered_idx]
            C_Sbasis = C_Sbasis[:, ordered_idx]
            
            self.C = np.dot(X, C_Sbasis)
            
            for i in range(M):
                for j in range(M):
                    for a in range(np.int(self.number_electrons / 2)):
                        self.density_matrix[i, j] = 2 * self.C[i, a] * self.C[j, a]
            
            diff = self._compare_matrices(density_matrix_old, self.density_matrix)
            density_matrix_old = self.density_matrix.copy()
    
    @staticmethod
    def _compare_matrices(P1: np.array,
                          P2: np.array) -> float:
        
        N = P1.shape[0]
        return frobenius_norm(P1 - P2) / N
            
    def print_matrices(self) -> None:
        
        print("\n")
        print("The system studied is composed of {}, and {} electrons. The basis of primitive Gaussians is {}.\n".format(self.molecule.atoms_list, self.number_electrons, self.molecule.basis))
        print("The overlap matrix S is:")
        print(self.S, "\n")
        print("The kinetic matrix T is:")
        print(self.T, "\n")
        print("The potential matrix V is:")
        print(self.V, "\n")
        print("The orbital coefficients matrix C is:")
        print(self.C, "\n")
        print("The density matrix P is:")
        print(self.density_matrix, "\n")
        
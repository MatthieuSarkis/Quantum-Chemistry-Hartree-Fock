import numpy as np

from src.atom import PrimitiveGaussian, Molecule

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
        self.C = np.zeros((self.total_number_orbitals, np.ceil(self.number_electrons / 2)))
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
        s = np.dot(U.T, np.dot(self.S, U))
        s_minushalf = np.diag(np.diagonal(s)**-0.5)
        X = np.dot(U, np.dot(s_minushalf, U.T))
        return X
    
    def self_consistent_field(self,
                              epsilon: float = 10e-4,
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
                            
            Fock = self.H_core + G

            Fock_Sbasis = np.dot(X.T, np.dot(Fock, X))
            eigenvalues_Fock_Sbasis, C_Sbasis = np.linalg.eig(Fock_Sbasis)
            
            idx = eigenvalues_Fock_Sbasis.argsort()
            eigenvalues_Fock_Sbasis = eigenvalues_Fock_Sbasis[idx]
            C_Sbasis = C_Sbasis[:, idx]
            
            self.C = np.dot(X, C_Sbasis)
            
            for i in range(M):
                for j in range(M):
                    for a in range(np.ceil(self.number_electrons / 2)):
                        self.density_matrix[i, j] = 2 * self.C[i, a] * self.C[j, a]
            
            diff = self._compare_matrices(density_matrix_old, self.density_matrix)
            density_matrix_old = self.density_matrix.copy()
    
    @staticmethod
    def _compare_matrices(P1: np.array,
                          P2: np.array) -> float:
        
        N = P1.shape[0]
        diff = 0
        for i in range(N):
            for j in range(N):
                diff += (1 / N)**2 * (P1[i, j] - P2[i, j])**2
        
        return np.sqrt(diff)
            
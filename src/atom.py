import numpy as np
from typing import List

class PrimitiveGaussian():
    """
    N * \exp(- \alpha |r - R_A|^2)
    """
    
    def __init__(self,
                 center: List[float],
                 alpha: float,
                 weight: float,
                 ) -> None:
        
        self.center = center
        self.alpha = alpha
        self.weight = weight
        self.normalization = (2.0 * self.alpha / np.pi)**0.75
        
    def __mul__(self, g):
        p = self.alpha + g.alpha
        
        return PrimitiveGaussian(center=)  

class AtomicOrbital():
    
    def __init__(self,
                 basis_weights: List[float]) -> None:
        self.basis_weights = basis_weights
        self.basis_functions = []
        
class Atom():
    
    def __init__(self,
                 atom_type: str,
                 coordinates: List[float],
                 zeta_scaling: List[float],
                 number_atomic_orbitals: int) -> None:
        
        self.atom_type = atom_type
        self.coordinates = coordinates
        self.zeta_scaling = zeta_scaling
        self.number_atomic_orbitals = number_atomic_orbitals
        self.atomic_orbitals = []
        
        self._initialize_atomic_orbitals()

    def _initialize_atomic_orbitals(self) -> None:
        pass




if __name__ == "__main__":
    atom = Atom('H', [0,0,0], [1,2], 2)
    orbital = AtomicOrbital([[1,2],[3,4]])
    g1 = PrimitiveGaussian([1,2,3], 1, 1)
    
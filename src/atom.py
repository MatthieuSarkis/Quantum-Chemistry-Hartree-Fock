import numpy as np
from typing import List

class PrimitiveGaussian():
    """
    weight * normalization * \exp(- \alpha |r - center|^2)
    """
    
    def __init__(self,
                 center: List[float],
                 alpha: float,
                 weight: float,
                 normalization: None,
                 ) -> None:
        
        self.center = np.array(center) 
        self.alpha = alpha # related to the standard deviation
        self.weight = weight # weight in the linear combination with other primitive gaussians to form an atomic orbital.
        
        if normalization is None:
            self.normalization = (2.0 * self.alpha / np.pi)**0.75 
        else:
            self.normalization = normalization
        
    def __mul__(self, 
                g: 'PrimitiveGaussian') -> 'PrimitiveGaussian':
        """Overload of the multiplication operator.
        
        The total normalizing factor comes from the normalizing factor of 
        both primitive gaussians, and from completing the square.
        """
        
        alpha = self.alpha + g.alpha
        diff = np.linalg.norm(self.center - g.center)**2
        N = self.normalization * g.normalization
        normalization = N * np.exp(-self.alpha * g.alpha / alpha * diff)
        weight = self.weight * g.weight
        center = (self.alpha * self.center + g.alpha * g.center) / alpha
        
        return PrimitiveGaussian(center=center,
                                 alpha=alpha,
                                 weight=weight,
                                 normalization=normalization)
        
    @staticmethod
    def overlap(g1: 'PrimitiveGaussian',
                g2: 'PrimitiveGaussian') -> float:
        
        g = g1 * g2
        prefactor = (np.pi / g.alpha)**1.5
        
        return prefactor * g.normalization
    
    @staticmethod
    def kinetic(g1: 'PrimitiveGaussian',
                g2: 'PrimitiveGaussian') -> float:
        
        g = g1 * g2
        prefactor = (np.pi / g.alpha)**1.5
        reduced_exponent = g1.alpha * g2.alpha / g.alpha

        return g.normalization * prefactor * reduced_exponent * (3 - 2 * reduced_exponent * np.linalg.norm(g1.center - g2.center)**2)
    
    
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
    
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

import math
import numpy as np
from typing import List

from src.utilities import boys, euclidean_distance2 

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
        diff = euclidean_distance2(self.center, g.center)
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

        return g.normalization * prefactor * reduced_exponent * (3 - 2 * reduced_exponent * euclidean_distance2(g1.center, g2.center))
    
    @staticmethod
    def potential(g1: 'PrimitiveGaussian',
                  g2: 'PrimitiveGaussian',
                  atom: 'Atom') -> float:
        
        assert isinstance(atom, Atom)
        g = g1 * g2
        R = atom.coordinates
        Z = atom.charge
    
        return (-2 * math.pi * Z / g.alpha) * g.normalization * boys(g.alpha * euclidean_distance2(g.center, R))
    
    @staticmethod
    def multi(g1: 'PrimitiveGaussian',
              g2: 'PrimitiveGaussian',
              g3: 'PrimitiveGaussian',
              g4: 'PrimitiveGaussian',
              ) -> float:
        
        ga = g1 * g2
        gb = g3 * g4
        prefactor = 2 * math.pi**2.5 / (ga.alpha * gb.alpha * math.sqrt(ga.alpha + gb.alpha))
    
        return prefactor * ga.normalization * gb.normalization * boys(ga.alpha * gb.alpha / (ga.alpha + gb.alpha) * euclidean_distance2(ga.center, gb.center))
class AtomicOrbital():
    
    def __init__(self,
                 basis_weights: List[float],
                 basis_exponents: List[float],
                 center: List[float],
                 zeta_weight: float,
                 ) -> None:
        
        assert len(basis_weights) == len(basis_exponents)
        self.basis_weights = basis_weights
        self.basis_exponents = basis_exponents
        self.center = center
        self.zeta_weight = zeta_weight 
        
        self.basis_primitive_gaussians = [PrimitiveGaussian(center=self.center,
                                                            alpha=self.basis_exponents[i],
                                                            weight=self.basis_weights[i],
                                                            normalization=None,
                                                            ) for i in range(len(self.basis_weights))]    
class Atom():
    
    zeta = {'H' : [1.24],
            'He': [2.0925],
            'Li': [2.69, 0.80],
            'Be': [3.68, 1.15],
            'B' : [4.68, 1.50],
            'C' : [5.67, 1.72]}

    max_quantum_number = {'H' : 1,
                          'He': 1,
                          'Li': 2,
                          'Be': 2,
                          'C' : 2}
    
    charge = {'H' : 1,
              'He': 2,
              'Li': 3,
              'Be': 4,
              'B' : 5,
              'C' : 6,
              'N' : 7,
              'O' : 8,
              'F' : 9,
              'Ne': 10}
    
    ### STOnG ###
    basis_weights = {'H' :[[0.444635, 0.535328, 0.154329], [0.700115, 0.399513, -0.0999672]]}
    basis_exponents = {'H' : [[0.109818, 0.405771, 2.22766], [0.0751386, 0.231031, 0.994203]]}
    #############
    
    def __init__(self,
                 atom_type: str,
                 coordinates: List[float] = [0.0, 0.0, 0.0],
                 ) -> None:
        
        self.atom_type = atom_type
        self.charge = Atom.charge[self.atom_type]
        self.zeta_scaling = Atom.zeta[self.atom_type]
        self.coordinates = coordinates
        self.number_atomic_orbitals = len(self.zeta_scaling)
        
        self.basis_weights = Atom.basis_weights[atom_type]
        self.basis_exponents = Atom.basis_exponents[atom_type]
        
        self.atomic_orbitals = [AtomicOrbital(basis_weights=self.basis_weights[i],
                                              basis_exponents=self.basis_exponents[i],
                                              center=self.coordinates,
                                              zeta_weight=self.zeta_scaling[i],
                                              ) for i in range(self.number_atomic_orbitals)]
        

class Molecule():
    
    def __init__(self,
                 atoms_list: List[str],
                 coordinates: List[List[float]],
                 number_electrons: int,
                 ) -> None:
        
        self.atoms_list = atoms_list
        self.coordinates = coordinates
        self.number_atoms = len(self.atoms_list)
        self.atoms = [Atom(atom_type=self.atoms_list[i],
                           coordinates=coordinates[i],
                           ) for i in range(self.number_atoms)]

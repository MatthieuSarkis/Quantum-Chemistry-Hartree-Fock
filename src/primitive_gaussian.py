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
from src.utilities import boys, double_factorial, euclidean_distance_squared

class PrimitiveGaussian():
    """
    weight * normalization * (x-center_x)^lx * (y-center_y)^ly * (z-center_z)^lz * \exp(- \alpha |r - center|^2)
    """
    
    def __init__(self,
                 center: List[float],
                 alpha: float,
                 weight: float, 
                 zeta_weight: float = 1.0,
                 normalization: float = None,
                 ) -> None:
        
        self.center = np.array(center)
        self.alpha = alpha * zeta_weight**2
        self.weight = weight
        
        if normalization is None:
            self.normalization = (2.0 * self.alpha / np.pi)**(3/4)
        else:
            self.normalization = normalization
        
    def __mul__(self, 
                g: 'PrimitiveGaussian') -> 'PrimitiveGaussian':
        """Overload of the multiplication operator.
        
        The total normalizing factor comes from the normalizing factor of 
        both primitive gaussians, and from completing the square.
        """
        
        alpha = self.alpha + g.alpha
        diff = euclidean_distance_squared(self.center, g.center)
        N = self.normalization * g.normalization
        normalization = N * np.exp(-(self.alpha * g.alpha / alpha) * diff)
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
        prefactor = (np.pi / g.alpha)**(3/2)

        return prefactor * g.normalization
    
    @staticmethod
    def kinetic(g1: 'PrimitiveGaussian',
                g2: 'PrimitiveGaussian') -> float:
        
        g = g1 * g2
        prefactor = (np.pi / g.alpha)**(3/2)
        reduced_exponent = g1.alpha * g2.alpha / g.alpha

        return g.normalization * prefactor * reduced_exponent * (3 - 2 * reduced_exponent * euclidean_distance_squared(g1.center, g2.center))
    
    @staticmethod
    def potential(g1: 'PrimitiveGaussian',
                  g2: 'PrimitiveGaussian',
                  atom: object) -> float:
        
        g = g1 * g2
        R = atom.coordinates
        Z = atom.charge
    
        return (-2 * math.pi * Z / g.alpha) * g.normalization * boys(g.alpha * euclidean_distance_squared(g.center, R))
    
    @staticmethod
    def multi(g1: 'PrimitiveGaussian',
              g2: 'PrimitiveGaussian',
              g3: 'PrimitiveGaussian',
              g4: 'PrimitiveGaussian',
              ) -> float:
        
        ga = g1 * g2
        gb = g3 * g4
        prefactor = 2 * math.pi**2.5 / (ga.alpha * gb.alpha * math.sqrt(ga.alpha + gb.alpha))
    
        return prefactor * ga.normalization * gb.normalization * boys(ga.alpha * gb.alpha / (ga.alpha + gb.alpha) * euclidean_distance_squared(ga.center, gb.center))
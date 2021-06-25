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

from typing import List

from src.primitive_gaussian import PrimitiveGaussian

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
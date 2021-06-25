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
from src.atomic_orbital import AtomicOrbital 
from src.info import * 
class Atom():
    
    def __init__(self,
                 atom_type: str,
                 coordinates: List[float] = [0.0, 0.0, 0.0],
                 basis: str = 'STO-3G',
                 ) -> None:
        
        self.atom_type = atom_type
        self.charge = charge[self.atom_type]
        self.zeta_scaling = zeta[self.atom_type]
        self.coordinates = coordinates
        self.number_atomic_orbitals = maximun_number_orbitals[self.atom_type]
        
        if basis == 'STOnG':
            self.basis_weights = weights_STOnG[atom_type]
            self.basis_exponents = exponents_STOnG[atom_type]
        elif basis == 'STO-3G':
            self.basis_weights = weights_STO_3G[atom_type]
            self.basis_exponents = exponents_STO_3G[atom_type]
        elif basis == '6-31G':
            self.basis_weights = weights_6_31G[atom_type]
            self.basis_exponents = exponents_6_31G[atom_type]
            
        self.atomic_orbitals = [AtomicOrbital(basis_weights=self.basis_weights[i],
                                              basis_exponents=self.basis_exponents[i],
                                              center=self.coordinates,
                                              zeta_weight=self.zeta_scaling[i],
                                              ) for i in range(self.number_atomic_orbitals)]
        

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
class Atom():
    
    zeta = {'H' : [1.24],
            'He': [2.0925],
            'Li': [2.69, 0.80]}

    maximal_orbital_quantum_number = {'H' : 1,
                                      'He': 1,
                                      'Li': 2}
    
    charge = {'H' : 1,
              'He': 2,
              'Li': 3}
    
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
        

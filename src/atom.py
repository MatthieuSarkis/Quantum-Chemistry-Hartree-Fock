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
    
    zeta = {'H' : [1.24], 'He': [2.0925], 'Li': [2.69, 0.80]}
    maximal_orbital_quantum_number = {'H' : 1, 'He': 1, 'Li': 2}
    charge = {'H' : 1, 'He': 2, 'Li': 3}
    
    ### STOnG basis ###
    weights_STOnG = {'H' :[[0.444635, 0.535328, 0.154329], [0.700115, 0.399513, -0.0999672]]}
    exponents_STOnG = {'H' : [[0.109818, 0.405771, 2.22766], [0.0751386, 0.231031, 0.994203]]} 
    ### STO-3G basis ###
    weights_STO_3G = {'H' :[[0.1543289673E+00, 0.5353281423E+00, 0.4446345422E+00]]}
    exponents_STO_3G = {'H' : [[0.3425250914E+01, 0.6239137298E+00, 0.1688554040E+00]]}
    ### 6-31G basis ###
    weights_6_31G = {'H' :[[0.3349460434E-01, 0.2347269535E+00, 0.8137573261E+00], [1.0000000]]}
    exponents_6_31G = {'H' : [[0.1873113696E+02, 0.2825394365E+01, 0.6401216923E+00], [0.1612777588E+00]]}
    
    def __init__(self,
                 atom_type: str,
                 coordinates: List[float] = [0.0, 0.0, 0.0],
                 basis: str = 'STO-3G',
                 ) -> None:
        
        self.atom_type = atom_type
        self.charge = Atom.charge[self.atom_type]
        self.zeta_scaling = Atom.zeta[self.atom_type]
        self.coordinates = coordinates
        self.number_atomic_orbitals = len(self.zeta_scaling)
        
        if basis == 'STOnG':
            self.basis_weights = Atom.weights_STOnG[atom_type]
            self.basis_exponents = Atom.exponents_STOnG[atom_type]
        elif basis == 'STO-3G':
            self.basis_weights = Atom.weights_STO_3G[atom_type]
            self.basis_exponents = Atom.exponents_STO_3G[atom_type]
        elif basis == '6-31G':
            self.basis_weights = Atom.weights_6_31G[atom_type]
            self.basis_exponents = Atom.exponents_6_31G[atom_type]
            
        self.atomic_orbitals = [AtomicOrbital(basis_weights=self.basis_weights[i],
                                              basis_exponents=self.basis_exponents[i],
                                              center=self.coordinates,
                                              zeta_weight=self.zeta_scaling[i],
                                              ) for i in range(self.number_atomic_orbitals)]
        

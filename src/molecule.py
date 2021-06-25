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

from src.atom import Atom

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

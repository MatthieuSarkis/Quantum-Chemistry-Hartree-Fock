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
from scipy.special import erf

from src.atom import Molecule

def xyz_to_molecule(file_name: str) -> Molecule:
    
    atoms_list = []
    coordinates = []

    with open(file_name, 'r') as file:
        for idx, line in enumerate(file):
            
            if idx == 0 or idx == 1:
                continue

            split = line.split()
            atom_type = split[0]
            coordinates = [float(split[1]),
                           float(split[2]),
                           float(split[3])]
            
            atoms_list.append(atom_type)
            coordinates.append(coordinates)

    return Molecule(atoms_list=atoms_list, coordinates=coordinates)

def euclidean_distance_squared(R1: np.array,
                               R2: np.array,
                               ) -> float:
    
    return np.linalg.norm(np.array(R1) - np.array(R2))**2
    
def frobenius_norm(A: np.array) -> float:
    
    N = A.shape[0]
    norm_squared = 0
    
    for i in range(N):
        for j in range(N):
            norm_squared += A[i, j]**2 
        
    return np.sqrt(norm_squared)

def boys(x: float) -> float:
    
    if x == 0:
        return 1.0
    else:
        return math.sqrt(math.pi / x) * erf(math.sqrt(x)) / 2
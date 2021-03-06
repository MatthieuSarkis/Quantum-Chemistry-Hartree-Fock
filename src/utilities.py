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
from typing import List, Tuple


def xyz_to_molecule(file_path: str) -> Tuple[List[str], List[List[float]]]:
    
    atoms_list = []
    coordinates = []

    with open(file_path, 'r') as file:
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

    return atoms_list, coordinates


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


def double_factorial(n: int) -> int:
    
    from functools import reduce
    return reduce(int.__mul__, range(n, 0, -2))


def boys(x: float) -> float:
    
    if x == 0:
        return 1.0
    else:
        return math.sqrt(math.pi / x) * erf(math.sqrt(x)) / 2
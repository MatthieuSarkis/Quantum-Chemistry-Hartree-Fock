import math
import numpy as np
from scipy.special import erf
from typing import List

from src.atom import Molecule

def xyz_reader(file_name: str) -> Molecule:
    
    atoms_list = []
    coordinates = []

    with open('file_name', 'r') as file:
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

def euclidean_distance2(R1: np.array,
                        R2: np.array,
                        ) -> float:
    
    return np.linalg.norm(np.array(R1) - np.array(R2))**2

def boys(x: float) -> float:
    
    if x == 0:
        return 1.0
    else:
        return (0.5 * math.sqrt((math.pi / x))) * erf(math.sqrt(x))
import math
import numpy as np
from scipy.special import erf
from typing import List

def xyz_reader(file_name):
    file = open(file_name, 'r')

    number_of_atoms = 0
    atom_type = []
    atom_coordinates = []

    for idx, line in enumerate(file):
        if idx == 0:
            try:
                number_of_atoms = int(line.split()[0])
            except:
                print("file not in correct format.")
        
        if idx == 1:
            continue

        if idx != 0:
            split = line.split()
            atom = split[0]
            coordinates = [float(split[1]),
                           float(split[2]),
                           float(split[3])]
            atom_type.append(atom)
            atom_coordinates.append(coordinates)

    file.close()

    return number_of_atoms, atom_type, atom_coordinates


def euclidean_distance2(R1: List[float],
                       R2: List[float],
                       ) -> float:
    
    return np.linalg.norm(np.array(R1) - np.array(R2))**2

def boys(x: float) -> float:
    
    if x == 0:
        return 1.0
    else:
        return (0.5 * math.sqrt((math.pi / x))) * erf(math.sqrt(x))
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

from argparse import ArgumentParser
from src.hartree_fock import HartreeFock
from src.utilities import xyz_to_molecule

def main(args):
    
    molecule = xyz_to_molecule(args.xyz_file)
    
    hf = HartreeFock(molecule=molecule,
                     number_electrons=args.number_electrons)
    
    hf.compute_Hamiltonian_contributions()
    hf.self_consistent_field(epsilon=10e-4)
    

if __name__ == '__main__':
    
    parser = ArgumentParser()

    parser.add_argument('-xyz_file', type=str, help='')
    parser.add_argument('--number_electrons', type=int, default=2, help='')

    args = parser.parse_args()
    main(args)

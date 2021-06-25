from argparse import ArgumentParser

from src.hartree_fock import HartreeFock
from src.utilities import xyz_reader

def main(args):
    
    molecule = xyz_reader(args.xyz_file)
    
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

# Hartree-Fock algorithm

This repository contains an implementation of the restricted Hartree-Fock algorithm in python, that I started as a step in my Quantum Chemistry study.

We restrict ourselves in a first step to s orbitals, i.e. l=0 orbital angular momentum quantum number, and will extend to p, d, f, etc. orbitals in a second step.

Relevant links:
* https://www.basissetexchange.org
* https://www.mathematica-journal.com/?s=Evaluation+of+Gaussian+Molecular+Integrals+I&x=0&y=0
## Requirements

* numpy
* scipy

```shell
pip install -r requirements.txt
python setup.py install
```
 ## Examples 
 
```shell
python main.py --xyz_file data/xyz/HeH.xyz --basis_name STO-3G --number_electrons 2
```

## License
[Apache License 2.0](https://github.com/MatthieuSarkis/Quantum-Chemistry-Hartree-Fock/blob/master/LICENSE)

## To do

* include higher orbitals
* more comments
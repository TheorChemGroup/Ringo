# Ringo –⁠ The Python Library for Kinematics-Driven Conformer Generation

`Ringo` is a Python library that uses inverse kinematics to analyze the conformational flexibility of (poly)cyclic molecules by identifying independently rotatable dihedral angles of the molecule, and generate conformers when these values are chosen. It provides a back-end for identification and manipulation of degrees of freedom of cyclic molecules, i.e. setting their values to generate corresponding conformations. `Ringo`'s algorithm processes one set of dihedral angle values in fractions of a millisecond and scales well with number of rings in the molecule, allowing for efficient and comprehensive conformational searches of polycyclic molecules.

Try it out by installing (Linux x86_64 only, Python >= 3.9):

```bash
pip install --upgrade ringo-ik
```

The PyPi release also includes Monte-Carlo sampling algorithm (MCR) based on inverse kinematics and a grid sampling algorithm. See [docs](https://knvvv.gitlab.io/ringo) for details on these. Source code is coming soon.

Usage example:

```python
import ringo
import random
import math
import numpy as np

NUM_TRIES = 10000

# Load and analyze molecular topology
mol = ringo.Molecule(sdf='my_molecule.sdf')
p = ringo.Confpool() # Initialize conformer storage

# Simple Monte-Carlo conformational sampling
dofs_list, dofs_values = mol.get_ps()
for i in range(NUM_TRIES):
    for i, cur_param in enumerate(dofs_list):
        newvalue = random.uniform(-math.pi, math.pi)
        dofs_values[i] = newvalue

    result = mol.prepare_solution_iterator()
    if result != 0:
        # If not successfull, then try another set of dihedral values
        continue

    # Request the list of all inverse kinematics solutions
    # that passed overlap checks
    sol_list: list[np.ndarray] = mol.get_solutions_list()
    for i, matr in enumerate(sol_list):
        p.include_from_xyz(matr, f"Conformation #{len(p) + 1}")
    p.atom_symbols = mol.get_symbols()

assert len(p) > 0, f'No conformers generated in {NUM_TRIES} trials'
p.save_xyz('result.xyz')
```

## Documentation

[![pipeline status](https://gitlab.com/knvvv/ringo/badges/master/pipeline.svg)](https://gitlab.com/knvvv/ringo/-/commits/master)

Available [here](https://knvvv.gitlab.io/ringo)


## Links and references

[PyPi page](https://pypi.org/project/ringo-ik)

Papers on inverse kinematics conformer generation using Ringo:

Krivoshchapov, N. V.; Medvedev, M. G. Ring Kinematics-Informed Conformation Space Exploration. WIREs Computational Molecular Science 2024, 14 (1), e1690. https://doi.org/10.1002/wcms.1690

Krivoshchapov, N. V.; Medvedev, M. G. Accurate and Efficient Conformer Sampling of Cyclic Drug-Like Molecules with Inverse Kinematics. J. Chem. Inf. Model. 2024, 64 (11), 4542–4552. https://doi.org/10.1021/acs.jcim.3c02040

Ringo is based on inverse kinematics solver called TLC:

Coutsias EA, Seok C, Jacobson MP, Dill KA. A kinematic view of loop closure. J Comput Chem. 2004 Mar;25(4):510-28. doi: https://doi.org/10.1002/jcc.10416

### Contact information

* Nikolai V. Krivoshchapov – nikolai.krivoshchapov@gmail.com
* Michael G. Medvedev – medvedev.m.g@gmail.com

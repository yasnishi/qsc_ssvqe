## A quantum surrogate circuit of the subspace-search variational quantum eigensolver
This repository contains the source code used to produce the results presented in "A quantum surrogate circuit of the subspace-search variational quantum eigensolver". The source code can reproduce the inferred potential energy surfaces for the ground and excited states of H2 molecule shown in Figure 5 of the above-mentioned paper.

## Requirements

The following packages is needed:

```bash
pip install qiskit
pip install pyscf
pip install openfermion
pip install tqdm
```

We have confirmed that this code works with qiskit==0.44.1.

## How to use

To use the scripts, simply set the input data in notebooks. However, only n=4 (number of qubits) trained quantum surrogate circuit (QSC) with D=6 (PQC depth) are available. Optimized parameters are stored in the opt_params directory. When reading them, one can change the code as follows

```python
filename = './opt_params/singlet.pickle'  # doublet.pickle, triplet.pickle
```

## Citaion

If you are doing any research using this source code, please cite the following paper:

> Y. Nishida and A. Fumihiko, "A surrogate circuit of the subspace-search variational quantum eigensolver", Journal_name, xxxxx (20xx).


## License

The source code is licensed MIT.



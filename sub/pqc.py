'''
  pqc

'''

import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Pauli

# ansatz
from sub.ansatz import pqc_ansatz_HE

# utils
from sub.utils import expval


def QNNcircuit0(num_qubits, X, reps_pqc, reps, params):
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    # rescaling of X
    _X = 2.0 * (X[0] - X_min) / (X_max - X_min) + 0.5
    sigma_x = [_X, 0, 0]

    # encode bond length of H2
    for i in range(num_qubits):
        qc.h(i)
        qc.ry(sigma_x[0], i)

    block_params_num = num_qubits * reps_pqc

    # 1st blovk layer of HQCNN
    pqc = pqc_ansatz_HE(qc, params[0:int(block_params_num)], reps_pqc)
    op0 = Pauli('IIIZ')
    op1 = Pauli('IIZI')
    op2 = Pauli('IZII')
    op3 = Pauli('ZIII')
    observable = [op0, op1, op2, op3]

    # measurement of <\pi Z> by pqc
    sigma_z = expval(observable, [pqc for i in range(len(observable))]) * np.pi

    # 2nd block layer of HQCNN
    qc2 = QuantumCircuit(num_qubits)

    # input <\pi Z> measured from 1st layer
    for i in range(num_qubits):
        qc2.h(i)
        qc2.ry(sigma_z[i], i)

    pqc2 = pqc_ansatz_HE(
        qc2, params[int(block_params_num):int(block_params_num * 2)], reps_pqc)

    return pqc2


def QNNcircuit1(num_qubits, X, reps_pqc, reps, params):
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    # normalization of X
    _X = 2.0 * (X[0] - X_min) / (X_max - X_min) + 0.5
    sigma_x = [_X, 0, 0]

    # encode bond length of H2
    for i in range(num_qubits):
        qc.h(i)
        qc.ry(sigma_x[0], i)

    block_params_num = num_qubits * reps_pqc

    # 1st layer of HQCNN
    pqc = pqc_ansatz_HE(qc, params[0:int(block_params_num)], reps_pqc)
    op0 = Pauli('IIIZ')
    op1 = Pauli('IIZI')
    op2 = Pauli('IZII')
    op3 = Pauli('ZIII')
    observable = [op0, op1, op2, op3]

    # measurement of <\pi Z> by pqc
    sigma_z = expval(observable, [pqc for i in range(len(observable))]) * np.pi    

    # 2nd layer of HQCNN
    qc2 = QuantumCircuit(num_qubits)

    # 1st excited state for SSVQE
    qc2.x(0)

    # input <\pi Z> measured from 1st layer
    for i in range(num_qubits):
        qc2.h(i)
        qc2.ry(sigma_z[i], i)

    pqc2 = pqc_ansatz_HE(
        qc2, params[int(block_params_num):int(block_params_num * 2)], reps_pqc)

    return pqc2


def QNNcircuit2(num_qubits, X, reps_pqc, reps, params):
    X_min = 0.45
    X_max = 3.0

    qc = QuantumCircuit(num_qubits)

    # normalization of X
    _X = 2.0 * (X[0] - X_min) / (X_max - X_min) + 0.5
    sigma_x = [_X, 0, 0]

    # encode bond length of H2
    for i in range(num_qubits):
        qc.h(i)
        qc.ry(sigma_x[0], i)

    block_params_num = num_qubits * reps_pqc

    # 1st block layer of HQCNN
    pqc = pqc_ansatz_HE(qc, params[0:int(block_params_num)], reps_pqc)
    op0 = Pauli('IIIZ')
    op1 = Pauli('IIZI')
    op2 = Pauli('IZII')
    op3 = Pauli('ZIII')
    observable = [op0, op1, op2, op3]

    # measurement of <\pi Z> by pqc
    sigma_z = expval(observable, [pqc for i in range(len(observable))]) * np.pi
    
    # 2nd block layer of HQCNN
    qc2 = QuantumCircuit(num_qubits)

    # 2nd excited state for SSVQE
    qc2.x(1)

    # input <\pi Z> measured from 1st layer
    for i in range(num_qubits):
        qc2.h(i)
        qc2.ry(sigma_z[i], i)

    pqc2 = pqc_ansatz_HE(
        qc2, params[int(block_params_num):int(block_params_num * 2)], reps_pqc)

    return pqc2

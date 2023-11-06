'''
  ansatz
'''


def pqc_ansatz_HE(qc, params, reps_pqc):
    num_qubits = qc.num_qubits

    for i in range(0, reps_pqc):
        for j in range(num_qubits - 1):
            if j % 2 == 0:
                qc.cx(j, j + 1)
            else:
                pass
        for j in range(num_qubits - 1):
            if j % 2 == 1:
                qc.cx(j, j + 1)
            else:
                pass
        for j in range(num_qubits):
            qc.ry(params[qc.num_qubits * i + j], j)

        qc.barrier()

    return qc

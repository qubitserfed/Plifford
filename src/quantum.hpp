#include "linear_algebra.hpp"

namespace gate_matrices {
    Matrix<Qj> Pauli_X = Matrix<Qj>({
        {0, 1},
        {1, 0}
    });

    Matrix<Qj> Pauli_Y = Matrix<Qj>({
        {0, -I},
        {I, 0}
    });

    Matrix<Qj> Pauli_Z = Matrix<Qj>({
        {1, 0},
        {0, -1}
    });

    Matrix<Qj> Pauli_I = Matrix<Qj>({
        {1, 0},
        {0, 1}
    });

    Matrix<Qj> H = Matrix<Qj>({
        {1, 1},
        {1, -1}
    }) / sqrt2;

    Matrix<Qj> S = Matrix<Qj>({
        {1, 0},
        {0, I}
    });

    Matrix<Qj> CNOT = Matrix<Qj>({
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1},
        {0, 0, 1, 0}
    });

    Matrix<Qj> SWAP = Matrix<Qj>({
        {1, 0, 0, 0},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {0, 0, 0, 1}
    });

    Matrix<Qj> CZ = Matrix<Qj>({
        {1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, -1}
    });
} // namespace gates

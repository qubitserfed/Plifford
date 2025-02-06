#include "pauli.hpp"

#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

Pauli::Pauli(int _N) {
    N = _N;
    paulis = std::vector<PauliGate>(_N, PauliGate::I);
    phase = 1;
}

Pauli::Pauli(int _N, std::vector<PauliGate> _paulis, Qj _phase = 1) {
    N = _N;
    paulis = _paulis;
    phase = _phase;
}

Pauli::Pauli(int _N, PauliGate _type, int pos, Qj _phase = 1) {
    N = _N;
    paulis = std::vector<PauliGate>(_N, PauliGate::I);
    paulis[pos] = _type;
    phase = _phase;
}

Pauli X(int n, int pos) {
    return Pauli(n, PauliGate::X, pos);
}

Pauli Y(int n, int pos) {
    return Pauli(n, PauliGate::Y, pos);
}

Pauli Z(int n, int pos) {
    return Pauli(n, PauliGate::Z, pos);
}

std::pair<PauliGate, Qj> base_mul(PauliGate a, PauliGate b) {
    if (a == PauliGate::I)
        return {b, 1};

    if (b == PauliGate::I)
        return {a, 1}; 

    if (a == b)
        return {PauliGate::I, 1};

    if (a == PauliGate::X) {
        if (b == PauliGate::Y)
            return {PauliGate::Z, 1};
        else // b == PauliGate::Z
            return {PauliGate::Y, -1};
    }

    if (a == PauliGate::Y) {
        if (b == PauliGate::Z)
            return {PauliGate::X, 1};
        else // b == PauliGate::X
            return {PauliGate::Z, -1};
    }
}

Pauli operator * (Pauli lhs, Pauli rhs) {
    if (lhs.N != rhs.N)
        throw std::runtime_error("Pauli operators have different sizes");

    std::vector<PauliGate> new_paulis(lhs.N);
    Qj new_phase = lhs.phase * rhs.phase;

    for (int i = 0; i < lhs.N; i++) {
        auto [new_gate, new_sign] = base_mul(lhs.paulis[i], rhs.paulis[i]);
        new_paulis[i] = new_gate;
        new_phase = new_phase * new_sign;
    }

    return Pauli(lhs.N, new_paulis, new_phase);
}

Pauli operator * (Qj lhs, Pauli rhs) {
    return Pauli(rhs.N, rhs.paulis, lhs * rhs.phase);
}

Pauli operator * (Pauli lhs, Qj rhs) {
    return rhs * lhs;
}

Matrix<Qj> Pauli::matrix(const Pauli& p) {
    const int N = p.N;
    const int matN = 1 << N;

    Matrix<Qj> X = Matrix<Qj>({
        {0, 1},
        {1, 0}
    });

    Matrix<Qj> Y = Matrix<Qj>({
        {0, -I},
        {I, 0}
    });

    Matrix<Qj> Z = Matrix<Qj>({
        {1, 0},
        {0, -1}
    });

    Matrix<Qj> I2 = Matrix<Qj>({
        {1, 0},
        {0, 1}
    });

    std::vector< Matrix<Qj> > small_paulis(N);
    for (int i = 0; i < N; ++i) {
        switch (p.paulis[i]) {
            case PauliGate::I:
                small_paulis[i] = I2;
                break;
            case PauliGate::X:
                small_paulis[i] = X;
                break;
            case PauliGate::Y:
                small_paulis[i] = Y;
                break;
            case PauliGate::Z:
                small_paulis[i] = Z;
                break;
        }
    }
    Matrix<Qj> res = identity<Qj>(matN);
    for (int i = 0; i < matN; i++) {
        for (int j = 0; j < matN; ++j) {
            Qj entry = 1;

            for (int bit = 0; bit < N; ++bit) {
                entry = entry * small_paulis[bit].get(i >> bit & 1, j >> bit & 1);
            }
        }
    }

    return res * p.phase;
}

void print(const Pauli& p) {
    print(p.phase);
    for (int i = 0; i < p.N; i++) {
        switch (p.paulis[i]) {
            case PauliGate::I:
                std::cout << "I";
                break;
            case PauliGate::X:
                std::cout << "X";
                break;
            case PauliGate::Y:
                std::cout << "Y";
                break;
            case PauliGate::Z:
                std::cout << "Z";
                break;
        }
    }
}

// in place operators

void Pauli::mulX(int pos) {
    paulis[pos] = PauliGate::X;
}

void Pauli::mulY(int pos) {
    paulis[pos] = PauliGate::Y;
}

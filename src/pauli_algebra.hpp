#pragma once

#include "linear_algebra.hpp"
#include "pauli.hpp"

#include <vector>

struct PauliVect {
    int N;
    std::vector<Qj> paulis;

    Matrix<Qj> matrix();

    PauliVect &operator += (const PauliVect &rhs);
    PauliVect &operator *= (const PauliVect &rhs);

    PauliVect(int _N);
    PauliVect(const Pauli &pauli);
    PauliVect(const Matrix<Qj> &m);
};

PauliVect operator * (const PauliVect& lhs, const PauliVect& rhs);
//PauliVect operator * (const Pauli& lhs, const PauliVect& rhs); // might not be necessary due to the existence of the constructor
//PauliVect operator * (const PauliVect& lhs, const Pauli &rhs);
PauliVect operator * (const PauliVect& lhs, const Qj& rhs);
PauliVect operator * (const Qj& lhs, const PauliVect& rhs);

PauliVect operator + (const PauliVect &lhs, const PauliVect &rhs);
PauliVect operator - (const PauliVect &lhs, const PauliVect &rhs);
PauliVect operator - (const PauliVect &p);

void print(const PauliVect &p);


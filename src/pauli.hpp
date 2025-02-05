#pragma once

#include "linear_algebra.hpp"

#include <vector>


enum class PauliGate {I, X, Y, Z};

struct Pauli {
    int N;

    std::vector<PauliGate> paulis;
    Qj phase;
    
    void mulX(int pos);
    void mulY(int pos);
    void mulZ(int pos);

    Pauli(int _N);
    Pauli(int _N, PauliGate _type, int pos, Qj _phase=1);
    Pauli(int _N, std::vector<PauliGate> _paulis, Qj _phase=1);

    Matrix<Qj> matrix(const Pauli& p);
};

Pauli X(int n, int pos);
Pauli Y(int n, int pos);
Pauli Z(int n, int pos);

Pauli operator * (Pauli lhs, Pauli rhs);
Pauli operator * (Qj lhs, Pauli rhs);
Pauli operator * (Pauli lhs, Qj rhs);

void print(const Pauli& p);

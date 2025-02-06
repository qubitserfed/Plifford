#include "pauli_algebra.hpp"
#include "pauli.hpp"

#include <bitset>
#include <iostream>
#include <vector>

// utility functions:

int pauli_to_index(const Pauli &p) {
    int index = 0;

    for (int i = 0; i < p.N; ++i) {
        int digit;
        switch ( p.paulis[i] ) {
            case PauliGate::I:
                digit = 0;
                break;
            case PauliGate::X:
                digit = 1;
                break;
            case PauliGate::Y:
                digit = 3;
                break;
            case PauliGate::Z:
                digit = 2;
                break;
        }
        index+= digit << 2 * i;
    }
    
    return index;
}

Pauli index_to_pauli(int index, int N) {
    Pauli p(N);
    for (int i = 0; i < N; ++i) {
        int digit = (index >> 2 * i) & 3;
        switch ( digit ) {
            case 0:
                p.paulis[i] = PauliGate::I;
                break;
            case 1:
                p.paulis[i] = PauliGate::X;
                break;
            case 2:
                p.paulis[i] = PauliGate::Z;
                break;
            case 3:
                p.paulis[i] = PauliGate::Y;
                break;
        }
    }
    return p;
}

// PauliVect constructors:

PauliVect::PauliVect(int _N) {
    N = _N;
    paulis.resize(1 << 2 * N);
}

PauliVect::PauliVect(const Pauli &p) {
    N = p.N;
    paulis.resize(1 << 2 * N);

    paulis[pauli_to_index(p)] = 1;
}

PauliVect operator * (const PauliVect &lhs, const PauliVect &rhs) { // Test, for the love of God
    std::vector<int> non_zero_lhs, non_zero_rhs;
    for (int i = 0; i < lhs.paulis.size(); ++i) {
        if (lhs.paulis[i] != 0) {
            non_zero_lhs.push_back(i);
        }
        if (rhs.paulis[i] != 0) {
            non_zero_rhs.push_back(i);
        }
    }

    PauliVect result(lhs.N);
    for (int i = 0; i < non_zero_lhs.size(); ++i) {
        for (int j = 0; j < non_zero_rhs.size(); ++j) {
            int fst_idx = non_zero_lhs[i];
            int snd_idx = non_zero_rhs[j];
            int index = non_zero_lhs[i] ^ non_zero_rhs[j];
            // need to include the Pauli commutation factor

            const int even_mask = 0xAAAAAAAA;
            const int odd_mask = 0x55555555;
            int commutation_factor =
                //               X positions on lhs    * Z positions on rhs
                std::bitset<32>( fst_idx & even_mask & ((snd_idx & odd_mask) >> 1) )
                .count() % 2;

            result.paulis[index] = result.paulis[index] + lhs.paulis[non_zero_lhs[i]] * rhs.paulis[non_zero_rhs[j]];
        }
    }
}

PauliVect operator * (const PauliVect &lhs, const Qj &rhs) {
    PauliVect result(lhs.N);
    for (int i = 0; i < lhs.paulis.size(); ++i) {
        result.paulis[i] = lhs.paulis[i] * rhs;
    }
    return result;
}

PauliVect operator * (const Qj &lhs, const PauliVect &rhs) {
    return rhs * lhs;
}

PauliVect operator + (const PauliVect &lhs, const PauliVect &rhs) {
    PauliVect result(lhs.N);
    for (int i = 0; i < lhs.paulis.size(); ++i) {
        result.paulis[i] = lhs.paulis[i] + rhs.paulis[i];
    }
    return result;
}

void print(const PauliVect &p) {
    for (int i = 0; i < p.paulis.size(); ++i) {
        if (p.paulis[i] != 0) {
            print(p.paulis[i]);
            std::cout << " ";
            print(index_to_pauli(i, p.N));
            if (i != p.paulis.size() - 1) {
                std::cout << " + ";
            }
        }
    }
}

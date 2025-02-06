#pragma once

#include "rings.hpp"

#include <algorithm>
#include <vector>

typedef long long i64;
typedef unsigned long long u64;

// General Vector Definition

template <typename T>
struct Vector {
    std::vector<T> vec;

    T get(int);
    bool is_zero();
    void set(int, T);
    int weight();

    Vector<T>& operator += (const Vector<T> &);
    Vector<T>& operator -= (const Vector<T> &);
    Vector<T>& operator *= (const T &);
    Vector<T>& operator /= (const T &);
    T  operator () (int);

    Vector();
    Vector(int);
    Vector(std::vector<T>);
};

// General Matrix Definition

template <typename T>
struct Matrix {
    int n, m;
    std::vector< std::vector<T> > vec;

    Vector<T> &row(int);
    Vector<T> &last_row();

    void append_row(Vector<T>);
    void pop_row();
    void set(int, int, T);
    T get(int, int);
    bool empty();

    void add_rows(int, int);
    void add_rows(int, int, T);
    void swap_rows(int, int);
    void swap_cols(int, int);
    void remove_zero_rows();
    void sort_rows();

    std::vector<Vector<T>> rows();
    std::vector<Vector<T>> columns();

    Matrix<T>& operator += (const Matrix<T> &);
    Matrix<T>& operator -= (const Matrix<T> &);
    Matrix<T>& operator *= (const T &);
    Matrix<T>& operator /= (const T &);
    Matrix<T>& operator *= (const Matrix<T> &);
    Matrix<T>& operator += (const T &);
    Matrix<T>& operator -= (const T &);
    Matrix<T>& operator /= (const T &);
    T  operator () (int, int);

    Matrix();
    Matrix(int);
    Matrix(int, int);
    Matrix(const std::vector< std::vector<T> > &);
    Matrix(const std::vector<Vector<T>> &);
    Matrix(const Vector<T> &);
};


// Utility Functions
int                                     popcount                        (u64 num);

// General Vector Operators
template<typename T> void               print                           (Vector<T> vec);
template<typename T> Vector<T>          operator +                      (const Vector<T> &, const Vector<T> &);
template<typename T> T                  operator *                      (const Vector<T> &, const Vector<T> &);
template<typename T> bool               operator <                      (const Vector<T> &, const Vector<T> &);
template<typename T> bool               operator ==                     (const Vector<T> &, const Vector<T> &);
template<typename T> bool               operator !=                     (const Vector<T> &, const Vector<T> &);
template<typename T> Vector<T>          kroeneker                       (Vector<T>, Vector<T>);

// General Matrix Operators
template<typename T> void               print                           (Matrix<T>);
//template<typename T> void               print_latex                     (Matrix<T>);
template<typename T> Matrix<T>          transpose                       (Matrix<T> &);
template<typename T> Vector<T>          operator *                      (Vector<T>, Matrix<T>);
template<typename T> Matrix<T>          operator *                      (Matrix<T>, Matrix<T>);
template<typename T> Matrix<T>          operator +                      (Matrix<T>, Matrix<T>);
template<typename T> Matrix<T>          identity                        (const int &);
template<typename T> Matrix<T>          zero                            (const int &, const int &);
template<typename T> Matrix<T>          kroeneker                       (Matrix<T>, Matrix<T>);
template<typename T> Matrix<T>          operator *                      (const T &, Matrix<T>);
template<typename T> Matrix<T>          operator *                      (Matrix<T>, const T &);
template<typename T> Matrix<T>          operator /                      (Matrix<T>, const T &);

// Specialized Matrix Operators
Matrix<Qj>                              adjoint                         (const Matrix<Qj> &);

// Linear Algebraic Algorithms
template<typename T> int                to_row_echelon                  (Matrix<T> &);
template<typename T> std::vector<int>   restricted_row_echelon          (Matrix<T> &, std::vector<int>);
template<typename T> Vector<T>          canonical_quotient              (Vector<T>, Matrix<T> &);
template<typename T> bool               in_span                         (Matrix<T>, const Vector<T> &);

template<typename T> int                rank                            (Matrix<T>);
template<typename T> T                  det                             (Matrix<T>);
template<typename T> bool               is_invertible                   (Matrix<T>);
template<typename T> bool               is_singular                     (Matrix<T>);
template<typename T> Matrix<T>          inv                             (Matrix<T>);

template<typename T> Matrix<T>          transposed_product              (const Matrix<T> &, const Matrix<T> &);
template<typename T> Vector<T>          transposed_product              (const Matrix<T> &, const Vector<T> &);

template<typename T> Matrix<T>          basis_completion                (Matrix<T>);
template<typename T> T                  sym_prod                        (Vector<T>, Vector<T>);
template<typename T> Matrix<T>          symplectic_complement           (const Matrix<T> &);
//template<typename T> Matrix<T>          isotropic_completion            (const Matrix<T> &);

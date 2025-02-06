#include "linear_algebra.hpp"
#include "rings.hpp"

#include <algorithm>
#include <bitset>
#include <vector>

#define SAFE


// Utility Functions
int popcount(u64 arg) {
    std::bitset<64>(arg).count();
}

// General Vector Methods

template <typename T>
T Vector<T>::get(int idx) {
#ifdef SAFE
    if (idx >= vec.size())
        throw std::runtime_error("Index out of bounds");
#endif
    return vec[idx];
}

template <typename T>
void Vector<T>::set(int idx, T val) {
    vec[idx] = val;
}

template <typename T>
bool Vector<T>::is_zero() {
    for (int i = 0; i < vec.size(); i++)
        if (vec[i] != 0)
            return false;
    return true;
}

template <typename T>
int Vector<T>::weight() {
    int res = 0;
    for (int i = 0; i < vec.size(); i++)
        if (vec[i] != 0)
            res++;
    return res;
}

template <typename T>
Vector<T> &Vector<T>::operator += (const Vector<T> &a) {
    *this = *this + a;
    return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator -= (const Vector<T> &a) {
    *this = *this - a;
    return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator *= (const T &a) {
    *this = *this * a;
    return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator /= (const T &a) {
    *this = *this / a;
    return *this;
}

template <typename T>
T Vector<T>::operator () (int idx) {
    return vec[idx];
}

// constructor
template <typename T>
Vector<T>::Vector() {
    vec = std::vector<T>();
}

template <typename T>
Vector<T>::Vector(int n) {
    vec = std::vector<T>(n);
}

template <typename T>
Vector<T>::Vector(std::vector<T> v) {
    vec = v;
}

// General Matrix Methods

template <typename T>
Vector<T> &Matrix<T>::row(int idx) {
#ifdef SAFE
    if (idx >= n)
        throw std::runtime_error("Index out of bounds");
#endif
    return vec[idx];
}

template <typename T>
Vector<T> &Matrix<T>::last_row() {
    return vec[n - 1];
}

template <typename T>
void Matrix<T>::append_row(Vector<T> v) {
#ifdef SAFE
    if (v.vec.size() != m)
        throw std::runtime_error("Inconsistent row size in append_row");
#endif
    vec.push_back(v);
    n++;
}

template <typename T>
void Matrix<T>::pop_row() {
#ifdef SAFE
    if (n == 0)
        throw std::runtime_error("Cannot pop from empty matrix");
#endif
    vec.pop_back();
    n--;
}

template <typename T>
void Matrix<T>::set(int i, int j, T val) {
#ifdef SAFE
    if (i >= n || j >= m)
        throw std::runtime_error("Index out of bounds");
#endif
    vec[i].set(j, val);
}

template <typename T>
T Matrix<T>::get(int i, int j) {
#ifdef SAFE
    if (i >= n || j >= m)
        throw std::runtime_error("Index out of bounds");
#endif
    return vec[i].get(j);
}

template <typename T>
bool Matrix<T>::empty() {
    return n == 0;
}

template <typename T>
void Matrix<T>::add_rows(int i, int j) {
#ifdef SAFE
    if (i >= n || j >= n)
        throw std::runtime_error("Index out of bounds");
    if (vec[i].vec.size() != vec[j].vec.size())
        throw std::runtime_error("Inconsistent row size in add_rows");
#endif
    for (int k = 0; k < m; k++)
        vec[i].set(k, vec[i].get(k) + vec[j].get(k));
}


template <typename T>
void Matrix<T>::add_rows(int i, int j, T lambda) {
#ifdef SAFE
    if (i >= n || j >= n)
        throw std::runtime_error("Index out of bounds");
    if (vec[i].vec.size() != vec[j].vec.size())
        throw std::runtime_error("Inconsistent row size in add_rows");
#endif
    for (int k = 0; k < m; k++)
        vec[i].set(k, vec[i].get(k) * lambda + vec[j].get(k));
}

template <typename T>
void Matrix<T>::swap_rows(int i, int j) {
#ifdef SAFE
    if (i >= n || j >= n)
        throw std::runtime_error("Index out of bounds");
#endif
    std::swap(vec[i], vec[j]);
}

template <typename T>
void Matrix<T>::swap_cols(int i, int j) {
#ifdef SAFE
    if (i >= m || j >= m)
        throw std::runtime_error("Index out of bounds");
#endif
    for (int k = 0; k < n; k++) {
        T ki = vec[k].get(i);
        T kj = vec[k].get(j);
        vec[k].set(i, kj);
        vec[k].set(j, ki);
    }
}

template <typename T>
void Matrix<T>::remove_zero_rows() {
    for (int i = 0; i < n; i++) {
        if (vec[i].is_zero()) {
            swap_rows(i, n - 1);
            pop_row();
            i--;
        }
    }
}

template <typename T>
void Matrix<T>::sort_rows() {
    std::sort(vec.begin(), vec.end(), [](Vector<T> &a, Vector<T> &b) {
        return a.vec < b.vec;
    });
}

template <typename T>
std::vector<Vector<T>> Matrix<T>::rows() {
    std::vector<Vector<T>> res;
    for (int i = 0; i < n; i++)
        res.emplace_back(vec[i]);
    return res;
}

template <typename T>
std::vector<Vector<T>> Matrix<T>::columns() {
    std::vector<Vector<T>> res(m, Vector<T>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            res[j].set(i, vec[i].get(j));
    return res;
}

template <typename T>
Matrix<T> &Matrix<T>::operator += (const Matrix<T> &a) {
    *this = *this + a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator -= (const Matrix<T> &a) {
    *this = *this - a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *= (const T &a) {
    *this = *this * a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator /= (const T &a) {
    *this = *this / a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *= (const Matrix<T> &a) {
    *this = *this * a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator += (const T &a) {
    *this = *this + a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator -= (const T &a) {
    *this = *this - a;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator /= (const T &a) {
    *this = *this / a;
    return *this;
}

template <typename T>
T Matrix<T>::operator () (int i, int j) {
    return vec[i](j);
}

// Matrix constructors:

template <typename T>
Matrix<T>::Matrix() {
    n = m = 0;
    vec = std::vector< Vector<T> >();
}

template <typename T>
Matrix<T>::Matrix(int n, int m) {
    this->n = n;
    this->m = m;
    vec = std::vector< Vector<T> >(n, Vector<T>(m));
}

template <typename T>
Matrix<T>::Matrix(int n) {
    this->n = n;
    this->m = n;
    vec = std::vector< Vector<T> >(n, Vector<T>(n));
}

template <typename T>
Matrix<T>::Matrix(const std::vector< std::vector<T> > &v) {
#ifdef SAFE
    for (int i = 0; i < v.size(); i++)
        if (v[i].size() != v[0].size())
            throw std::runtime_error("Inconsistent row size in Matrix constructor");
#endif

    n = v.size();
    m = v[0].size();
    vec = std::vector< Vector<T> >(n);
    for (int i = 0; i < n; i++)
        vec[i] = Vector<T>(v[i]);
}

template <typename T>
Matrix<T>::Matrix(const std::vector<Vector<T>> &v) {
#ifdef SAFE
    for (int i = 0; i < v.size(); i++)
        if (v[i].size() != v[0].size())
            throw std::runtime_error("Inconsistent row size in Matrix constructor");
#endif
    n = v.size();
    m = v[0].vec.size();
    vec.resize(n);
    for (int i = 0; i < n; i++)
        vec[i] = v[i];
}

template <typename T>
Matrix<T>::Matrix(const Vector<T> &v) {
    n = 1;
    m = v.vec.size();
    vec = std::vector< Vector<T> >(1, v);
}

// General Vector Operators

template <typename T>
void print(Vector<T> v) {
    for (int i = 0; i < v.vec.size(); i++)
        std::cout << v.vec[i] << " ";
    std::cout << std::endl;
}

template <typename T>
Vector<T> operator + (const Vector<T> &a, const Vector<T> &b) {
#ifdef SAFE
    if (a.vec.size() != b.vec.size())
        throw std::runtime_error("Inconsistent vector size in operator +");
#endif
    Vector<T> res(a.vec.size());
    for (int i = 0; i < a.vec.size(); i++)
        res.vec[i] = a.vec[i] + b.vec[i];
    return res;
}

template <typename T>
T operator * (const Vector<T> &a, const Vector<T> &b) {
#ifdef SAFE
    if (a.vec.size() != b.vec.size())
        throw std::runtime_error("Inconsistent vector size in operator *");
#endif
    T res = 0;
    for (int i = 0; i < a.vec.size(); i++)
        res += a.vec[i] * b.vec[i];
    return res;
}

template <typename T>
bool operator < (const Vector<T> &a, const Vector<T> &b) {
    return a.vec < b.vec;
}

template <typename T>
bool operator == (const Vector<T> &a, const Vector<T> &b) {
    return a.vec == b.vec;
}

template <typename T>
bool operator != (const Vector<T> &a, const Vector<T> &b) {
    return a.vec != b.vec;
}

// General Matrix Operators

template <typename T>
void print(Matrix<T> mat) {
    for (int i = 0; i < mat.n; i++) {
        for (int j = 0; j < mat.m; j++)
            std::cout << mat.get(i, j) << " ";
        std::cout << std::endl;
    }
}

template <typename T>
Matrix<T> transpose(Matrix<T> &mat) {
    Matrix<T> res(mat.m, mat.n);
    for (int i = 0; i < mat.n; i++)
        for (int j = 0; j < mat.m; j++)
            res.set(j, i, mat.get(i, j));
    return res;
}

template <typename T>
Vector<T> operator * (Vector<T> a, Matrix<T> mat) {
#ifdef SAFE
    if (a.vec.size() != mat.n)
        throw std::runtime_error("Inconsistent vector size in operator *");
#endif
    Vector<T> res(mat.m);
    for (int i = 0; i < mat.m; i++)
        for (int j = 0; j < mat.n; j++)
            res.vec[i] += a.vec[j] * mat.get(j, i);
    return res;
}

template <typename T>
Matrix<T> operator * (Matrix<T> a, Matrix<T> b) {
#ifdef SAFE
    if (a.m != b.n)
        throw std::runtime_error("Inconsistent matrix size in operator *");
#endif
    Matrix<T> res(a.n, b.m);
    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < b.m; j++)
            for (int k = 0; k < a.m; k++)
                res.set(i, j, res.get(i, j) + a.get(i, k) * b.get(k, j));
    return res;
}

template <typename T>
Matrix<T> operator + (Matrix<T> a, Matrix<T> b) {
#ifdef SAFE
    if (a.n != b.n || a.m != b.m)
        throw std::runtime_error("Inconsistent matrix size in operator +");
#endif
    Matrix<T> res(a.n, a.m);
    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < a.m; j++)
            res.set(i, j, a.get(i, j) + b.get(i, j));
    return res;
}

template <typename T>
Matrix<T> identity(const int &n) {
    Matrix<T> res(n, n);
    for (int i = 0; i < n; i++)
        res.set(i, i, 1);
    return res;
}

template <typename T>
Matrix<T> zero(const int &n, const int &m) {
    return Matrix<T>(n, m);
}

template <typename T>
Matrix<T> kroeneker(Matrix<T> a, Matrix<T> b) {
    Matrix<T> res(a.n * b.n, a.m * b.m);

    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < a.m; j++)
            for (int k = 0; k < b.n; k++)
                for (int l = 0; l < b.m; l++)
                    res.set(i * b.n + k, j * b.m + l, a.get(i, j) * b.get(k, l));

    return res;
}

template <typename T>
Vector<T> kroeneker(Vector<T> a, Vector<T> b) {
    Vector<T> res(a.vec.size() * b.vec.size());

    for (int i = 0; i < a.vec.size(); i++)
        for (int j = 0; j < b.vec.size(); j++)
            res.vec[i * b.vec.size() + j] = a.vec[i] * b.vec[j];

    return res;
}

template <typename T>
Matrix<T> operator * (const T &a, Matrix<T> mat) {
    Matrix<T> res(mat.n, mat.m);
    for (int i = 0; i < mat.n; i++)
        for (int j = 0; j < mat.m; j++)
            res.set(i, j, a * mat.get(i, j));
    return res;
}

template <typename T>
Matrix<T> operator * (Matrix<T> mat, const T &a) {
    return a * mat;
}

template <typename T>
Matrix<T> operator / (Matrix<T> mat, const T &a) {
    return a / mat;
}

Matrix<Qj> adjoint(Matrix<Qj> mat) {
    Matrix<Qj> res = transpose(mat);
    const int N = res.n;
    const int M = res.m;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j)
            res.set(i, j, conj(res.get(i, j), 7));
    return res;
}

// Matrix algorithms


// puts the argument matrix in row echelon form and returns the rank of the matrix
template <typename T>
int to_row_echelon(Matrix<T> &mat) {
    int rank = 0;
    for (int i = 0; i < mat.m; i++) {
        int pivot = -1;
        for (int j = rank; j < mat.n; j++) {
            if (mat.get(j, i) != 0) {
                pivot = j;
                break;
            }
        }
        if (pivot == -1)
            continue;
        mat.swap_rows(rank, pivot);
        for (int j = rank + 1; j < mat.n; j++) {
            T lambda = mat.get(j, i) / mat.get(rank, i);
            mat.add_rows(j, rank, -lambda);
        }
        rank++;
    }
    return rank;
}

template <typename T>
std::vector<int> restricted_row_echelon(Matrix<T> &mat, std::vector<int> columns) {

}

template <typename T>
Vector<T> canonical_quotient(Vector<T> v, Matrix<T> &mat) {

}

template <typename T>
bool in_span (Matrix<T> mat, const Vector<T> &v) {
    return canonical_quotient(v, mat).is_zero();
}

template <typename T>
int rank(Matrix<T> mat) {
    return to_row_echelon(mat);
}

template <typename T>
int det(Matrix<T> mat) {
    if (mat.n != mat.m)
        throw std::runtime_error("Determinant of non-square matrix");
    T res = 1;
    mat.to_row_echelon();
    for (int i = 0; i < mat.n; i++)
        res *= mat.get(i, i);
    return res;
}

template <typename T>
bool is_invertible(Matrix<T> mat) {
    return mat.n == mat.m && rank(mat) == mat.n;
}

template <typename T>
bool is_singular(Matrix<T> mat) {
    return !is_invertible(mat);
}

template <typename T>
Matrix<T> inv(Matrix<T> mat) {
    if (mat.n != mat.m)
        throw std::runtime_error("Taking inverse of non-square matrix");

    const int N = mat.n;
    Matrix<T> tmat(N, 2 * N), res(N, N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            tmat.set(i, j, mat.get(i, j));
        tmat.set(i, i + N, 1);
    }

    to_row_echelon(tmat);
    for (int i = 0; i < N; ++i) {
        T lambda = tmat.get(i, i);
        if (lambda == 0)
            throw std::runtime_error("Matrix is singular");
        for (int j = 0; j < N; ++j)
            tmat.set(i, j, tmat.get(i, j) / lambda);
    }

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            res.set(i, j, mat.get(i, j + N));
    return res;
}

// computes A*B^T
template <typename T>
Matrix<T> transposed_product(const Matrix<T> &a, const Matrix<T> &b) {
    Matrix res(a.n, b.n);
    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < b.n; j++)
            for (int k = 0; k < a.m; k++)
                res.set(i, j, res.get(i, j) + a.get(i, k) * b.get(j, k));
    return res;
}


// computes A*v^T
template <typename T>
Vector<T> transposed_product(const Matrix<T> &a, const Vector<T> &b) { // same as a*b tho
    Vector<T> res(a.n);
    for (int i = 0; i < a.n; i++)
        for (int j = 0; j < a.m; j++)
            res.vec[i] += a.get(i, j) * b.vec[j];
    return res;
}

template<typename T>
Matrix<T> basis_completion(Matrix<T> mat) {
    const int N = mat.n;
    const int M = mat.m;
    int rank = to_row_echelon();

    Matrix<T> res(N, M - rank);
    std::set<int> pivot_columns;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            if (mat.get(i, j) != 0) {
                pivot_columns.insert(j);
                break;
            }
        }
    }

    assert( pivot_columns.size() == rank );
    
    int ptr = 0;
    for (int j = 0; j < M; ++j) {
        if (pivot_columns.find(j) == pivot_columns.end()) {
            res.set(ptr, j, 1);
            ++ptr;
        }
    }

    return res;
}

template<typename T>
T sym_prod(Vector<T> a, Vector<T> b) {
    if (a.vec.size() != b.vec.size())
        throw std::runtime_error("Inconsistent vector size in sym_prod");
    if (a.vec.size() % 2 == 1)
        throw std::runtime_error("Symmetric product of odd size vectors");

    T res = 0;
    for (int i = 0; i < a.vec.size(); i+= 2)
        res += a.vec[i] * b.vec[i + 1] - a.vec[i + 1] * b.vec[i];
    return res;
}

/*
    Given an isotropic subspace A of T^n spanned by the rows of a matrix, this function
    returns a matrix whose rows span a subspace B of T^n such that A + B is the symplectic orthogonal complement of A
    and A and B intersect trivially. In particular, if the rows of the input matrix are the stabiliizers of some
    stabilizer code, the rows of the output matrix will span the logical operators of the code.
*/
//template<typename T>
//Matrix<T> isotropic_completion(Matrix<T> basis) {
//    std::vector<Vector<T>> rows = basis.rows();
//}

template<typename T>
Matrix<T> symplectic_complement(const Matrix<T> &basis_mat) {
    const int N = basis.m;
    if (N % 2 == 1)
        throw std::runtime_error("Symplectic complement of odd size subspace");

    std::vector< Vector<T> > basis = basis_mat.rows();

    Matrix<T> res = identity(N);

    for (int i = 0; i < res.size(); ++i) {
        bool commuting = true;
        for (const auto &bv: basis) {
            if (sym_prod(res.row(i), bv) != 0) {
                for (int j = i + 1; j < res.size(); ++j) {
                    if (sym_prod(res.row(j), bv) != 0)
                        res.add_rows(i, j);
                }
                commuting = false;
                break;
            }
            if (!commuting){
                res.swap_cols(i, N - 1);
                res.pop_row();
                i-= 1;
            }
        }
    }

    return res;
}

#include "linear_algebra.hpp"
#include "rings.hpp"

#include <algorithm>
#include <bitset>
#include <vector>

#define SAFE

// General Vector Methods

template <typename T>
T Vector<T>::get(int idx) {
#ifdef SAFE
    if (idx >= vec.size())
        throw "Index out of bounds";
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
}

template <typename T>
Vector<T> &Vector<T>::operator -= (const Vector<T> &a) {
    *this = *this - a;
}

template <typename T>
Vector<T> &Vector<T>::operator *= (const T &a) {
    *this = *this * a;
}

template <typename T>
Vector<T> &Vector<T>::operator /= (const T &a) {
    *this = *this / a;
}

template <typename T>
Vector<T> Vector<T>::operator () (int idx) {
    return Vector<T>(vec[idx]);
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
        throw "Index out of bounds";
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
        throw "Inconsistent row size in append_row";
#endif
    vec.push_back(v);
    n++;
}

template <typename T>
void Matrix<T>::pop_row() {
#ifdef SAFE
    if (n == 0)
        throw "Cannot pop from empty matrix";
#endif
    vec.pop_back();
    n--;
}

template <typename T>
void Matrix<T>::set(int i, int j, T val) {
#ifdef SAFE
    if (i >= n || j >= m)
        throw "Index out of bounds";
#endif
    vec[i].set(j, val);
}

template <typename T>
T Matrix<T>::get(int i, int j) {
#ifdef SAFE
    if (i >= n || j >= m)
        throw "Index out of bounds";
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
        throw "Index out of bounds";
    if (vec[i].vec.size() != vec[j].vec.size())
        throw "Inconsistent row size in add_rows";
#endif
    for (int k = 0; k < m; k++)
        vec[i].set(k, vec[i].get(k) + vec[j].get(k));
}

template <typename T>
void Matrix<T>::swap_rows(int i, int j) {
#ifdef SAFE
    if (i >= n || j >= n)
        throw "Index out of bounds";
#endif
    std::swap(vec[i], vec[j]);
}

template <typename T>
void Matrix<T>::swap_cols(int i, int j) {
#ifdef SAFE
    if (i >= m || j >= m)
        throw "Index out of bounds";
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
Matrix<T> &Matrix<T>::operator += (const Matrix<T> &a) {
    *this = *this + a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator -= (const Matrix<T> &a) {
    *this = *this - a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *= (const T &a) {
    *this = *this * a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator /= (const T &a) {
    *this = *this / a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator *= (const Matrix<T> &a) {
    *this = *this * a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator += (const T &a) {
    *this = *this + a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator -= (const T &a) {
    *this = *this - a;
}

template <typename T>
Matrix<T> &Matrix<T>::operator /= (const T &a) {
    *this = *this / a;
}

template <typename T>
Vector<T> Matrix<T>::operator () (int i, int j) {
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

//Matrix(std::vector< std::vector<T> >);
//Matrix(std::vector<Vector<T>>);

template <typename T>
Matrix<T>::Matrix(std::vector< std::vector<T> > v) {
#ifdef SAFE
    for (int i = 0; i < v.size(); i++)
        if (v[i].size() != v[0].size())
            throw "Inconsistent row size in Matrix constructor";
#endif

    n = v.size();
    m = v[0].size();
    vec = std::vector< Vector<T> >(n);
    for (int i = 0; i < n; i++)
        vec[i] = Vector<T>(v[i]);
}

template <typename T>
Matrix<T>::Matrix(std::vector<Vector<T>> v) {
#ifdef SAFE
    for (int i = 0; i < v.size(); i++)
        if (v[i].size() != v[0].size())
            throw "Inconsistent row size in Matrix constructor";
#endif
    n = v.size();
    m = v[0].vec.size();
    vec.resize(n);
    for (int i = 0; i < n; i++)
        vec[i] = v[i];
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
        throw "Inconsistent vector size in operator +";
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
        throw "Inconsistent vector size in operator *";
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
        throw "Inconsistent vector size in operator *";
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
        throw "Inconsistent matrix size in operator *";
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
        throw "Inconsistent matrix size in operator +";
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

/*

template<typename T> Matrix<T>          operator *                      (const T &, Matrix<T>);
template<typename T> Matrix<T>          operator *                      (Matrix<T>, const T &);
template<typename T> Matrix<T>          operator /                      (const T &, Matrix<T>);
template<typename T> Matrix<T>          operator /                      (Matrix<T>, const T &);
*/

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


// Utility Functions
int popcount(u64 arg) {
    std::bitset<64>(arg).count();
}

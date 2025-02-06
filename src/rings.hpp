typedef long long i64;
typedef unsigned long long u64;

// Ring and Fields definitions

struct GF2 {
    bool val;

    GF2 &operator += (GF2);
    GF2 &operator -= (GF2);
    GF2 &operator *= (GF2);
    GF2 &operator /= (GF2);

    GF2();
    GF2(i64 a);
    GF2(bool a);
};

struct GF4  {
    bool e1, e2; // represents: e1 + e2 * \alpha where \alpha is the multiplicative group generator of GF4

    GF4 &operator += (GF4);
    GF4 &operator -= (GF4);
    GF4 &operator *= (GF4);
    GF4 &operator /= (GF4);

    GF4();
    GF4(i64 a);
    GF4(GF2 a);
    GF4(bool a, bool b);
};

struct Q {
    i64 num, den;

    Q &operator += (Q);
    Q &operator -= (Q);
    Q &operator *= (Q);
    Q &operator /= (Q);

    Q();
    Q(i64 a);
    Q(i64 a, i64 b); // a/b
};

struct Qj { // Field extension of Q by j (where j is an elementary 8th root of unity)
    Q projs[4]; // represents: projs[0] * 1 + projs[1] * j + projs[2] * j^2 + projs[3] * j^3 + projs[4] * j^4

    Qj &operator += (Qj);
    Qj &operator -= (Qj);
    Qj &operator *= (Qj);
    Qj &operator /= (Qj);

    Qj();
    Qj(i64 a);
    Qj(Q a);
    Qj(Q a, Q b, Q c, Q d);
};

struct Z4 {
    char val;
};


// GF2 Operators
GF2     operator +  (GF2, GF2);
GF2     operator -  (GF2, GF2);
GF2     operator -  (GF2);
GF2     operator *  (GF2, GF2);
GF2     inv         (GF2);
GF2     operator /  (GF2, GF2);
bool    operator == (GF2, GF2);
bool    operator != (GF2, GF2);
bool    operator <  (GF2, GF2);
void    print       (GF2);
void    print_latex (GF2);

// GF4 Operators
GF4     operator +  (GF4, GF4);
GF4     operator -  (GF4, GF4);
GF4     operator -  (GF4);
GF4     operator *  (GF4, GF4);
GF4     inv         (GF4);
GF4     operator /  (GF4, GF4);
bool    operator == (GF4, GF4);
bool    operator != (GF4, GF4);
bool    operator <  (GF4, GF4);
void    print       (GF4);
void    print_latex (GF4);

// Z4 Operators
Z4      operator +  (Z4, Z4);
Z4      operator -  (Z4, Z4);
Z4      operator -  (Z4);
Z4      operator *  (Z4, Z4);
bool    operator == (Z4, Z4);
bool    operator != (Z4, Z4);
bool    operator <  (Z4, Z4);
void    print       (Z4);
void    print_latex (Z4);
// Q Operators
Q       operator +  (Q, Q);
Q       operator -  (Q, Q);
Q       operator -  (Q);
Q       conj        (Q, int e=7); // applies to first argument the automorphism mapping j to j^e, defaults to complex conjugate
Q       operator *  (Q, Q);
Q       inv         (Q);
Q       operator /  (Q, Q);
bool    operator == (Q, Q);
bool    operator != (Q, Q);
bool    operator <  (Q, Q);
bool    operator <= (Q, Q);
bool    operator >= (Q, Q);
bool    operator >  (Q, Q);
void    print       (Q);
void    print_latex (Q);
// Qj Operators
Qj      operator +  (Qj, Qj);
Qj      operator -  (Qj, Qj);
Qj      operator -  (Qj);
Qj      operator *  (Qj, Qj);
Qj      operator *  (Q, Qj);
Qj      operator *  (Qj, Q);
Qj      operator /  (Qj, Q);
Qj      conj        (Qj, int e); // applies to first argument the automorphism mapping j to j^e
Q       norm        (Qj);
Q       trace       (Qj);
Qj      inv         (Qj);
Qj      operator /  (Qj, Qj);
bool    operator == (Qj, Qj);
bool    operator != (Qj, Qj);
bool    operator <  (Qj, Qj);
void    print       (Qj);
void    print_latex (Qj);

const Qj J = Qj(Q(0, 1), Q(1, 1), Q(0, 1), Q(0, 1));
const Qj I = J * J;
const Qj sqrt2 = (I - Qj(1)) * J;


// J = (I + 1) / sqrt 2

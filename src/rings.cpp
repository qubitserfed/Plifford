#include "rings.hpp"

#include <format>
#include <iostream>
#include <string>

// utility functions
i64 gcd(i64 a, i64 b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b) {
        i64 t = b;
        b = a % b;
        a = t;
    }
    return a;
}

// Positive characteristic rings -------------
// GF2 constructors
GF2::GF2() {
    val = false;
}

GF2::GF2(i64 a) {
    val = (a % 2 != 0);
}

GF2::GF2(bool a) {
    val = a;
}

// GF2 operators

GF2 operator + (GF2 a, GF2 b) {
    return GF2 {a.val != b.val};
}

GF2 operator - (GF2 a, GF2 b) {
    return GF2 {a.val != b.val};
}

GF2 operator - (GF2 a) {
    return a;
}

GF2 operator * (GF2 a, GF2 b) {
    return GF2 {a.val && b.val};
}

GF2 inv(GF2 arg) {
    if (!arg.val)
        throw "Division by zero";
    return arg;
}

GF2 operator / (GF2 a, GF2 b) {
    if (!b.val) {
        throw "Division by zero";
    }
    return GF2 {a.val && b.val};
}

bool operator == (GF2 a, GF2 b) {
    return a.val == b.val;
}

bool operator != (GF2 a, GF2 b) {
    return a.val != b.val;
}

bool operator < (GF2 a, GF2 b) {
    return !a.val && b.val;
}

void print(GF2 a) {
    std::cout << a.val;
}

void print_latex(GF2 a) {
    std::cout << a.val << "_2";
}

// GF4 constructors
/*
    GF4();
    GF4(int a);
    GF4(GF2 a);
    GF4(bool a, bool b);
*/

GF4::GF4() {
    e1 = false;
    e2 = false;
}

GF4::GF4(i64 a) {
    e1 = (a % 2 != 0);
    e2 = false;
}

GF4::GF4(GF2 a) {
    e1 = a.val;
    e2 = false;
}

GF4::GF4(bool a, bool b) {
    e1 = a;
    e2 = b;
}

// GF4 Operators
GF4 operator + (GF4 a, GF4 b) {
    return GF4 {a.e1 != b.e1, a.e2 != b.e2};
}

GF4 operator - (GF4 a, GF4 b) {
    return GF4 {a.e1 != b.e1, a.e2 != b.e2};
}

GF4 operator - (GF4 a) {
    return a;
}

GF4 operator * (GF4 a, GF4 b) {
    return GF4 {
        (a.e1 != b.e1) != (a.e2 && b.e2), // 1 part
        (a.e1 && b.e2) != (a.e2 && b.e1)  // \alpha part
    };
}

GF4 inv(GF4 arg) {
    if (arg.e1 == 0 && arg.e2 == 0)
        throw "Division by zero";
    else if (arg.e1 == 1 && arg.e2 == 0)
        return GF4 {1, 0};
    else if (arg.e1 == 0 && arg.e2 == 1)
        return GF4 {1, 0};
    else if (arg.e1 == 1 && arg.e2 == 1)
        return GF4 {0, 1};
    else
        throw "Invalid GF4 element";
}

GF4 operator / (GF4 a, GF4 b) {
    return a * inv(b);
}

bool operator == (GF4 a, GF4 b) {
    return a.e1 == b.e1 && a.e2 == b.e2;
}

bool operator != (GF4 a, GF4 b) {
    return a.e1 != b.e1 || a.e2 != b.e2;
}

bool operator < (GF4 a, GF4 b) {
    if (a.e1 < b.e1)
        return true;
    else if (a.e1 == b.e1)
        return a.e2 < b.e2;
    else
        return false;
}

void print(GF4 a) {
    std::cout << '(' << a.e1 << ',' << a.e2 << ')';
}

void print_latex(GF4 a) {
    std::cout << a.e1 << " + " << a.e2 << "\\alpha";
}

// Z4 operators
Z4 operator + (Z4 a, Z4 b) {
    return Z4 {char( (a.val + b.val) % 4 ) }; // trivially optimizable
}

Z4 operator - (Z4 a, Z4 b) {
    return Z4 {char( (a.val - b.val + 4) % 4 ) }; // trivially optimizable
}

Z4 operator - (Z4 a) {
    return Z4 {char( (4 - a.val) % 4 ) }; // trivially optimizable
}

Z4 operator * (Z4 a, Z4 b) {
    return Z4 {char( (a.val * b.val) % 4 ) }; // trivially optimizable
}

bool operator == (Z4 a, Z4 b) {
    return a.val == b.val;
}

bool operator != (Z4 a, Z4 b) {
    return a.val != b.val;
}

bool operator < (Z4 a, Z4 b) {
    return a.val < b.val;
}

void print(Z4 a) {
    std::cout << a.val;
}

void print_latex(Z4 a) {
    std::cout << a.val << "_4";
}




// Zero characteristic fields -------------
// Q constructors
Q::Q() {
    num = 0;
    den = 1;
}

Q::Q(i64 a) {
    num = a;
    den = 1;
}

Q::Q(i64 a, i64 b) {
    if (b == 0)
        throw "Division by zero";
    i64 g = gcd(a, b);
    num = a / g;
    den = b / g;
    if (den < 0) {
        num = -num;
        den = -den;
    }
}
// Q operators
Q operator + (Q a, Q b) {
    i64 up = a.num * b.den + b.num * a.den;
    i64 down = a.den * b.den;
    i64 g = gcd(up, down);
    return Q {up / g, down / g};
}

Q operator - (Q a, Q b) {
    i64 up = a.num * b.den - b.num * a.den;
    i64 down = a.den * b.den;
    if (down < 0) {
        up = -up;
        down = -down;
    }
    i64 g = gcd(up, down);
    return Q {up / g, down / g};
}

Q operator - (Q a) {
    return Q {-a.num, a.den};
}

Q operator * (Q a, Q b) {
    i64 up = a.num * b.num;
    i64 down = a.den * b.den;
    i64 g = gcd(up, down);
    return Q {up / g, down / g};
}

Q inv(Q a) {
    if (a.num == 0)
        throw "Division by zero";

    Q tmp = Q {a.den, a.num};
    if (tmp.den < 0) {
        tmp.num = -tmp.num;
        tmp.den = -tmp.den;
    }
    return tmp;
}

Q operator / (Q a, Q b) {
    return a * inv(b);
}

bool operator == (Q a, Q b) {
    return a.num == b.num && a.den == b.den;
}

bool operator != (Q a, Q b) {
    return a.num != b.num || a.den != b.den;
}

bool operator < (Q a, Q b) {
    return a.num * b.den < b.num * a.den;
}

bool operator > (Q a, Q b) {
    return a.num * b.den > b.num * a.den;
}

bool operator <= (Q a, Q b) {
    return a.num * b.den <= b.num * a.den;
}

bool operator >= (Q a, Q b) {
    return a.num * b.den >= b.num * a.den;
}

void print(Q a) {
    std::cout << a.num << '/' << a.den;
}

void print_latex(Q a) {
    std::cout << "\\frac{" << a.num << "}{" << a.den << "}";
}


// Qj constructors

Qj::Qj() {
    for (int i = 0; i < 4; i++)
        projs[i] = Q();
}

Qj::Qj(Q a) {
    projs[0] = a;
    for (int i = 1; i < 4; i++)
        projs[i] = Q();
}

Qj::Qj(i64 a) {
    projs[0] = Q(a, 1);
    for (int i = 1; i < 4; i++)
        projs[i] = Q();
}

Qj::Qj(Q a, Q b, Q c, Q d) {
    projs[0] = a;
    projs[1] = b;
    projs[2] = c;
    projs[3] = d;
}

// Qj operators

Qj operator + (Qj a, Qj b) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = a.projs[i] + b.projs[i];
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator - (Qj a, Qj b) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = a.projs[i] - b.projs[i];
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator - (Qj a) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = -a.projs[i];
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator * (Qj a, Qj b) {
    Q res[4], temp[8];
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            temp[i + j] = a.projs[i] * b.projs[j];
    for (int i = 0; i < 4; i++)
        res[i] = temp[i] - temp[i + 4];
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator * (Q a, Qj b) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = a * b.projs[i];
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator * (Qj a, Q b) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = a.projs[i] * b;
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj operator / (Qj a, Q b) {
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = a.projs[i] / b;
    return Qj(res[0], res[1], res[2], res[3]);
}

Qj conj(Qj a, int e=7) {
    if (e < 0 || e > 4)
        throw "Invalid automorphism";
    Q res[4];
    for (int i = 0; i < 4; i++)
        res[i] = res[i] + a.projs[(i * e) % 4];
    return Qj(res[0], res[1], res[2], res[3]);
}
Q norm(Qj a) {
    Qj temp = a * conj(a, 3) * conj(a, 5) * conj(a, 7);
    return temp.projs[0];
}

Q trace(Qj a) {
    Qj temp = a + conj(a, 3) + conj(a, 5) + conj(a, 7);
    return temp.projs[0];
}

Qj inv(Qj a) {
    Q n = norm(a);
    if (n.num == 0)
        throw "Division by zero";
    return conj(a, 3) * conj(a, 5) * conj(a, 7) / n;
}

Qj operator / (Qj a, Qj b) {
    return a * inv(b);
}

bool operator == (Qj a, Qj b) {
    for (int i = 0; i < 4; i++)
        if (a.projs[i] != b.projs[i])
            return false;
    return true;
}

bool operator != (Qj a, Qj b) {
    return !(a == b);
}

bool operator < (Qj a, Qj b) {
    for (int i = 0; i < 4; i++)
        if (a.projs[i] < b.projs[i])
            return true;
        else if (a.projs[i] != b.projs[i])
            return false;
    return false;
}

void print(Qj a) {
    std::cout << '(';
    for (int i = 0; i < 4; i++) {
        print(a.projs[i]);
        if (i < 3)
            std::cout << ", ";
    }
    std::cout << ')';
}

void print_latex(Qj a) {
    std::cout << '(';
    for (int i = 0; i < 4; i++) {
        std::string term;
        std::string versor;
        switch (i) {
            case 0:
                versor = "";
                break;
            case 1:
                versor = "j^1";
                break;
            case 2:
                versor = "j^2";
                break;
            case 3:
                versor = "j^3";
                break;
        }

        if (a.projs[i].num == 0) {
            term = "0";
        } else if (a.projs[i].den == 1) {
            term = std::to_string(a.projs[i].num) + versor;
        } else {
            term = std::format("\\frac{{ {} {} }}{{{}}}{}", a.projs[i].num, versor, a.projs[i].den);
        }

        if (i < 3)
            std::cout << " + ";
    }
    std::cout << ')';
}

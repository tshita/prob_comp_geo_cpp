[🏠 Home](../index.md)  |  [🔗 AOJ CGL_2_A](https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/all/CGL_2_A)

# 問題：平行・垂直（Parallel / Orthogonal）
2 直線 $l_1, l_2$ が与えられる。 $l_1$ と $l_2$ が平行、直交、またはそれ以外かを判定せよ。

# 解法
$l_1$ と $l_2$ が平行・垂直かどうかは各々の直線上を通る異なる任意の 2 点から構成される 2 つのベクトルが平行・垂直かどうかの問題に帰着される。

2 つのベクトル $\mathbf{v_1}, \mathbf{v_2}$ が平行かどうかは外積の絶対値 $| \mathbf{v_1} \times \mathbf{v_2}|$ で判定する（`abs_cross` 関数）。外積の絶対値は 2 つのベクトルがなす平行四辺形の符号付面積に等しいので、その値が 0 なら 2 つのベクトルが平行ということになる。  
平行かどうかは `is_parallel` 関数で判定する。

2 つのベクトルが垂直かどうかは内積で判定する（`dot` 関数）。内積は $\mathbf{v_1} \cdot \mathbf{v_2} = |\mathbf{v_1}| |\mathbf{v_2}| \cos \theta$ である。ただし、 $|\mathbf{v_1}|$ と $|\mathbf{v_2}|$ はそれぞれ $\mathbf{v_1}$ と $\mathbf{v_2}$ の大きさで、 $\theta$ は $\mathbf{v_1}, \mathbf{v_2}$ のなす角である。2 つのベクトルが直角の場合は $\cos \theta$ が $0$ となるので、内積が $0$ となる。  
垂直かどうかは `is_orthogonal` 関数で判定する。

---------------------------------------------------------------------------------------------

<details>
<summary>src/aoj/CGL_2_A.cc を表示</summary>

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <algorithm>
#include <cassert>
#include <vector>
#include <stack>

using Real = long double;

constexpr Real EPS = 1e-10;
const Real PI = acos(static_cast<Real>(-1.0)); // GCC 4.6.1 以上で acos() は constexpr の場合がある

inline int sign(Real a) { return (a < -EPS) ? -1 : (a > EPS) ? +1 : 0; }
inline bool eq(Real a, Real b)  { return sign(a - b) == 0; }  // a = b
inline bool neq(Real a, Real b) { return !eq(a, b); }         // a != b
inline bool lt(Real a, Real b)  { return sign(a - b) == -1; } // a < b
inline bool leq(Real a, Real b) { return sign(a - b) <= 0; }  // a <= b
inline bool gt(Real a, Real b)  { return sign(a - b) == 1; }  // a > b
inline bool geq(Real a, Real b) { return sign(a - b) >= 0; }  // a >= b

/**
 * Point in two dimensional
 */
struct Point2 {
    Real x{0.0}, y{0.0};

    explicit Point2() {}
    Point2(Real x, Real y) : x(x), y(y) {}

    // Arithmetic operator between Point2s
    Point2 operator+(const Point2 &rhs) const { return Point2(x + rhs.x, y + rhs.y); }
    Point2 operator-(const Point2 &rhs) const { return Point2(x - rhs.x, y - rhs.y); }
    Point2 operator*(const Point2 &rhs) const { // cross product between Point2s
        return Point2(x * rhs.x - y * rhs.y, x * rhs.y + y * rhs.x);
    }

    // Unary operator and compound assignment operator
    Point2 operator-() const { return {-x, -y}; }
    Point2& operator+=(const Point2 &rhs) { return *this = *this + rhs; }
    Point2& operator-=(const Point2 &rhs) { return *this = *this - rhs; }

    // Arithmetic operator between Point2 and Real
    Point2 operator*(Real rhs) const { return Point2(x * rhs, y * rhs); }
    Point2 operator/(Real rhs) const { return Point2(x / rhs, y / rhs); }

    // Comparison operation
    bool operator==(const Point2 &rhs) const { return eq(x, rhs.x) && eq(y, rhs.y); }
    bool operator!=(const Point2 &rhs) const { return neq(x, rhs.x) || neq(y, rhs.y); }
    bool operator<(const Point2 &rhs) const { return lt(x, rhs.x) || (eq(x, rhs.x) && lt(y, rhs.y)); }
    bool operator>(const Point2 &rhs) const { return gt(x, rhs.x) || (eq(x, rhs.x) && gt(y, rhs.y)); }

    // Other operator
    // ユークリッド距離を返す
    Real abs(void) const { return std::hypot(x, y); }

    // ユークリッド距離の二乗を返す
    Real abs2(void) const { return x * x + y * y; }

    // 単位はラジアンで範囲 [-PI, PI] で x 軸の正の方向となす角度を返す
    // atan2(y, x) は y / x の逆正接を返す（arctan(y / x)）
    // atan(z) と異なりどの象限に属しているか分かるので正しい符号を返す
    Real arg(void) const { return atan2(y, x); }

    // 内積
    Real dot(const Point2 &rhs) const { return x * rhs.x + y * rhs.y; }

    // 原点を中心に反時計回りに90度回転する
    Point2 rotate90(void) { return *this = Point2(-y, x); }

    // 原点を中心に反時計回りに angle [rad] だけ回転する
    void rotate(Real angle) {
            *this = Point2(cos(angle) * x - sin(angle) * y, sin(angle) * x + cos(angle) * y);
    }
};

Point2 operator*(Real a, Point2 p) { return p * a; }

// Output and input of a Point2
std::ostream& operator<<(std::ostream &os, const Point2 &p) { return os << p.x << ' ' << p.y; }
std::istream& operator>>(std::istream &is, Point2 &p) { return is >> p.x >> p.y; }

// ベクトル p1 と p2 の内積： dot(p1, p2) = |a| |b| cos(theta)
inline Real dot(const Point2 &p1, const Point2 &p2) { return p1.x * p2.x + p1.y * p2.y; }

// ベクトル p1 と p2 の外積の絶対値 |p1 x p2| ： |p1 x p2| = |p1| |p2| sin(theta)
// 原点, p1, p2 を頂点とする符号付き平行四辺形の面積（ p1 から p2 へ反時計回りで符号が正）
inline Real abs_cross(const Point2 &p1, const Point2 &p2) { return p1.x * p2.y - p1.y * p2.x; }

/**
 * Line in two dimensional
 */
class Line : public std::array<Point2, 2> {
public:
    explicit Line() {}
    Line(const Point2 &p1, const Point2 &p2) {
        (*this)[0] = p1;
        (*this)[1] = p2;
    }
};

// Input and output of a Line
std::istream& operator>>(std::istream &is, Line &l) { return is >> l[0] >> l[1]; }
std::ostream& operator<<(std::ostream &os, const Line &l) { return os << l[0] << ' ' << l[1]; }


// --------------------8<------- start of main part of library -------8<--------------------

inline bool is_parallel(const Line &l1, const Line &l2) {
    return eq(abs_cross(l1[0] - l1[1], l2[0] - l2[1]), 0.0);
}

inline bool is_orthogonal(const Line &l1, const Line &l2) {
    return eq(dot(l1[0] - l1[1], l2[0] - l2[1]), 0.0);
}

// --------------------8<------- end of main part of library   -------8<--------------------


int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    unsigned q;
    std::cin >> q;

    for (auto i = 0u; i < q; ++i) {
        Line line1, line2;
        std::cin >> line1 >> line2;

        if (is_parallel(line1, line2)) std::cout << 2 << std::endl;
        else if (is_orthogonal(line1, line2)) std::cout << 1 << std::endl;
        else std::cout << 0 << std::endl;
    }

    return 0;
}
```

</details>

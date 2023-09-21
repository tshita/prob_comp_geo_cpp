[🏠 Home](../../README.md)

# 射影（Projection）: [AOJ CGL_1_A](https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/1/CGL_1_A)
ユークリッド平面上で 3 点 $p_0, p_1, p$ が与えられる。 $p_1$ と $p_2$ を通る直線に対する $p$ の射影を求めよ。

# 解法
直線 $l$ に対する点 $p$ の射影とは、 $p$ から $l$ への垂線を引いた交点のことである。

`projection` 関数に直線 $l$ と点 $p$ を与えると射影が返される。  

$l$ を通る 2 点を $p_0, p_1$ 、求める射影を $q$ とする。このとき、直角三角形 $\triangle p_0 p q$ を考えて射影を導出する。
2 つのベクトル $\overrightarrow{p_0 p_1}, \overrightarrow{p_0 p}$ のなす角度を $\theta$ とすると、  

```math
| \overrightarrow{p_0 q} | = | \overrightarrow{p_0 p} | \, \cos{\theta} = | \overrightarrow{p_0 p} | \frac{\overrightarrow{p_0 p_1} \cdot \overrightarrow{p_0 p}}{|\overrightarrow{p_0 p_1}| |\overrightarrow{p_0 p}|}  = \frac{\overrightarrow{p_0 p_1} \cdot \overrightarrow{p_0 p}}{|\overrightarrow{p_0 p_1}|}
```

となる。ここで、 $\overrightarrow{p_0 p_1} \cdot \overrightarrow{p_0 p}$ はベクトル $\overrightarrow{p_0 p_1}$ と $\overrightarrow{p_0 p}$ の内積（`dot` 関数）で、 $|\overrightarrow{p_0 q}|$ はベクトル $|\overrightarrow{p_0 q}|$ のノルム（`Point2::abs` 関数）を表している。  以上で求める射影は、

```math
q = p_0 + |\overrightarrow{p_0 q}| \frac{\overrightarrow{p_0 p_1}}{|\overrightarrow{p_0 p_1}|} = p_0 + \frac{|\overrightarrow{p_0 p_1}| \, \overrightarrow{p_0 p_1} \cdot \overrightarrow{p_0 p}}{|\overrightarrow{p_0 p_1}|^2}
```

となる。あとはこの式の通りに実装すればよい。


# 参考文献
- [高校数学の美しい物語「正射影ベクトルの公式の証明と使い方」](https://manabitimes.jp/math/933) （最終アクセス 2023年9月21日）  
  実装後に見つけたが導出が同じでより分かりやすい説明のためおススメ。


----------------------------------------------------------------------------------------------------------------------------------------------

<details>
<summary>src/CGL_1_A.cc を表示</summary>

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <cassert>

using Real = double;

constexpr Real EPS = 1e-10;
const Real PI = acos(static_cast<Real>(-1.0)); // GCC 4.6.1 以上で acos() は constexpr の場合がある

inline int sign(Real a) { return (a < -EPS) ? -1 : (a > EPS) ? +1 : 0; }
inline bool eq(Real a, Real b)  { return sign(a - b) == 0; }  // a = b
inline bool neq(Real a, Real b) { return !eq(a, b); }         // a != b
inline bool lt(Real a, Real b)  { return sign(a - b) == -1; } // a < b
inline bool leq(Real a, Real b) { return sign(a - b) <= 0; }  // a <= b
inline bool gt(Real a, Real b)  { return sign(a - b) == 1; }  // a > b
inline bool geq(Real a, Real b) { return sign(a - b) >= 0; }  // a >= b

// change between degree and radian
inline Real to_radian(const Real degree) { return degree * PI / 180.0; }
inline Real to_degree(const Real radian) { return radian * 180.0 / PI; }

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

/**
 * Intersection testing
 */
Point2 projection(const Line &l, const Point2 &p) {
    Point2 dir = l[1] - l[0];
    Real t = dot(p - l[0], dir) / dir.abs2();
    return l[0] + dir * t;
}


int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    Line line;
    std::cin >> line;

    unsigned int q;
    std::cin >> q;

    for (unsigned int i = 0; i < q; ++i) {
        Point2 p;
        std::cin >> p;
        
        std::cout << projection(line, p) << "\n";
    }

    return 0;
}
```

</details>

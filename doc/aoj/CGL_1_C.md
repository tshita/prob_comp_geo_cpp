# 反時計回り: [AOJ CGL_1_C](https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/1/CGL_1_C)
ユークリッド平面上に 3 点 $p_0, p_1, p_2$ が与えられる。$p_0, p_1, p_2$ をこの順に訪れるとき $p_1$ でどの方向に折れるのかを判定せよ。

# 解法
`enum class CCW` で 3 点の向きの結果を表す。
- `CCW::COUNTER_CLOCKWISE`:  
    $p_1$ で反時計回りに折れる
- `CCW::CLOCKWISE`:  
    $p_1$ で時計回りに折れる
- `CCW::ONLINE_FRONT`:  
    $p_0, p_1, p_2$ の順番で同一直線上に乗る
- `CCW::ONLINE_BACK`:  
    $p_2, p_0, p_1$ の順番で同一直線状に乗る
- `CCW::ON_SEGMENT`:  
    「$p_0$ と $p_1$ を結ぶ線分上に $p_2$ が乗るか」、「$p_0$ と $p_1$ が等しくなく、かつ、 $p_1$ と $p_2$ が等しいか」、または、「3 点が等しいか」のいずれか
- `CCW::OTHER`:  
    それ以外の場合（そのような場合は無いので例外扱い）

判定は `ccw` 関数で行う。詳しくは関数を参照。


<details>
<summary>src/CGL_1_C.cc を表示</summary>

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <vector>

using Real = double;

constexpr Real EPS = 1e-8;
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
    Real abs(void) const { return std::hypot(x, y); } // ユークリッド距離を返す
    // ユークリッド距離の二乗を返す
    Real abs2(void) const { return x * x + y * y; }
    // 単位はラジアンで範囲 [-PI, PI] で x 軸の正の方向となす角度を返す
    // atan2(y, x) は y / x の逆正接を返す（arctan(y / x)）
    // atan(z) と異なりどの象限に属しているか分かるので正しい符号を返す
    Real arg(void) const { return atan2(y, x); }
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
// 原点, p1, p2 を頂点とする符号付き三角形の面積（ p1 から p2 へ反時計回りで符号が正）
inline Real abs_cross(const Point2 &p1, const Point2 &p2) { return p1.x * p2.y - p1.y * p2.x; }

// ベクトル p1 から p2 への角度を返す（単位はラジアン）
// p1 と p2 のなす角度で小さい方で p1 から p2 へ反時計回りなら符号は正，時計回りなら符号は負
inline Real arg(const Point2 &p1, const Point2 &p2) { return atan2(abs_cross(p1, p2), dot(p1, p2)); }

// Whether one object contains the other
enum class CONTAIN {
    IN,
    ON,
    OUT,
};

// Counter-Clockwise predicate (a, b, c)
enum class CCW : int {
    COUNTER_CLOCKWISE = 1,     // counter clockwise
    CLOCKWISE         = -1,    // clockwise
    ONLINE_FRONT      = 2,     // a--b--c on line or (a == b and b != c)
    ONLINE_BACK       = -2,    // c--a--b on line
    ON_SEGMENT        = 0,     // a--c--b on line or (a != b and b == c) or (a == b == c)
    OTHER             = -3,
};

CCW inv(const CCW dir) {
    switch (dir) {
        case (CCW::COUNTER_CLOCKWISE): return CCW::CLOCKWISE;
        case (CCW::CLOCKWISE): return CCW::COUNTER_CLOCKWISE;
        case (CCW::ONLINE_FRONT): return CCW::ONLINE_BACK;
        case (CCW::ONLINE_BACK): return CCW::ONLINE_FRONT;
        default: return dir;
    }
}

CCW ccw(const Point2 &a, Point2 b, Point2 c) {
    b -= a;  c -= a;
    if (sign(abs_cross(b, c)) == 1) return CCW::COUNTER_CLOCKWISE;
    if (sign(abs_cross(b, c)) == -1) return CCW::CLOCKWISE;
    if (sign(dot(b, c)) == -1)       return CCW::ONLINE_BACK;
    if (sign(b.abs2() - c.abs2()) == -1)   return CCW::ONLINE_FRONT;
    return CCW::ON_SEGMENT;
}

auto ccw_t(const Point2 &a, Point2 b, Point2 c) {
    return static_cast<std::underlying_type<CCW>::type>(ccw(a, std::move(b), std::move(c)));
}


/**
 * Distance and Intersection Point2
 */
inline Real distance(const Point2 &p1, const Point2 &p2) {
    return (p1 - p2).abs();
}


int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    Point2 p0, p1;
    unsigned int q;
    std::cin >> p0 >> p1 >> q;

    for (unsigned i = 0; i < q; ++i) {
        Point2 p2;
        std::cin >> p2;

        switch (ccw(p0, p1, p2)) {
            case CCW::COUNTER_CLOCKWISE:
                std::cout << "COUNTER_CLOCKWISE\n";
                break;
            case CCW::CLOCKWISE:
                std::cout << "CLOCKWISE\n";
                break;
            case CCW::ONLINE_BACK:
                std::cout << "ONLINE_BACK\n";
                break;
            case CCW::ONLINE_FRONT:
                std::cout << "ONLINE_FRONT\n";
                break;
            case CCW::ON_SEGMENT:
                std::cout << "ON_SEGMENT\n";
                break;
            case CCW::OTHER:
                exit(EXIT_FAILURE);
        }
    }

    return 0;
}
```

</details>

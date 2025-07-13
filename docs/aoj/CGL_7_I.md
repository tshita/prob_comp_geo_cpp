[🏠 Home](../index.md)  |  [🔗 AOJ CGL_7_I](https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/all/CGL_7_I)

# ２つ円の共通部分の面積（Area of Intersection between Two Circles）
2 つの円 $c_1, c_2$ の共通部分の面積を求めよ。

# 解法
2 つの円の共通部分の面積は `area_intersection` で求める。 $c_1$ と $c_2$ の位置関係で場合分けをする。

初めに、 $c_1$ が $c_2$ を包含する、または、$c_2$ が $c_1$ を包含する場合の共通部分はそれぞれ $c_1$ と $c_2$ に等しい。2 つの円の包含関係の判定は `Circle::contain` 関数を使う（c.f. [高校数学の美しい物語「https://manabitimes.jp/math/745」](https://manabitimes.jp/math/745)）。

他方を包含しない場合は 2 円の交点を通る直線と $c_1$ と $c_2$ の中心の位置関係で 2 通りの場合分けをする。どの場合も円板の **扇形（circular sector）** と **弓形（circular segment）** の面積を用いるので先にそれらの求め方を説明する。  
円板の扇形の面積は `Circle::area_circular_sector(Point2 p1, Point p2)` 関数で求める。円周上の 2 点 $p_1, p_2$ のなす角を $\theta$ （`arg` 関数）としたとき面積は、 $\frac{1}{2} r^2 \theta$ となる。  
次に、弓形の面積は `Circle::area_circular_segment(Point2 p1, Point2)` 関数で求める。弓形の面積は、扇形の面積から、中心と $p_1, p_2$ を頂点とする三角形の面積を引いた値に等しい。  

それでは、 $c_1$ と $c_2$ の共通部分の面積を求める。 $c_1$ と $c_2$ の交点を $p_1, p_2$ とする。
$c_1$ と $c_2$ の中心が直線 $p_1 p_2$ が分割する異なる領域にある場合は、共通部分の面積は 2 円の弓形の面積の和に等しい。

それ以外の、 $c_1$ と $c_2$ の中心が直線 $p_1 p_2$ が分割する同じ領域にある場合を考える。ここで、一般性を失わずに $c_2$ の中心が $c_1$ の中心よりも直線 $p_1 p_2$ に近いとする。このときの共通部分の面積は、 $c_2$ の面積に $c_1$ の弓型の面積を足して $c_2$ の弓型の面積を引いた値に等しい。

実装上の注意として、`cross_point` で求めた 2 円の交点は定義からそれぞれの円周上の点であることを期待されるが、数値誤差のため中心からの距離が半径に等しくない場合がある。

----------------------------------------------------------------------------------------------------------------------------------------------

<details>
<summary>src/CGL_7_I.cc を表示</summary>

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <cassert>
#include <vector>

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
inline Real to_radian(const Real degree) { return degree * (PI / 180.0); }
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

// ベクトル p1 と p2 の外積の絶対値 |p1 x p2| ： |p1 x p2| = |p1| |p2| sin(theta)
// 原点, p1, p2 を頂点とする符号付き平行四辺形の面積（ p1 から p2 へ反時計回りで符号が正）
inline Real abs_cross(const Point2 &p1, const Point2 &p2) { return p1.x * p2.y - p1.y * p2.x; }

// ベクトル p1 から p2 への角度を返す（単位はラジアン）
// p1 と p2 のなす角度で小さい方で p1 から p2 へ反時計回りなら符号は正，時計回りなら符号は負
inline Real arg(const Point2 &p1, const Point2 &p2) { return atan2(abs_cross(p1, p2), dot(p1, p2)); }

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
 * Circle in two dimensional
 */
class Circle : public Point2 {
public:
    Real r;
    explicit Circle() {}
    Circle(const Point2 &p, Real r = 0.0) : Point2(p), r(r) {}
    Circle& operator=(const Point2 &p) {
        this->x = p.x; this->y = p.y;
        return *this; 
    }

    Real area() const { return PI * r * r; }
    Real area_circular_sector(const Point2 &p1, const Point2 &p2, const bool strict_check = false) const; // 扇形の面積
    Real area_circular_segment(const Point2 &p1, const Point2 &p2, const bool strict_check = false) const; // 弓形の面積

    // Is p contained or on segment or otherwise? : O(n)
    // CONTAIN contain(const Point2 &p) const;
    bool contain(const Circle &rhs) const;
};

// Input of a circle
std::istream& operator>>(std::istream &is, Circle &c) { return is >> c.x >> c.y >> c.r; }
std::ostream& operator<<(std::ostream &os, const Circle &c) { return os << c.x << ' ' << c.y << ' ' << c.r; }


/**
 * Intersection testing
 */
Point2 projection(const Line &l, const Point2 &p) {
    Point2 dir = l[1] - l[0];
    Real t = dot(p - l[0], dir) / dir.abs2();
    return l[0] + dir * t;
}

// 2円が内接・2点で交わる・外接するなら true、含まれる・離れいているなら false
inline bool is_intersect(const Circle &c1, const Circle &c2) {
    return sign(c1.r + c2.r - (c1 - c2).abs()) >= 0 &&
        sign((c1 - c2).abs() - std::abs(c1.r - c2.r) >= 0);
}


/**
 * Distance and Intersection Point2
 */
inline Real distance(const Point2 &p1, const Point2 &p2) {
    return (p1 - p2).abs();
}

inline Real distance(const Line &l, const Point2 &p) {
    return (p - projection(l, p)).abs();
}

std::vector<Point2> cross_point(const Circle &c1, const Circle &c2) {
    if (!is_intersect(c1, c2))
        return std::vector<Point2>();
    Real d = distance(c1, c2);

    // Herbie による提案: Real r1_cos = (d * d + c1.r * c1.r - c2.r * c2.r) / (2.0 * d);
    Real r1_cos = 0.5 * (d + ((c1.r + c2.r) / d) * (c1.r - c2.r));
    Real h = std::sqrt(c1.r * c1.r - r1_cos * r1_cos);
    Point2 base = c1 + (c2 - c1) * r1_cos / d;
    Point2 dir = (c2 - c1).rotate90() * h / d;
    if (dir == Point2(0, 0))
        return {base};
    return {base + dir, base - dir};
}

Real Circle::area_circular_sector(const Point2 &p1, const Point2 &p2, const bool strict_check) const {
    // p1 または p2 が円周上の点ではない場合（数値誤差のため corss_point 関数で求めた点が円周上の点ではない場合がある）
    if (strict_check) {
        if (neq(r, distance(*this, p1)) || neq(r, distance(*this, p2))) return 0.0;
    }
    if (p1 == p2) return 0.0;
    return 0.5 * r * r * std::abs(::arg(p1 - *this, p2 - *this));
}

Real Circle::area_circular_segment(const Point2 &p1, const Point2 &p2, const bool strict_check) const {
        // p1 または p2 が円周上の点ではない場合（数値誤差のため corss_point 関数で求めた点が円周上の点ではない場合がある）
    if (strict_check) {
        if (neq(r, distance(*this, p1)) || neq(r, distance(*this, p2))) return 0.0;
    }
    Real area = this->area_circular_sector(p1, p2);
    if (eq(area, 0.0)) return 0.0;
    return area - 0.5 * std::abs(abs_cross(p1 - *this, p2 - *this));
}

// this が rhs を含む・内接するなら true、それ以外（2点で交わる・外接する・離れている）なら false を返す
bool Circle::contain(const Circle &rhs) const {
    return leq(rhs.r, this->r) && leq((*this - rhs).abs(), std::abs(this->r - rhs.r));
}


// --------------------8<------- start of main part of library -------8<--------------------
// 2円の共通部分の面積を求める
Real area_intersection(const Circle &c1, const Circle &c2) {
    if (c1.contain(c2)) return c2.area(); // c2 が c1 に含まれる場合
    if (c2.contain(c1)) return c1.area(); // c1 が c2 に含まれる場合
    if (!is_intersect(c1, c2)) return 0.0; // 離れている場合

    const auto ps = cross_point(c1, c2);
    if (ps.size() != 2) return 0.0; // 内接または外接する場合

    // 2点で交わる場合
    const Real dist_c1_ps = distance(Line(ps[0], ps[1]), c1);
    const Real dist_c2_ps = distance(Line(ps[0], ps[1]), c2);
    const Real dist_c1_c2 = distance(c1, c2);

    if (eq(dist_c1_ps + dist_c2_ps, dist_c1_c2) && neq(0.0, dist_c1_ps) && neq(0.0, dist_c2_ps)) {
        return c1.area_circular_segment(ps[0], ps[1]) + c2.area_circular_segment(ps[0], ps[1]);
    }
    else if (leq(dist_c1_ps, dist_c2_ps)) {
        return c1.area() - (c1.area_circular_segment(ps[0], ps[1]) - c2.area_circular_segment(ps[0], ps[1]));
    }
    else {
        return c2.area() - (c2.area_circular_segment(ps[0], ps[1]) - c1.area_circular_segment(ps[0], ps[1]));
    }
}
// --------------------8<------- end of main part of library -------8<--------------------


int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    Circle c1, c2;
    std::cin >> c1 >> c2;
    std::cout << area_intersection(c1, c2) << std::endl;

    return 0;
}
```

</details>

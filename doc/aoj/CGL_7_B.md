[🏠 Home](../../README.md)

# 内接円（Incircle of a Triangle）: [AOJ CGL_7_B](https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/all/CGL_7_B)
三角形 $\triangle ABC$ 内接円を求めよ。

三角形の内接円（ *incircle* ）とは、三角形の内部にあり三角形の 3 辺に接する円のことである。内接円の中心を内心（ *incenter* ）と呼ぶ。

# 解法
三角形の内接円は存在して、その半径 $r$ はヘロンの公式から導出され、

$ r = \frac{2S}{a + b + c}$

である。ただし、 $S$ は $\triangle ABC$ の面積で、$a, b, c$ はそれぞれ $\angle A, \angle B, \angle C$ の対辺の長さとする。すなわち、 $a + b + c$ は $\triangle ABC$ の周長に等しい。

内心は 3 頂点の重み付き平均

$\frac{a}{a + b + c} (x_a, y_a) + \frac{b}{a + b + c} (x_b, y_b) + \frac{c}{a + b + c} (x_c, y_c)$

に等しい。ここで、 $(x_a, y_a), (x_b, y_b), (x_c, y_c)$ はそれぞれ $\angle A, \angle B, \angle C$ の座標とする。

内接円は `incircle_triangle` 関数で求める。  
第二引数の `bool strict_definition` は十分条件をチェックしている。すなわち、上の式で求めた円が三角形の 3 辺に接するのかを確かめている。  
実験的に十分条件を満たさないインスタンスが存在するかは確認できていないが、数値誤差の影響を受けずに正しい内接円を返すかは確信できなかったので追加した。


## 参考文献
- [Wikipedia: 三角形の内接円と傍接円](https://ja.wikipedia.org/wiki/%E4%B8%89%E8%A7%92%E5%BD%A2%E3%81%AE%E5%86%85%E6%8E%A5%E5%86%86%E3%81%A8%E5%82%8D%E6%8E%A5%E5%86%86)



----------------------------------------------------------------------------------------------------------------------------------------------

<details>
<summary>src/CGL_7_B.cc を表示</summary>

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <array>
#include <cassert>
#include <vector>
#include <optional>

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

// dir と逆の向きを返す
CCW inv(const CCW dir) {
    switch (dir) {
        case (CCW::COUNTER_CLOCKWISE): return CCW::CLOCKWISE;
        case (CCW::CLOCKWISE): return CCW::COUNTER_CLOCKWISE;
        case (CCW::ONLINE_FRONT): return CCW::ONLINE_BACK;
        case (CCW::ONLINE_BACK): return CCW::ONLINE_FRONT;
        default: return dir;
    }
}

// Counter-Clockwise predicaste (a, b, c)
CCW ccw(const Point2 &a, Point2 b, Point2 c) {
    b -= a;  c -= a;
    if (sign(abs_cross(b, c)) == 1) return CCW::COUNTER_CLOCKWISE;
    if (sign(abs_cross(b, c)) == -1) return CCW::CLOCKWISE;
    if (sign(dot(b, c)) == -1)       return CCW::ONLINE_BACK;
    if (sign(b.abs2() - c.abs2()) == -1)   return CCW::ONLINE_FRONT;
    return CCW::ON_SEGMENT;
}

// ccw 関数の戻り値を CCW の規定型（int 型）で返す
auto ccw_t(const Point2 &a, Point2 b, Point2 c) {
    return static_cast<std::underlying_type<CCW>::type>(ccw(a, std::move(b), std::move(c)));
}

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

CCW ccw(const Line &l, const Point2 &p) { return ccw(l[0], l[1], p); }
auto ccw_t(const Line &l, const Point2 &p) { return ccw_t(l[0], l[1], p); }


/**
 * Segment in two dimensional
 */
class Segment : public Line {
public:
    explicit Segment() {}
    Segment(const Point2 &p1, const Point2 &p2) : Line(p1, p2) {}
};


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

inline bool is_intersect(const Circle &c, const Point2 &p) { // p is in interior or boundary
    return leq((c - p).abs(), c.r);
}

inline bool is_intersect(const Circle &c, const Segment &s) {
    return is_intersect(c, s[0]) || is_intersect(c, s[1]) ||
        (is_intersect(c, projection(s, c))
         && ccw(s[0], projection(s, c), s[1]) == CCW::ONLINE_FRONT);
}

/**
 * Distance and Intersection Point2
 */
inline Real distance(const Point2 &p1, const Point2 &p2) {
    return (p1 - p2).abs();
}

// 円 c と線分 s の交点を求める（s の端点も含む可能性がある）
std::vector<Point2> cross_point(const Circle &c, const Segment &s) {
    if (!is_intersect(c, s)) return std::vector<Point2>();

    const Point2 mid = projection(s, c), e = (s[1] - s[0]) / (s[1] - s[0]).abs();
    if (eq(c.r, (mid - c).abs())) return { mid };

    const Real len = sqrt(c.r * c.r - (mid - c).abs2());
    const Point2 p1 = mid + e * len, p2 = mid - e * len;
    const CCW ccw1 = ccw(s[0], p1, s[1]); 

    if (p1 == p2 && ccw1 == CCW::ONLINE_FRONT) return {p1};

    const CCW ccw2 = ccw(s[0], p2, s[1]);
    std::vector<Point2> ps;
    if (ccw1 == CCW::ONLINE_FRONT || p1 == s[1]) ps.emplace_back(p1);
    if (ccw2 == CCW::ONLINE_FRONT || p2 == s[1]) ps.emplace_back(p2);

    if (ps.size() == 2 && ccw(s[0], ps.back(), ps.front()) == CCW::ONLINE_FRONT)
        std::swap(ps.front(), ps.back());
    return ps;
}


/**
 * Polygon: The order of the points in the polygon is counter-clockwise.
 */
class Polygon {
public:
    explicit Polygon() {}
    explicit Polygon(int size) : points(size){}
    explicit Polygon(std::initializer_list<Point2> p) : points(p.begin(), p.end()) {}
    // explicit Polygon(std::vector<Point2> p) : points(std::move(p)) {}
    explicit Polygon(std::vector<Point2> &&p) : points(p) {}

    Real area() const; // area of polygon : O(n)

    std::vector<Point2> points;
};

// Output of a polygon
std::ostream& operator<<(std::ostream &os, const Polygon &poly) {
    for (auto p : poly.points) os << p << ", ";
    return os;
}

Real Polygon::area() const {
    const int n = points.size();
    assert(1 < n);

    Real area = abs_cross(points[n - 1], points[0]);
    for (int i = 0; i < n - 1; ++i)
        area += abs_cross(points[i], points[i + 1]);
    return 0.5 * area;
}


// --------------------8<------- start of main part of library -------8<--------------------

// Get incircle of a triangle
std::optional<Circle> incircle_triangle(const Polygon &triangle, const bool strict_definition = false) {
    if (triangle.points.size() != 3) return std::nullopt;

    const Real area = std::abs(triangle.area());
    const Real len_01 = distance(triangle.points[0], triangle.points[1]);
    const Real len_02 = distance(triangle.points[0], triangle.points[2]);
    const Real len_12 = distance(triangle.points[1], triangle.points[2]);

    Circle incircle;
    const Real perimeter = len_01 + len_02 + len_12; // perimeter of the triangle
    incircle.r = 2.0 * area / perimeter;
    incircle = (len_12 * triangle.points[0] + len_02 * triangle.points[1] + len_01 * triangle.points[2]) / perimeter;

    // Check that the incircle must touch the three sides of the triangle
    if (strict_definition) {
        auto p01 = cross_point(incircle, Segment(triangle.points[0], triangle.points[1]));
        auto p02 = cross_point(incircle, Segment(triangle.points[0], triangle.points[2]));
        auto p12 = cross_point(incircle, Segment(triangle.points[1], triangle.points[2]));

        if (p01.size() != 1 || p02.size() != 1 || p12.size() != 1) return std::nullopt;

        if (neq(incircle.r, (p01.front() - incircle).abs())) return std::nullopt;
        if (neq(incircle.r, (p02.front() - incircle).abs())) return std::nullopt;
        if (neq(incircle.r, (p12.front() - incircle).abs())) return std::nullopt;
    }

    return incircle;
}
// ---------------------8<------- end of main parto of library -------8<---------------------

int main() {
    std::cout << std::fixed << std::setprecision(10);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    Polygon triangle(3);
    for (auto &p : triangle.points) {
        std::cin >> p;
    }

    auto incircle = incircle_triangle(triangle);
    if (incircle) {
        std::cout << incircle.value() << std::endl;
    }
    else {
        exit(EXIT_FAILURE);
    }

    return 0;
}
```

</details>

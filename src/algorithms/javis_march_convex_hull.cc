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

// Counter-Clockwise predicate (a, b, c)
enum class CCW : int {
    COUNTER_CLOCKWISE = 1,     // counter clockwise
    CLOCKWISE         = -1,    // clockwise
    ONLINE_FRONT      = 2,     // a--b--c on line or (a == b and b != c)
    ONLINE_BACK       = -2,    // c--a--b on line
    ON_SEGMENT        = 0,     // a--c--b on line or (a != b and b == c) or (a == b == c)
    OTHER             = -3,
};

CCW ccw(const Point2 &a, Point2 b, Point2 c) {
    b -= a;  c -= a;
    if (sign(abs_cross(b, c)) == 1) return CCW::COUNTER_CLOCKWISE;
    if (sign(abs_cross(b, c)) == -1) return CCW::CLOCKWISE;
    if (sign(dot(b, c)) == -1)       return CCW::ONLINE_BACK;
    if (sign(b.abs2() - c.abs2()) == -1)   return CCW::ONLINE_FRONT;
    return CCW::ON_SEGMENT;
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

     std::vector<Point2> points;
};


// --------------------8<------- start of main part of library -------8<--------------------

// Jarvis's March: O(n h) (n = #ps, h = #Polygon)
Polygon convex_hull(std::vector<Point2> ps) {
    if (ps.size() < 3) return Polygon(std::move(ps));

    const int n = ps.size();
    // find the smallest y-coordinate and set on ps[0]
    for (int i = 1; i < n; ++i) {
        if (ps[i].y < ps.front().y || ((ps[i].y == ps.front().y) && (ps[i].x < ps.front().x))) {
            std::swap(ps.front(), ps[i]);
        }
    }

    int next = 1;
    for ( ; next < n; ++next) {
        for (int i = next + 1; i < n; ++i) {
            const auto check = ccw(ps[next - 1], ps[next], ps[i]);
            if (CCW::CLOCKWISE == check || CCW::ON_SEGMENT == check) {
                std::swap(ps[next], ps[i]);
            }
        }
        if (CCW::CLOCKWISE == ccw(ps[next - 1], ps[next], ps.front())) {
            --next;
            break;
        }
    }

    if (next + 1 < n) ps.resize(next + 1);

    return Polygon(std::move(ps));
}

// --------------------8<------- end of main part of library   -------8<--------------------

int main() {
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    // [AOJ: CGL_4_A](https://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_4_A&lang=ja)
    int n;
    std::cin >> n;

    std::vector<Point2> ps(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> ps[i];
    }

    Polygon cv = convex_hull(ps);
    std::cout << cv.points.size() << '\n';
    for (const auto &p: cv.points) {
        std::cout << p << '\n';
    }

    return 0;
}

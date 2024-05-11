#include <iostream>
#include <vector>
#include <algorithm>

struct Point {
    int x = 0, y = 0;
    Point() {}
    Point(int _x, int _y) : x(_x), y(_y) {}
    Point operator-(const Point &rhs) { return Point(x - rhs.x, y - rhs.y); }
    Point operator-=(const Point &rhs) { return *this = *this - rhs; }
};

std::istream& operator>>(std::istream &is, Point &p) { return is >> p.x >> p.y; }
int distance(const Point &p1, const Point &p2) { return std::abs(p1.x - p2.x) + std::abs(p1.y - p2.y); }
int abs_cross(const Point &p1, const Point &p2) { return p1.x * p2.y - p1.y * p2.x; }

bool IsCounterClockwise(const Point &a, Point b, Point c) {
    b -= a; c -= a;
    return abs_cross(b, c) >= 0;
}

std::vector<int> GetCodedPolygonalLine() {
    int m;
    std::cin >> m;
    std::vector<int> poly;
    poly.reserve(m - 1 + m - 2);
    Point prev, mid, cur;

    std::cin >> prev >> mid;
    poly.push_back(distance(prev, mid));

    for (int i = 2; i < m; ++i) {
        std::cin >> cur;
        poly.push_back(IsCounterClockwise(prev, mid, cur));
        poly.push_back(distance(mid, cur));
        std::swap(prev, mid);
        std::swap(mid, cur);
    }

    return poly;
}

bool SamePolygonalLine(const std::vector<int> &pl1, const std::vector<int> &pl2) {
    if (pl1.size() != pl2.size()) return false;

    const bool forward = std::equal(pl1.begin(), pl1.end(), pl2.begin(), pl2.end());
    if (forward) return true;

    int idx = 0;
    return std::equal(pl1.begin(), pl1.end(), pl2.rbegin(), pl2.rend(),
                    [&idx](int p1, int p2) {
                        if (idx++ % 2 == 0) { return p1 == p2; }
                        else { return p1 != p2; }
                    });
}

void PolygonalLineSearch(const int n) {
    std::vector<int> original = GetCodedPolygonalLine();

    for (int i = 0; i < n; ++i) {
        auto poly = GetCodedPolygonalLine();
        if (SamePolygonalLine(original, poly)) {
            std::cout << i + 1 << "\n";
        }
    }
}

int main() {
    int n;
    while (std::cin >> n, n) {
        PolygonalLineSearch(n);
        std::cout << "+++++\n";
    }

    return 0;
}
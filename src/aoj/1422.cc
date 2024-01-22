#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>


// --------------------8<------- start of main part of solution -------8<--------------------
const double PI = acos(static_cast<double>(-1.0));

class SolveLoopOfChocolate {
public:
    SolveLoopOfChocolate(const int _n, const int _r) : n(_n), r(_r) { spheres.reserve(n); }

    struct Point3 { int x, y, z; };
    const int n; // #sphere (4 <= n <= 100)
    const int r; // radius of spheres; all of them are same (2 <= r <= 100)
    std::vector<Point3> spheres; // center of spheres

    void AddSphere(int x, int y, int z) { spheres.emplace_back(Point3{x, y, z}); }

    double GetVolumeUnionSpheres() const {
        double volume = n * VolumeSphere();

        for (int i = 0; i + 1 < n; ++i) {
            volume -= VolumeIntersectSphere(Distance(spheres[i], spheres[i + 1]));
        }
        volume -= VolumeIntersectSphere(Distance(spheres[0], spheres[n - 1]));

        return volume;
    };

private:
    // distance of two sphere's center
    double Distance(const Point3 &s1, const Point3 &s2) const { 
        return std::hypot(s1.x - s2.x, s1.y - s2.y, s1.z - s2.z); 
    }

    // volume of a sphere of radius r
    double VolumeSphere() const { return PI * r * r * r * 4.0 / 3.0; }

    // volume of the intersection of two spheres when distance between their center is d
    double VolumeIntersectSphere(const double d) const {
        return PI * 2.0 / 3.0 * std::pow(r - 0.5 * d, 2) * (2.0 * r + 0.5 * d);
    }
};
// --------------------8<------- end of main part of solution   -------8<--------------------


int main() {
    std::cout << std::fixed << std::setprecision(8);
    std::cin.tie(0); std::ios::sync_with_stdio(false);

    int n, r;
    std::cin >> n >> r;
    auto solver = SolveLoopOfChocolate(n, r);
    for (int i = 0; i < n; ++i) {
        int x, y, z;
        std::cin >> x >> y >> z;
        solver.AddSphere(x, y, z);
    }

    std::cout << solver.GetVolumeUnionSpheres() << std::endl;

    return 0;
}
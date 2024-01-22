#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>

using ll = long long;

template<class T>
std::vector<std::pair<T, T>> PrimeFactorization(T n) {
    std::vector<std::pair<T, T>> pf;
    T m = n;
    for (T i = 2; i * i <= n; ++i) {
        if (m % i != 0) continue;
        T cnt = 0;
        while (m % i == 0) { ++cnt; m /= i; }
        pf.emplace_back(std::make_pair(i, cnt));
    }
    if (1 < m) pf.emplace_back(std::make_pair(m, 1));

    return pf;
}

int main() {
  ll p, q;
  std::cin >> p >> q;

  q /= std::gcd(p, q);
  const auto pf = PrimeFactorization(q);

  ll ans = std::accumulate(pf.begin(), pf.end(), (ll)1, [](ll acc, const auto &pf_i) {
    return acc * pf_i.first;
  });

  std::cout << ans << std::endl;

  return 0;
}
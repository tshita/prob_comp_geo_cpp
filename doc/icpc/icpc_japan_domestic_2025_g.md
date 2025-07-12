[🏠 Home](../../README.md)  |  [🔗 Problem URL](https://storage.googleapis.com/files.icpc.jp/domestic2025/problems/problems_ja.pdf#page=11)

# 問題：面の数
3次元ユークリッド空間中に 2 つの平面 $H_1: z = 1$ と $H_2: z = 2$ がある。  
実数列 $d = (d_1, \ldots, d_n)$ と $d' = (d'_1, \ldots, d'_m)$ が与えられる。平面 $H_1, H_2$ に頂点の内角が原点から見て反時計周り順にそれぞれ $d, d'$ となるように凸多角形を構成して、それら2 つの凸多角形を含む $n + m$ 個の頂点からなる凸多面体の面の数としてあり得るものを全列挙せよ。

**制約 :** $3 \le n, m \le 50$、$10^{-9} \le d_i, d_i' < 180$

# 解法
$H_1$ と $H_2$ 上にある凸多角形をそれぞれ $P_1, P_2$ とし、面の数を求める凸多角形を $P = \rm{conv}(P_1 \cup P_2)$ とする。  
$P$ は明らかに有界なので凸多面体である。凸多面体の Minkowski-Weyl の定理から $P$ の H 表現がある。この問題は H 表現で考察するとよい。
$P_1$ と $P_2$ の任意の辺に対してその辺を含む $P$ のファセットが存在し、構成方法からファセットは 3 角形か 4 角形であることが分かる。特に、4 角形のとき $P_1$ と $P_2$ 上にある辺の法線ベクトルの向きは互いに等しい。  

平面上にある凸多角形で内角の和が指定されているとき、辺を構成する 1 つの超平面をひとつ固定すれば残りの辺の法線ベクトルの向きは一意に定まる。すべてのファセットが 3 角形となる凸多角形は必ず存在する。なぜならば、$P_1$ と $P_2$ のどの辺の法線ベクトルの向きも異なるものを構成できるからである。そのときの面の数は $n + m + 2$ である。一般に 4 角形の面が $k$ 個ある凸多角形の面の数は $n + m + 2 - k$ である。よって、4 角形の面が 1 つ以上ある凸多角形を考えればよく、それは $H_1$ と $H_2$ に垂直な超平面をひとつ固定して、$P_1$ と $P_2$ の辺の組合せすべてに対して、各辺の組がその超平面によって構成されるとして凸多角形はを構成すればよい。

詳しくは [忘れても大丈夫](https://kyopro.hateblo.jp/entry/2025/07/07/125501) に書いた。

---------------------------------------------------------------------------------------------

<details>
<summary>src/icpc/icpc_japan_domestic_2025_g.cc を表示</summary>

```cpp
#include <iostream>
#include <vector>
#include <algorithm>

// 加算時の丸め誤差を考慮
bool is_same(double lhs, double rhs) {
    return std::abs(lhs - rhs) < 1e-10;
}

int main() {
    int n, m;

    while (std::cin >> n, n) {
        // Input H1
        std::vector<double> in_ang_h1(n);
        for (auto &in_ang: in_ang_h1) {
            std::cin >> in_ang;
        }

        // Input H2
        std::cin >> m;
        std::vector<double> in_ang_h2(m);        
        for (auto &in_ang: in_ang_h2) {
            std::cin >> in_ang;
        }

        // 面の数の上限（面がすべて三角形の場合）
        const int MAX_NUM_FACE = n + m + 2;
        // exist_num_face[i] = true の場合は面の数が i となる凸多面体が存在する
        std::vector<bool> exist_num_face(MAX_NUM_FACE + 1, false);

        // 面の数が MAX_NUM_FACE を達成する凸多面体が存在する
        exist_num_face[MAX_NUM_FACE] = true;

        // 多角形 H1 と H2 の各辺の法線ベクトル
        std::vector<double> norm1(n + 1), norm2(m + 1);

        // 多角形 H1 と H2 のそれぞれ i 番目と j 番目の辺を対応させる
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                // H1 の各辺の法線ベクトルを設定（i 番目の辺の法線ベクトルは 0°）
                for (int idx = 0; idx < n; ++idx) {
                    norm1[idx + 1] = norm1[idx] + 180.0 - in_ang_h1[(i + idx) % n];
                }
                // H2 の各辺の法線ベクトルを設定（j 番目の辺の法線ベクトルは 0°）
                for (int idx = 0; idx < m; ++idx) {
                    norm2[idx + 1] = norm2[idx] + 180.0 - in_ang_h2[(j + idx) % m];
                }

                // H1 と H2 の辺の法線ベクトルが等しい組数を求める
                int num_match_pair = 0;
                int idx2 = 0;
                for (int idx1 = 0; idx1 < n; ++idx1) {
                    while (idx2 < m && !is_same(norm1[idx1], norm2[idx2]) && norm2[idx2] < norm1[idx1]) {
                        ++idx2;
                    }

                    if (idx2 < m && is_same(norm1[idx1], norm2[idx2])) {
                        ++num_match_pair;
                    }
                }

                // 現在の凸多面体の面の数を設定
                exist_num_face[n + m + 2 - num_match_pair] = true;
            }
        }

        int num = std::count(exist_num_face.begin(), exist_num_face.end(), true); // 出力フォーマットのため
        for (int i = 1; i <= MAX_NUM_FACE; ++i) {
            if (exist_num_face[i]) {
                std::cout << i;
                if (--num) std::cout << " ";
            }
        }
        std::cout << std::endl;
    }

    return 0;
}
```

</details>

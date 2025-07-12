#include <iostream>
#include <vector>
#include <algorithm>

// 加算時の丸目誤差を考慮
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
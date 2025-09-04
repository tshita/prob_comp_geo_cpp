> Back to [Table of Contents](../index.md)

- **線分交差問題に対する平面走査法（plane sweep algorithm）**  
    線分交差問題とは平面上にある線分の集合が与えられたときに交点をすべて答える問題。平面走査法は $O(n \log n + I \log n)$ 時間で解くアルゴリズム。ただし、$I$ は線分の間の交点数で、このように実行時間が出力サイズに依存するアルゴリズムを **出力サイズに敏感なアルゴリズム（output-sensitive algorithm）** という。詳細は **[2, pp. 21-33]** を参照。
- **マーチング・キューブ法（Marching cubes）**  
    三次元ボリュームデータから等値面のポリゴンメッシュを抽出する、 Lorensen と Cline によるアルゴリズム。  
    閾値によって各ボクセルの値を 0 または 1 に二値化して、その値を頂点に持つ立方体に対して、あらかじめ定義された 15 種類のパターンに基づき、0 と 1 の頂点を分割するように辺上にメッシュの頂点を配置する。このとき、辺に配置される頂点の位置は、その端点のボクセル値を用いた線形補間によって決定される。
    - William E. Lorensen, Harvey E. Cline, ``Marching cubes: A high resolution 3D surface construction algorithm''. ACM SIGGRAPH Comput. Graph., 21 (1987) pp. 163-169. [DOI:10.1145/37402.37422](https://dl.acm.org/doi/pdf/10.1145/37402.37422) (Open Access)
    - [tatsy. Programming for Beginners Programming for everyone!](https://tatsy.github.io/programming-for-beginners/cpp/march-cubes/)：とても分かりやすい解説。
    - Nintendo Switch 2 のゲームソフト「[ドンキーコング バナンザ](https://www.nintendo.com/jp/games/switch2/aaaca/index.html)」では、破壊可能なオブジェクトの管理・描画において、マーチング・キューブ法に類似した手法が用いられている（特許出願 [特願2024-211588](https://www.j-platpat.inpit.go.jp/c1801/PU/JP-2025-113168/11/ja) を参考にした私見）。


# References
[1] T.コルメン，C.ライザーソン，R.リベスト，C.シュタイン（浅野哲夫・岩野和生・梅尾博司・山下雅史・和田幸一 訳）：『アルゴリズムイントロダクション 第3版 総合版：世界標準MIT教科書』．近代科学社，東京，2013．  
[2] M.ドバーグ，O.チョン，M.ファンクリベルド，M.オーバマーズ（浅野哲夫 訳）：『コンピュータ・ジオメトリ 計算幾何学：アルゴリズムと応用 第3版』. 近代科学社, 東京, 2010.
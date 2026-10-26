# combkit：組合數學 × ManimGL 動畫工具包

從 catalan930 專案（卡塔蘭數、Chu 的廣義卡塔蘭數、MIT 18.212、Fulton《Young Tableaux》的演示）整理出來的程式，
可以直接拿到別的計畫用；1.6.0 起也收了 spglass（SK 自旋玻璃）專案的隨機矩陣組合。版本 1.10.1（2026-10-16）。
整台電腦的 ManimGL 環境、依賴與檔案結構，見同資料夾的 `環境與交接_ENVIRONMENT.md`（副本；正本在 `C:\mypython\manim`）。

**路徑的寫法：** 本文件裡的相對路徑都相對於本文件所在的資料夾 `C:\mypython\manim\projects\combkit`
（例如 `examples/ex01_dyck_paths.py` 就是 `C:\mypython\manim\projects\combkit\examples\ex01_dyck_paths.py`）。

---

## 1. 核心想法：先算，再畫

所有演示都分成三層，數學和畫面分開：

```text
 數學層（純 Python，不需要 manim）          畫圖層                         成品
 combkit.counting ─┐                    ┌─ combkit.draw.*（ManimGL）──→ 動畫場景
 combkit.tableaux ─┼─→ 資料（樹、字、表）─┤
 combkit.raney    ─┤
 combkit.rmt_moments ┤
                   │   每個結果都用暴力   └─ combkit.tikz（TikZ 字串）──→ LaTeX 講義的圖
                   │   列舉交叉驗證
```

例如 Dyck 路徑：

1. `counting.binary_trees(3)` 列出 5 棵二元樹；`counting.tree_to_dyck(t)` 把每棵樹變成 `"UUDUDD"` 這樣的字；
   `counting.dyck_words(3)` 用另一種方法（直接暴力列舉）再算一次，兩者必須相同。
2. `draw.catalan_draw.dyck_item(t, colored=True)` 把同一棵樹畫成 manim 物件，依「根／左子樹／右子樹」上色；
   `tikz.dyck_prims(...)` 則畫成講義用的 TikZ。
3. 場景只負責排版與動畫，畫面上出現的每個數字都 `assert` 和數學層一致。

好處：數學錯誤在渲染之前就會被抓到；同一份資料可以同時產生動畫與講義的圖；
因為所有卡塔蘭物件都從同一棵樹產生，三角剖分、括號、Dyck 路徑、圓上配對……的「第 k 個」彼此真的一一對應。

---

## 2. 需求與使用方式

| 部分 | 需要 |
|---|---|
| `combkit.counting`、`combkit.tableaux`、`combkit.rep_data`、`combkit.raney`、`combkit.raney_examples`、`combkit.schubert`、`combkit.fulton_appendix`、`combkit.rmt_moments`、`combkit.tikz`、`combkit.tikz_trees` | 只要 Python 3.10 以上，不需要其他套件 |
| `combkit.draw.*` | 已安裝的 **ManimGL**（3Blue1Brown 版，不是 Manim Community；用 3b1b/manim commit fafa083 = v1.7.2 測過）、LaTeX（MiKTeX＋`cm-super`）、ffmpeg、中文字型「微軟正黑體」 |

兩種用法擇一：

- **A. 不安裝**：把整個 `combkit` 資料夾放在任何地方，在程式開頭加
  ```python
  import sys; sys.path.insert(0, r"C:\mypython\manim\projects\combkit")
  ```
  （`examples/` 裡的範例就是這樣寫的。）
- **B. 安裝成套件**（可編輯模式，之後改程式不用重裝）：
  ```bash
  python -m pip install -e C:/mypython/manim/projects/combkit
  ```
  `pyproject.toml` 刻意**沒有**列出 ManimGL，所以 pip 不會去換掉你固定版本的 ManimGL。

使用 `combkit.draw` 時，`manimgl` 必須在放有 `custom_config.yml` 的資料夾裡執行（`examples/custom_config.yml` 可以直接複製；它只能寫英文）。

---

## 3. 共用的資料表示法

| 物件 | 表示法 | 例子 |
|---|---|---|
| 二元樹（卡塔蘭的「樹」） | `None` = 空（葉子）；`(L, R)` = 一個內部節點 | `((None,(None,None)),None)` |
| Dyck 路徑 | `'U'`（上坡）、`'D'`（下坡）組成的字串 | `"UUDUDD"` |
| 格子路徑 | `'R'`（右）、`'U'`（上） | `"RURRUU"` |
| 圓上的點、配對 | 0 起算，從正上方順時針；配對是 `(a, b)` 的串列 | `[(0,5),(1,2),(3,4)]` |
| 凸多邊形的頂點 | 0 起算，底邊是 `0`–`(m−1)`，第 0 點在左下、順時針 | 三角剖分 = 對角線串列 |
| 分割 λ | 非遞增 tuple | `(4,3,1)` |
| Young 表 | 列 tuple 的 tuple；格子 `(i, j)` 從 0 起算，第 0 列在最上面 | `((1,2,4),(3,5,6))` |
| 斜表 λ/μ | 同上，第 i 列前面放 μ_i 個 `None` | `((None,None,2),(1,3))` |
| 多項式 | `{指數 tuple: 係數}` | Schur 多項式 |
| 機率 | 精確的 `Fraction` | 鉤上漫步的角落機率 |
| Raney 的字串／字列 | 自然數 tuple，第 x 項是第 x 個符號的元數（0 = 引數），前序記法；組成 $(n;a_1,a_2,\ldots)$ 寫成 `n` 與 tuple `a` | `(2,3,0,2,0,0,0,1,0)`、`L_formula(2,(0,1,1))` |
| 1/N 的多項式（隨機矩陣動差） | `{N 的冪次: 係數}` | GOE 的 E Tr J⁴ = 2N + 5 + 5/N 寫成 `{1: 2, 0: 5, -1: 5}` |
| 閉路徑 | 頂點序列（list），頂點按第一次到訪編號、從 0 起算 | `[0,1,0,2,3,2,0]` |
| 根／左／右 | 字串 `'root'`、`'left'`、`'right'`，畫圖時對應橘／青／粉紅（TikZ：`cTop`／`cLeft`／`cRight`） | `dyck_parts(t)` |

---

## 4. 模組說明

完整的函式清單（簽名＋說明第一行）見 `API.md`，由 `tools/make_api.py` 自動產生。

### 4.1 `combkit.counting`（純 Python）

直接執行 `python -m combkit.counting` 會跑完全部自我檢查（約 3 秒）。

| 節 | 內容 | 主要函式 |
|---|---|---|
| 1 | 卡塔蘭數、Chu 的廣義卡塔蘭數公式 | `catalan(n)`、`inner(n_vec, v_vec)`、`closed_form(n_vec, v_vec)` |
| 2 | 圓周上的不交叉分組（暴力列舉） | `blocks_cross`、`is_noncrossing`、`enumerate_noncrossing`、`brute_force` |
| 3 | Chu 的遞迴（固定一點的分解） | `recurrence(n_vec, v_vec)` |
| 4 | Raney 公式的證明：路徑、好路徑、循環引理、切兩段 | `step_weight`、`path_heights`、`is_good_path`、`raney_A`、`path_words`、`good_paths`、`rotate`、`good_rotations`、`cut_index` |
| 5 | 經典卡塔蘭物件，全部由同一棵樹產生 | `binary_trees(n)`、`tree_to_paren_tokens`、`tree_to_triangulation`、`triangulations`、`tree_to_dyck`、`dyck_parts`、`dyck_words`、`dyck_to_lattice`、`lattice_words`、`reflect_tail`、`dyck_to_brackets`、`dyck_to_ballot`、`bertrand_strict_count`、`dyck_to_matching`、`matching_parts`、`series_square_shift` |
| 6 | MIT 18.212 L1–4 補充 | 生成函數：`series_mul`、`sqrt_1m4x_coeffs`、`catalan_from_sqrt`；醉漢走路：`drunkard_fall_words`、`drunkard_partial_sums`；循環移位：`cycle_classes`；平面樹：`binary_to_plane`、`plane_trees`、`plane_tree_dfs`；對偶樹：`triangulation_dual`；弧圖：`perfect_matchings`、`matching_is_noncrossing`、`matching_is_nonnesting`；排序：`stack_sort`、`queue_bypass_sort`；模式：`contains_pattern`、`avoiders`；分割與 Young 表：`integer_partitions`、`hook_length`、`hook_product`、`hook_formula`、`syt`、`syt_to_dyck`；鉤上漫步：`hook_walk_probs`、`corner_prob_formula`、`hook_identity` |

### 4.2 `combkit.tableaux`（純 Python）

Fulton《Young Tableaux》Notation 與第 1 章。`python -m combkit.tableaux` 用暴力列舉驗證 Row Bumping Lemma、
乘積的結合律、rectification 與挑法無關、Rect(T∗U) = T·U、Schur 多項式的對稱性等（約 4 秒）。

- 分割：`partitions`、`conjugate`、`corners`、`skew_boxes`、`is_horizontal_strip`、`is_vertical_strip`
- 表：`is_tableau`、`is_standard`、`shape`、`transpose`、`content`、`ssyt(lam, m)`、`standard_tableaux`
- 列插入：`row_insert(T, x)` → (新表, 擠壓路徑, 新格子)、`bump_steps`（逐步動畫用）、`reverse_bump`、`tableau_word`、`insert_word_steps`、`product(T, U)`
- jeu de taquin：`inside_corners`、`outside_corners`、`forward_slide(S, corner)` → (新斜表, 洞的路徑)、`reverse_slide`、`jdt_steps`、`rectify`、`star(T, U)`
- Schur：`schur_poly(lam, m)`、`complete_h`、`elementary_e`、`kostka`
- 書上的例子：`FULTON_TABLEAU`、`FULTON_STANDARD`、`FULTON_INSERT_T`、`FULTON_PRODUCT_U`、`FULTON_SKEW`
- 第 2–4 章（第 7 節）：
  - 字：`word_runs`、`word_to_tableau`、`skew_from_word`、`column_word`
  - Knuth 等價：`knuth_moves(w)` → [(位置, 規則, 新字)]（K′、K″ 與反向）、`knuth_class`、`knuth_path(w1, w2)`、`insertion_knuth_chain(T, x)`（把一次列插入拆成基本變換）
  - 遞增子序列：`max_increasing(w, k, decreasing)` → (L(w,k), 一組子序列)、`greene_numbers`、`remove_largest`、`tableau_from_greene`（第 3 章唯一性證明當演算法）
  - R-S-K：`rsk_steps(top, bottom=None)`、`rsk`、`inverse_rsk(P, Q)`、`matrix_to_array`、`array_to_matrix`、`rsk_matrix`、`transpose_matrix`
  - 矩陣球：`ball_numbering(A)`、`derived_matrix(A)`、`matrix_ball(A)` → (P, Q, 各輪)
  - 其他：`pieri_shapes`、`poly_mul`、`poly_sum`、`lex_greater`、`kostka_matrix(n)`、`hook_content_count(λ, m)`、`count_involutions(n)`
  - 書上的例子：`CH2_WORD`、`CH2_SLIDE_*`、`CH3_WORD`、`CH3_TABLEAU`、`CH4_WORD`、`CH4_ARRAY`、`CH4_MATRIX`、`CH4_P`、`CH4_Q`、`CH4_HOOK_SHAPE`
- 第 5–6 章（第 8 節）：
  - LR 規則：`factorizations(V, λ, μ)`（T·U = V 的所有寫法）、`rectifying_skews(ν, λ, U0)`、`is_reverse_lattice`、`lattice_tableau(μ)`、`lr_tableaux(ν, λ, μ)`、`lr_coefficient(λ, μ, ν)`、`lr_product(λ, μ)`、`skew_decomposition(ν, λ)`、`reverse_numbering`、`remmel_whitney`
  - 對稱多項式：`monomial_sym`、`power_sum`、`h_lambda`、`e_lambda`、`p_lambda`、`z_lambda`、`schur_expand(f, m)`（寫成 Schur 函數的組合）、`poly_det`、`jacobi_trudi(λ, m, dual)`、`alternant`、`principal_specialization`、`character_table(n)`、`schur_inner`、`omega`
  - 書上的例子：`CH5_RL_GOOD`、`CH5_RL_BAD`、`CH5_NU`、`CH5_LAM`、`CH5_DECOMP`、`CH5_S21_SQUARED`、`CH5_ZERO`

### 4.2b `combkit.rep_data`（純 Python）

Fulton 第 7 章（對稱群的表示）的固定算例與小工具。`load(key)` 讀套件內 `data/fulton7_math.json`（已驗證的算例：主例 λ=(2,1)、(3,2)、內容 (2,2,1) 的 Young's rule、特徵標、兩輪 straightening、模 3…）；`act`、`compose`、`sign`、`cycles`、`cycle_type`（置換用一行記法）；`tabloid_key`、`polytabloid`（$v_T$）、`act_vec`、`specht21`、`rep_matrix21`（第 $j$ 欄 = 第 $j$ 個基底的像）；`gr_mul`、`symmetrizers`（$a_T,b_T,c_T$）。`python -m combkit.rep_data` 會用 `tableaux` 獨立重算並核對資料。

### 4.2c `combkit.raney`（純 Python）

G. N. Raney, *Functional composition patterns and power series reversion* (1960) 的計算；catalan930 的 `notes/raney1960_導讀.pdf`（33 頁）的數字與圖都出自這裡。記號照 Raney：字串 σ 是元數的 tuple，秩 $r(\sigma)=\sum(\sigma(x)-1)$，字＝前序寫出的平面樹，$n$ 字列＝$n$ 棵樹的森林。
- 字與字列：`rank`、`head_ranks`（秩路徑）、`is_word`、`is_list`（Rosenbloom 判別法；`is_list_recursive` 照定義）、`split_words`、`word_to_tree`、`tree_to_word`、`composition`
- 計數：`L_formula(n, a)`（定理 2.2：$mL=nM$）、`L_brute`、`lists_with`、`strings_with`、`multinomial`、`a0_of`、`length_of`、`thm23_rhs`、`thm24_rhs`、`split_head`
- 循環引理：`cyclic_variant`、`list_shifts`、`raney_alpha`、`raney_proof_trace`（照 Raney 的證明走一遍：α、σ₁、切成的字、n 種輪換）
- 和 Chu 講義（`counting` 第 4 節的路徑）互換：`to_chu_path`／`from_chu_path`（倒過來讀、符號 j 走 1−j）、`chu_heights`、`a_from_chu`；不交叉分法 ↔ 字：`partition_to_word`、`word_to_partition`；Gessel 的後序碼 `gessel_suffix_code`、`prefix_code`
- 冪級數（係數用 `Fraction`）：`ps_mul`、`ps_pow`、`ps_compose`、`solve_w`（解 $w=x\phi(w)$）、`lagrange_coeff(phi, n, t)`、`exp_series`；多變數的 $L^*$：`lstar_coeffs`、`mv_mul`、`mv_pow`
- 第 6–7 節：`trinomial_coeff`（Fuss–Catalan）、`kelly_P`／`kelly_P_brute`（Problem 4983）、`polygon_dissections_brute`、`a_table`、`b_from_a`、`c_from_b`、`gould_G`、`gould_phi`
- 第 9 節（整數序列的 Raney 引理，catalan930 的 `notes/ran_note015_改寫.pdf` 與動畫 R01–R06）：`heights`、`reading(seq, i)`（從第 i 點開始讀一圈）、`positive_count`、`is_good_reading`、`good_starts`、`last_minimum`、`positive_count_order`（加強版的排序法）、`periodic_heights`；出棧：`stack_sequences`、`stack_ops_to_steps`、`pop_order`；摸牌（洛谷 P6672）：`draws_left`、`deck_ok`、`card_game_brute`、`card_game_formula`；主例 `NOTE015_SEQ = (2,-1,-3,1,2)`

`python -m combkit.raney` 用暴力列舉驗證每個定理，並和 `counting` 交叉核對。

### 4.2d `combkit.rmt_moments`（純 Python）

高斯變數數「所有配對」，大隨機矩陣的特徵值數「不交叉配對」——spglass 專案（`C:\mypython\spglass`）先修講義練習解答的組合部分。
- 配對：`double_factorial_odd(k)`（$(2k-1)!!$）、`pairing_counts(k)` → (全部配對數, 不交叉配對數)
- 有限 $N$ 的精確動差（Wick 配對＋數獨立指標）：`goe_moment_poly(k)`、`gue_moment_poly(k)`（Harer–Zagier）、`normalized_moment(poly, N)`（精確 `Fraction`）。
  約定：GOE 非對角方差 $1/N$、對角 $2/N$；GUE $\mathbb{E}|H_{ab}|^2=1/N$。例：GOE $\frac1N\mathbb{E}\,\mathrm{Tr}\, J^4=2+\frac5N+\frac5{N^2}$、$\frac1N\mathbb{E}\,\mathrm{Tr}\, J^6=5+\frac{22}N+\frac{52}{N^2}+\frac{41}{N^3}$
- 存活的路徑：`tree_walks(k)`（$C_k$ 條樹形閉路徑）、`walk_to_dyck`、`walk_pairing`（同一條邊的兩次經過配成一對＝不交叉配對）

`python -m combkit.rmt_moments` 用逐指標暴力計算（直接對獨立高斯元素取動差）核對 GOE 多項式，並核對路徑 ↔ Dyck ↔ 配對的一一對應。

### 4.2d `combkit.raney_examples`（純 Python）

用循環引理「製造」例子：把符號隨機打亂、再隨機挑一個好的起點轉過去，得到的字列在同組成的字列中是**均勻分布**的（每個字列恰好被 $m$ 個（字串, 起點）讀到）。
- 抽樣：`random_string`、`random_list(n, a, rng, trace=False)`、`random_word`、`random_composition`、`random_sum_one`、`random_noncrossing(n_vec, v_vec, rng)`（均勻隨機的字 → 「固定點」的一一對應 → 均勻隨機的不交叉分法）、`multiset`
- 字 → 運算式：`word_to_expression(word)`（回傳前序寫法與一般寫法，例如 `('FTaGbcdfe', 'F(T(a,G(b,c),d),f(e))')`）
- 例題集：`examples_markdown(seed, count)`，八類題目（運算式與字、循環引理、數字列、不交叉分法、Lagrange、卷積、整數序列、摸牌），每題答案都由程式算出並 assert；命令列 `python -m combkit.raney_examples --md --seed 7 --count 3 --out 檔名`

`python -m combkit.raney_examples`（不加參數）是自我檢查：用計數檢查均勻性。

### 4.2e `combkit.schubert`（純 Python）

Fulton 第 9–10 章：旗流形、Schubert 簇與 Schubert 多項式。置換用一行記法的 tuple，多項式沿用 `tableaux` 的 `{指數 tuple: 係數}`。
- 置換與 Bruhat 序：`length`、`descents`、`act_right`、`w0`、`perms`、`rank_matrix`、`bruhat_leq`（秩矩陣）、`bruhat_leq_tableau`（前綴排序）、`bruhat_chain`（Fulton Lemma 11 的下降鏈）、`bruhat_levels`、`bruhat_covers`
- 胞腔與點數：`flag_cell_pattern(w)`、`grass_cell_pattern(lam, r, n)`（列梯形矩陣的 1／*／0）、`gaussian_binomial`、`q_factorial`、`count_subspaces_F2`、`count_flags_F2`（在 $\mathbb F_2$ 上暴力數）
- Plücker：`plucker(A)`（最大子式）、`plucker_from_coords`（由座標還原矩陣）、`flag_ring_dim`；Schubert 計算：`box_product(lam, mu, r, n)`（方框裡的 LR 乘法）、`rectangle_count`
- 多項式：`divided_difference(f, i)`（差商 $\partial_i$）、`schubert_poly(w)`、`schubert_path(m, steps)`、`pipe_dreams(w)`、`pipe_exits`、`pipe_weight`、`grassmannian_partition`、`grassmannian_perm`、`staircase_monomials`、`poly_str`

`python -m combkit.schubert` 是自我檢查：pipe dream 的總和 = Schubert 多項式（$S_4$ 全部）、Grassmann 置換的 Schubert 多項式 = Schur 多項式（$S_5$ 全部）、秩矩陣判準 = 前綴判準 = 下降鏈、$\mathbb F_2$ 上的點數 = $q$-二項式／$q$-階乘、Fulton 書上的例子逐字重現。

### 4.2f `combkit.fulton_appendix`（純 Python）

Fulton 附錄 A（組合的變形）與附錄 B（拓撲的計數）。表用列的 tuple；對偶字母在 $[N]$ 上用 $x^*=N+1-x$。
- 對偶與 evacuation：`dual_word(w, N)`、`evacuation(T, N)`、`evacuation_steps(T, N)`（逐步：剩下的表、空出的格、填進的 $x^*$）
- 行插入：`column_insert(T, x)`（回傳新表、路徑、新格子）、`column_insert_word(w)`
- 陣列與矩陣：`matrix_to_pairs`、`lex_sort`、`antilex_sort`、`rsk_array`、`burge(A)`（Burge 對應 $(R,S)$）、`corners(A)`（四個角的四對表）
- 0-1 矩陣：`dual_rsk(A)`（Knuth 的對應，形狀互為共軛）、`zero_one_matrices(rows, cols)`、`gale_ryser(lam, mu)`
- 鑰匙：`compact_skew(col_lengths)`、`frank_tableau(T, col_lengths)`、`keys(T)` → $(K_-(T),K_+(T))$
- 附錄 B：`betti_from_cells`、`chern_projective(n)`（$(1+h)^{n+1}$）、`plucker_degree(r, m)`

`python -m combkit.fulton_appendix` 是自我檢查：書上 p.184 的 $T^*$、p.209 的 frank 斜表、p.210 的鑰匙逐格重現；四個角的定理 1 用所有 $2\times2$ 矩陣核對；0-1 對應在 511 個 $3\times3$ 矩陣上一一對應；練習 19、20 在 $n\le5$ 全部核對。

### 4.2g `combkit.walk1d`（純 Python）

一維簡單隨機漫步（步長 1、左右各 1/2），J. Rudnick & G. Gaspari《Elements of the Random Walk》第 1 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk1d.py`），用 `tools/sync_from_walk017.py` 複製過來。
- 計數與機率：`walk_count(N, d)`（$C(N,d)$，不同奇偶為 0）、`walk_prob`（`Fraction`）、`distribution(N)`、`biased_prob(N, n, p)`、`all_walks`、`sequences_with_heads`、`positions`；主例 `MAIN_COINS = "HTHHTHTTHH"`
- 模擬（固定種子）：`ensemble(M, N, seed)`、`finals`、`mean`、`mean_square`、`histogram`
- 動差：`moment_exact`（對分布加總）、`moment_formula`（$N,\ 3N^2-2N,\ 15N^3-30N^2+16N,\dots$）、`moment_via_cosh`（$(-2)^k\,\partial_y^k\cosh(y/2)^N$ 的冪級數）、`fourth_moment_pairings`
- Stirling 與高斯：`stirling`、`ln_stirling`、`stirling_ratio`、`gamma_integrand`／`gamma_gauss`（$w^Ne^{-w}$ 與它的高斯近似）、`gauss_density`、`gauss_prob`（含奇偶性的因子 2）、`ln_binom_approx`
- 遞迴與擴散：`pascal_rows`、`diffusion_step`（`Fraction` 精確）、`second_difference`、`heat_kernel`
- 傅立葉與生成函數：`chi_power_coeffs`、`fourier_count`（均分求積，$M>N+|d|$ 時精確）、`saddle_q0`／`saddle_both`（$q=0$ 與 $q=\pm\pi$ 兩個鞍點）、`gz0_coeffs`、`d2_gen_coeffs`、`d4_gen_coeffs`、`cauchy_coeff(f, N, r, M)`（圓周上數值取係數）、`walk_saddle_z`
- 最陡下降（$f=(z_c-z)^{-\alpha}$）：`sd_exact`、`sd_approx`（書上 (1:S-15) 的更正版，前因子 $1/\sqrt{2\pi}$）、`sd_book`（照印的 $1/2\pi$）、`sd_saddle_exact`、`sd_p`、`sd_p2`、`sd_full`

`python -m combkit.walk1d` 是自我檢查：公式 = 暴力列舉、三種動差算法一致、書上的數字（$\ln10!=15.1044$、15.0961…）與 walk017 場景卡上的每個數字。

### 4.2h `combkit.walk_gf`（純 Python）

格子上的隨機漫步與生成函數，《Elements of the Random Walk》第 2 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_gf.py`），用 `tools/sync_from_walk017.py` 複製過來。不依賴 `walk1d`。
- 生成函數入門：`exp_coeffs`、`legendre_gf`（由 $1/\sqrt{1-2xz+z^2}$ 展開，`Fraction` 精確）、`legendre_bonnet`、`poly_int11`、`ex22_integral`、`ex23_coeffs`、
  `bessel_series`／`bessel_integral`、理想氣體 `ideal_gas_Z`／`grand_probs`、會消失的走路人 `mortal_walkers`／`mortal_trail`／`mortal_origin_exact`（$(1-z)\frac2\pi K(z)$，`ellipk` 用 AGM）
- 格子計數：`SQ_STEPS`、`SC_STEPS`、`BCC_STEPS`、`FCC_STEPS`（整數座標）、`count_grid(N, steps)`（遞迴「最後一步來自鄰居」）、`count_sq`（兩枚硬幣公式）、
  `returns_sq`／`returns_sc`／`returns_bcc`／`returns_fcc`、`uv_walk`、`trail_sq`
- 傅立葉：`chi_1d`／`chi_sq`／`chi_sc`／`chi_bcc`／`chi_fcc`、`chi_from_steps`、`CHI_SMALL_Q`、`CELL_VOLUME`、`BIPARTITE`、`g_fourier`、`fourier_count_sq`／`fourier_count_sc`（均分求積，精確）、`chi2_terms`
- 高斯與普適性：`gauss_sc`（書上 (2.20)）、`gauss_sq`、`gauss_lattice(N, kind, r2, parity)`（含單胞體積與奇偶因子）、`asym_returns_sq`／`asym_returns_1d`、
  `eq222`（(2.22) 的 $q^4$ 修正）、`coin_quartic`
- 離開格子：`p_gauss`、`p_cauchy`、`chi_cauchy`（$e^{-qr_0}$）、`P_cauchy_N`、`cauchy_radial_median`、`steps_gauss_3d`／`steps_cauchy_3d`、`median_distance`
- 擴散：`prob_grid_sq`、`diffuse_sq`、`variance`、`ex28_sides`、`g_cont`／`g_cont_book`（(2.38) 的正負號）、`G1_exact`／`G1_continuum`
- 補充：`gauss_moment`、`gauss_gen`、`square_complete`、`reciprocal`、`SC_BASIS`／`BCC_BASIS`／`FCC_BASIS`、`roots_sum`、
  `series_coeffs_power`、`coeff_power_exact`／`coeff_power_asym`、`coeff_two_poles`、`coeff_log_exact`／`coeff_log_asym`（(2:S3-11) 的更正版）、`digamma`

`python -m combkit.walk_gf` 是自我檢查：遞迴＝兩枚硬幣＝傅立葉積分、書上的數字、16 處筆誤的更正與 walk017 第 2 章場景卡上的每個數字。

### 4.2i `combkit.walk_rec`（純 Python）

回家（常返）、維度與走過幾個點，《Elements of the Random Walk》第 3 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_rec.py`），用 `tools/sync_from_walk017.py` 複製過來。依賴同資料夾的 `walk_gf`。
- 到訪次數的生成函數：`walks_1d`、`visit_counts_1d`（暴力列舉）＝`visit_counts_gf_1d`（(3.9)–(3.10) 的閉式）、`visit_poly`／`counts_from_poly`（(3.1)(3.2)）、`virial_terms_1d`（Fig. 3.1 的排容）
- 一維：`never_return_1d`（$2\binom{N-1}{\lfloor(N-1)/2\rfloor}$）＝`never_return_1d_brute`＝`never_return_1d_gf`、`escape_1d`、`escape_1d_asym`（$\sqrt{2/\pi N}$）／`escape_1d_book`、
  偏硬幣 `escape_biased`／`escape_biased_series`／`escape_biased_dp`
- 各維度：`returns_count(kind, N)`（kind = '1d'、'sq'、'sc'、'bcc'、'fcc'）、`return_probs`、`first_return_probs`、`no_return_probs`（更新方程）、
  `never_return_dp`（逐步遞迴核對）、`return_count_dist`、`no_return_sq_book`／`no_return_sq_refined`、`K_d`、`tail_integral`、`G0`（$G(z;0)$；sc 用 Bessel 積分）、`G0_series`
- Watson 積分：`watson_gamma`（Gamma 閉式）、`watson_sc_bessel`、`watson_bcc_series`、`watson_fcc_series`、`G_C`、`escape`、`f_m`、`sc_correction_coeff`（練習 3.6）
- 足跡：`new_site_probs_brute`、`distinct_sites_brute`、`distinct_sites_1d`（閉式）、`distinct_sites`（生成函數）、`footprint_1d_asym`／`footprint_2d_book`／`footprint_2d_mw`
- 模擬（固定種子）：`simulate_first_returns`、`simulate_return_counts`、`simulate_distinct`、`trail`；`SIM_*`、`TRAIL_*` 是 `check()` 鎖定的結果

`python -m combkit.walk_rec` 是自我檢查：生成函數＝暴力列舉＝逐步遞迴、Watson 積分三種算法、書上的數字、筆誤的更正與 walk017 第 3 章場景卡上的每個數字（約 2 秒）。

### 4.2j `combkit.walk_bc`（純 Python）

邊界條件、穩定態與靜電類比，《Elements of the Random Walk》第 4 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_bc.py`），用 `tools/sync_from_walk017.py` 複製過來。依賴同資料夾的 `walk_gf`、`walk_rec`。
約定：吸收牆是「踩到就消失」的格點；一面牆在 0（鏡像點 $-x_0$），兩面牆在 0 與 $L$（內點 $1..L-1$，特徵值 $2\cos(m\pi/L)$）。
- 反射原理：`walks`、`survivors_brute`＝`survivors_image`（$C(N;x-x_0)-C(N;x+x_0)$）＝`survivors_dp`、`reflect_first`／`reflection_pairs`（Fig. 4.1 的一一對應）、
  `survival_prob`／`survival_asym`（$x_0\sqrt{2/\pi N}$）、三維吸收平面 `survival_3d`／`survival_3d_dp`／`survival_3d_asym`、`wall_reduction_1d`（(4.8) 的一維版）
- 兩牆之間：`box_dp`＝`box_modes`（正弦模態）＝`box_images`（無窮多個鏡像，Ex. 4.3）、`box_eigs`、`box_gf`／`box_amplitude`（(4.12)(4.13) 的更正）、
  `box_total`、`box_survival`、`box_total_modes`、`box_total_cont`（(4.19)）、`box_long`、`box_survival_long`、`decay_rate`、`odd_inv_squares`（$\pi^2/8$）
- 球與立方體：`j0`、`j0_zero`、`sphere_weight`（$8R^3/\pi n^2$）／`sphere_weight_numeric`、`sphere_total_cont`、`sphere_rate`、格子球 `ball_rate`、
  立方體 `cube_total`／`cube_total_modes`（Ex. 4.5）、`PREFACTOR`
- 穩定態：`laplacian_coeff`（1/6）、`tent_lattice`＝`tent_formula`（(4.49)）、`absorb_left`（賭徒破產，Ex. 4.6）、`wall_fluxes`；點源 `time_integral`（$1/4\pi Dr$）、
  sc 格子 Green 函數 `sc_green`（$\to3/2\pi r$）；鏡像源 `conc_wall`／`plane_flux`；吸收球 `sphere_flux`、`kelvin_image`、`hit_prob`（$r_0/R$）、`cell_rmax`
- 書外補充：圓盤電容 `disk_capacitance`（$2a/\pi$）、Berg & Purcell (1977) `berg_purcell`（$na/(na+\pi r_0)$）、`bp_half`；§4.3 的吸收邊界 `half_counts`／`half_gf`（Catalan 數）

`python -m combkit.walk_bc` 是自我檢查：鏡像＝暴力列舉＝逐步遞迴＝模態和、穩定態的精確解、Green 函數、16 處筆誤的更正與 walk017 第 4 章場景卡上的每個數字（約 1 秒）。

### 4.2k `combkit.walk_var`（Python＋numpy／scipy.special）

隨機漫步的變形，《Elements of the Random Walk》第 5 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_var.py`），用 `tools/sync_from_walk017.py` 複製過來。依賴同資料夾的 `walk_bc`（→ `walk_gf`、`walk_rec`）、numpy 與 scipy.special。
約定：密度一律標準歸一化（書上二維是 $2\pi\rho$、三維是 $2\pi^2\rho$）；電報方程的 $\gamma$ 是反轉率，二三維的 $\gamma$ 是轉成隨機新方向的速率。
- 偏向漫步：`biased_pmf`＝`biased_dp`、`biased_gauss`（(5.7)）、`taylor_residual`（(5.3) 的號）、`diffusion_D`（Ex. 5.1：$2p(1-p)$ vs $\frac12$）、`width_ratio`、
  「先走一步再亂走」`shift_pmf`／`shift_dp`／`shift_fourier`（(5.9)(5.10)）
- 有慣性的漫步：`persistent_dp`（左右兩個狀態，Fraction）／`persistent_np`、`persistent_var`＝`persistent_var_formula`、`persistent_var_asym`（$NQ/(1-Q)$）、`persistence_length`、
  `T_matrix`／`T_total`、`eigs`（(5.17)(5.18)）、`eig_residual`（(5.37) 的更正）、`total_fourier`＝`total_eig`、`recursion2_residual`（(5.35)）、`gauss_persistent`
- 電報方程：`tele_cont`（更正後的 (5.61)，$I_0$、$I_1$）、`tele_spike`（$e^{-\gamma t}/2$）、`tele_mass`、`tele_var`、`tele_fourier`＝`tele_fourier_numeric`＝`talbot(tele_laplace)`、
  `tele_gauss`（$D=v^2/2\gamma$）、`tele_lattice`（Ex. 5.3 的格子極限）
- 二三維的「跑一跑、轉個彎」：`angle_avg_2d`／`angle_avg_3d`（(5.77)(5.93)）、`sm_residual`（Ex. 5.7）、`rt2_cont`／`rt2_radial`／`rt2_mass`／`rt2_hankel`（更正後的 (5.84)+(5.89)）、
  `rt3_pole_book`／`rt3_cut_book`／`rt3_once_book`（(5.101)(5.103)）、`rt3_radial`（Fig. 5.6）、`rt3_mass`、`rt3_fourier`＝`talbot(rt3_laplace)`、固定種子的模擬 `rt_mc`、`radial_fraction`
- 連續時間漫步：`survival`、`AB`（$A=1-\lambda B$）、`steps_prob`＝`steps_prob_laplace`（$A^mB$）＝`steps_prob_conv`（卷積）、`mean_steps`、`ctrw_lattice`＝`ctrw_bessel`（Ex. 5.8：$e^{-p_0t}I_x(p_0t)$）、
  `ctrw_gauss`（(5.133)）、`ctrw_transform`（(5.128)）、`mortal_total`＝`mortal_series`＝`mortal_laplace`（Ex. 5.9）；通用的數值 Laplace 反轉 `talbot(F, t)`

`python -m combkit.walk_var` 是自我檢查：精確＝逐步遞迴＝轉移矩陣、電報方程與二三維的公式＝Laplace 反轉＝模擬、總數守恆、10 處筆誤的更正與 walk017 第 5 章場景卡上的每個數字（約 1 秒）。

### 4.2l `combkit.walk_shape`（Python＋numpy／scipy.integrate）

隨機漫步的形狀，《Elements of the Random Walk》第 6 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_shape.py`），用 `tools/sync_from_walk017.py` 複製過來。不依賴其他 combkit 模組（numpy、scipy.integrate）。
約定：$N$ 步走過 $N+1$ 個點，$\overleftrightarrow T=\frac1{N+1}\sum_j(\vec r_j-\bar r)(\vec r_j-\bar r)^{\mathsf T}$，特徵值＝主迴轉半徑的平方 $R_i^2$；模擬一律固定種子。
- 迴轉半徑張量：`tensor`（定義式）＝`tensor_pairs`（兩兩差）＝`tensor_kramers`（Kramers 恆等式 (6.21)，`kramers` 是精確的 $a_{\alpha\beta}$，含對角線）、`tensor_frac`（整數點的精確值）、
  `principal`、`rotate`、`extent`（$\vec n\cdot T\cdot\vec n$）、主例 `small_example`（「東東北西」：$\frac35,\frac15$）、`cloud`／`fig63`／`fig64`（Fig. 6.1、6.3、6.4 的固定種子漫步）
- 非球度：`asph`（個別值）、`A_d`＝$\frac{4+2d}{4+5d}$、`moments_exact`（高斯：Wick 配對；格子：再加四次累積量）＝`moments_enum`（暴力列舉）、`A_exact`（有限 $N$）、
  `Aprime`（(6.17)）＝`Aprime_laplace`＝`Aprime_closed`（$\zeta(3)$）、`sim3d`（三維 2 萬條 100 步）
- 無限多維：`staircase`、`eigs_exact`（$\frac1{4(N+1)\sin^2(k\pi/2(N+1))}$）vs `eigs_asym`（$\frac{N+1}{\pi^2k^2}$）、`sine_mode`、`staircase_asph`（→ $\frac25$）、`orth_fraction`（Ex. 6.2）
- 高斯步伐與 1/d 展開：`wick_pairings`／`wick_count`（$(2m-1)!!$）、`sim_gauss_moments`、`sum_identity`（$\frac34$）、`mean_1d`／`var_1d`（(6.61)(6.62)）、
  `eig_sim(d)`（$d=100,20,3$ 的特徵值模擬，$\xi=U^{\mathsf T}\eta$ 的技巧）、`onedee_table`
- Fixman：`gamma_density`（(6.73)）、`fixman_K`（閉式 $[\sqrt z/\sinh\sqrt z]^{d/2}$）＝`fixman_K_product`、`fixman_density`（數值反轉）、`fixman_mean_var`＝`fixman_moments_numeric`；`A_saw`（(6.77)）

`python -m combkit.walk_shape` 是自我檢查：每個結果兩種以上獨立算法、書上的數字（0.394274、(6.63)–(6.71) 的比值）、12 處筆誤的更正與 walk017 第 6 章場景卡上的每個數字（約 2 秒）。

### 4.2m `combkit.walk_saw`（Python＋numpy／scipy）

路徑積分與自迴避，《Elements of the Random Walk》第 7 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_saw.py`），用 `tools/sync_from_walk017.py` 複製過來。不依賴其他 combkit 模組（numpy、scipy.integrate／special／optimize）。
約定：二維正方格子（E、N、W、S）；「嚴格」自迴避＝$N+1$ 個點兩兩不同；平均場的格點數 $\mathcal N$（程式裡 `Ns`）、步數 $m$。
- 路徑和：`count1`、`count2`、`count_brute`、`split_terms`／`split_sum`（(7.3)：一維 $1\cdot1+2\cdot2+1\cdot1=6$、二維 $36$）、`multi_split`（(7.5)）、`gauss1`／`gauss_conv`／`conv_ratio`（Ex. 7.1：變異數相加；書上的前因子少了 3）
- 自迴避漫步：`saw_enum`（DFS＋對稱，$c_N$ 與 $\sum R^2$，$N\le14$ 約 0.2 秒）＝`SAW_KNOWN_C`（文獻）＝`saw_brute`（暴力列舉 $4^N$ 條）、`saw_R2`、`saw_frac`、`saw_list`（逐條點列）、`sample_saw`／`sample_rw`、`two_nu_eff`（→ $\frac32$）
- $1-f$ 展開：`ie_terms`（逐條漫步）＝`ie_terms_subsets`（逐個子集）＝`order1_formula`（一階）、`ie_partial`（部分和，最後＝$c_N$）、`weighted_count`、`order2_types`；
  弧線圖 `matchings`、`chord_count`（$(2k-1)!!$）、`is_irreducible`／`irreducible_count`（1, 2, 10, 74）、`dyson_counts`（$1/(1-I)$ 的係數＝$(2k-1)!!$）
- 生成函數與 Dyson：`convolve`、`gf_check`、`fourier_conv_error`、`dyson_partial`／`dyson_closed`（$1/(g_0^{-1}-\Sigma)$）、`sigma1_1d`（一維一階自能 $1/\sqrt{1-4z^2}$）、`g_dyson_1d`
- 冪次與平均場：`power_exp`（$1-l(4-d)/2$）＝`power_exp_count`、`loop_slope`（數值積分量冪次）、`G_MF`（erfcx 閉式）＝`G_MF_quad`、`cmf_ratio`（(7.59) 更正版，quad／梯形／拋物柱函數／鞍點）vs `cmf_asym`（$e^{-fm^2/\mathcal N}$）、`ex75`
- Flory：`flory_F`、`flory_Rstar`（解析）＝`flory_Rstar_num`（數值最大化）、`flory_slope`（＝$\frac3{d+2}$）、`NU_QUOTED`

`python -m combkit.walk_saw` 是自我檢查：每個結果兩種以上獨立算法、13 處筆誤的更正與 walk017 第 7 章場景卡上的每個數字（約 1 秒）。

### 4.2n `combkit.walk_scale`（Python＋numpy／scipy）

標度入門，《Elements of the Random Walk》第 8 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_scale.py`），用 `tools/sync_from_walk017.py` 複製過來。用相對 import 引用同套件的 `walk_gf`（Cauchy、冪級數係數、回原點條數）與 `walk_saw`（$c_N$、Flory）。
約定：長尾巴（Pareto）步伐 $P(\text{步長}>r)=r^{-y}$（三維 $p\propto r^{-3-y}$，書上 $r^{-5+\delta}$，$y=2-\delta$）；§8.1.2 的 $a=\chi(0)$（簡單立方 $a=2d$、$b=1$、$z_c=\frac1{2d}$）；重新標度 $s'=l_2s$、$\vec r\,'=l_1\vec r$。
- 非高斯：`moment2_trunc`（(8.2) 截斷的二階動差）、`chi_pareto`（QAWF 數值）vs `chi_small`（$1-Bq^y$，`B_pareto` 閉式＝`B_pareto_quad`）、`median_width`／`width_pred`／`width_slope`（固定種子；寬度 $\propto N^{1/y}$）、`trail`、
  `levy_density`（數值逆 Fourier）vs `levy_tail`（換變數的尾巴 $r^{-5+\delta}$；$y=1$ 即第 2 章 Cauchy）、`subst_integral`／`subst_closed`（(8.6)(8.7)，$w=kr$）
- 標度形式：`g_lattice`、`g_88b`／`g_88c`、`collapse`／`collapse_limit`（一維 $(z_c-z)g\to\frac2{4+x^2}$）、`g_general`／`gamma_numeric`（$\gamma=(2-\eta)\nu$）、
  `G_lattice_sc`（三維格子 $G$ 的 Bessel 積分）vs `G_yukawa`（$\frac3{2\pi}e^{-r/\xi}/r$）、`xi_88c`／`xi_lattice`
- 最陡下降：`ratio_exact`（一維二項式）vs `sd_ratio`（(8.20) 更正版；`book=True` 照書上少了 $1-\nu$）、`sd_delta`＝`sd_delta_num`、`fisher`（$\frac1{1-\nu}$）、`gamma_exp`、`growth`／`p_numeric`、
  `ratio_method`（二維 $c_N$ 的比值法：$\hat\mu\to2.638$、$\hat\gamma\to\frac{43}{32}$）、`theta`、`return_slope`
- 重新標度：`kin_factor`、`int_factor`（＝`rescale_check` 的數值）、`sa_amplitude`（$N^{(4-d)/2}$，上臨界維度 4）、`flory_l1`＝`flory_l1_num`、`nu_rescale`（$\frac3{d+2}$）、`ex88_*`（冪次型交互作用）

`python -m combkit.walk_scale` 是自我檢查：每個結果兩種以上獨立算法、7 處筆誤的更正與 walk017 第 8 章場景卡上的每個數字（約 1 秒）。

### 4.2o `combkit.walk_crit`（Python＋numpy／scipy）

標度與臨界現象，《Elements of the Random Walk》第 9 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_crit.py`），用 `tools/sync_from_walk017.py` 複製過來。用相對 import 引用同套件的 `walk_saw`（$c_N$、$\langle R^2\rangle_N$、重合點對數、配對）與 `walk_scale`（$\mu$、$\gamma$、$\eta$、$\nu$ 的引用值）。
約定：$J=k_B=1$、$K=J/T$、無因次外場 $H=h/T$；平均場 $m=\tanh(m/T)$（$T_c=1$）、二維 Ising $T_c=2/\ln(1+\sqrt2)$；Landau $F=\frac r2m^2+\frac u4m^4-hm$；O($n$) 小圖用 $|\vec S|^2=n$、權重 $\prod(1+K\vec S_i\cdot\vec S_j)$。
- 自發磁化：`mf_m`＝`mf_m_newton`、`onsager_M`（＝低溫展開 `onsager_lowT`）、`beta_slope`（0.5000、0.1250）、`ms_table`
- 對稱破缺：`rotation_spins`／`rotation_check`（(9.4) 在整體旋轉下不變）、`ising_exact`（逐一列舉：$Z$、$\langle|m|\rangle$、$P(M)$、外場）＝`ring4_closed`、
  `mc_run`／`flip_tau`／`mc_trace`／`mc_snapshot`（棋盤式 Metropolis、固定種子；$T=2$ 的翻轉等待 133、608、2151 sweep）、`mc_check_L4`（取樣器＝精確列舉）
- 指數與標度律：一維 Ising `ising1d_G`／`ising1d_xi`／`ising1d_chi`（＝轉移矩陣＝環的列舉）、`EXPONENTS`（平均場、二維 Ising、二維自迴避）＋`laws`（Rushbrooke、Fisher、hyperscaling、Widom；Fraction）、
  `book3d_alpha`（書上 Table 9.1 的 $\alpha=0.12$ 應為 0.19）、`widom_exponents`／`widom_G`／`widom_slopes`（Ex. 9.1–9.2）
- GLW：`fourier_check`（(9.30)–(9.32) 的離散版）、`landau_F`／`landau_min`、`ex94_*`（$\beta=\frac12$、$\Delta=\frac32$、$\alpha=0$）、`zmf_avg`（(9.34) 的有限體積積分）
- §9.4：`saw_gf`＝`saw_gf_realspace`（$\xi^2$ 的兩種列舉）、`coef_log`＝`coef_log_rec`、`coef_ratio`、`r2_ratio`（(9.41) 是 $z_c^{2\nu}$，書上印成 $z_c^\gamma$）
- $n\to0$：`hs_rhs`（(9.46)）、`hs_book72`、`visit_weight`／`intersection_factor`（每相交一次 $e^{-u/2}$）、`weak_saw_count`、`heat_kernel*`、`gauss_I`＝`gauss_I_quad`、`gauss_moment`、`wick4`、`wick_count`、
  `trotter_err`／`trotter_global`（Ex. 9.6）、`ex95_check`、`ladder_G`（＝`ladder_G_ising`、`ladder_G_xy`；$n\to0$＝`ladder_saw_sum`，只剩自迴避路徑）

`python -m combkit.walk_crit` 是自我檢查：每個結果兩種以上獨立算法、10 處筆誤的更正與 walk017 第 9 章場景卡上的每個數字（約 4 秒，其中 Metropolis 約 2 秒）。

### 4.2p `combkit.walk_cond`（Python＋numpy／scipy）

O($n$) 模型的平均場、自旋波與關在盒子裡的漫步，《Elements of the Random Walk》第 10 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_cond.py`），用 `tools/sync_from_walk017.py` 複製過來。用相對 import 引用同套件的 `walk_saw`（無限平面的 $c_N$、全部 $4^N$ 條漫步、(7.57) 的 $G_{MF}$）、`walk_crit`（$z_c=1/\mu$、指數表）、`walk_rec`（sc 的 $G_c$、回家機率）與 `walk_gf`（sc 的 $C(N;0)$）。
約定：$e^{-r}=z/z_c$；平均場 $f(\lambda)=e^{-N_0(\frac r2\lambda+\frac u4\lambda^2)}$、$\Gamma_{MF}=\frac{N_0}2\int_0^\infty f$（書上的 $\frac{N_0}4$ 是筆誤）；環面格點 $s=x+Ly$。
- 環面上的自迴避漫步：`torus_counts`（位元遮罩 DFS、8 重對稱；＝`torus_counts_plain`＝`torus_counts_brute`）、`torus_vs_inf`、`torus_first_diff`（$N=L$ 時少 4 條）、`GRAY_Q4`（$16\times5712=91392$）、
  畫圖用的 `torus_walks`／`torus_example`／`torus_segments`（越界的一步拆兩個半段）／`torus_lift`、`torus_mean_N`
- 平均場：`mf_I`（erfcx＝quad）、`gamma_mf`、`gamma_n`（有限 $n$ 的徑向積分＝拋物柱函數；$n\to0$ 的極限與極限順序）、`s_integral`（$2\pi i$ 不是 $i\pi$）、漸近、`fn3_check`（＝`walk_saw.G_MF`）
- 自旋波：`wall_M`／`wall_bvp`（牆邊的 $\tanh$）、`ex101_grad_check`、`ex102_coeffs`、`sw_X`／`sw_Xinv`／`sw_lambdas`（Fraction 精確）、`sw_check`、`mode_table`（$n\to0$ 時橫向重數 $-1$）、`gauss_n2_check`
- 相關洞：`D_box`（三維週期格子，FFT）＝`D_direct`＝`D_inf`＋`box_shift`（Ewald 常數 2.837297）、`coulomb_fit`（$-D\approx c_0-c/\Delta$、$c\approx\frac1{4\pi}$）、`D_cont`（Ex. 10.5）、`case1_check`；
  §10.4：`sc_prob`／`sc_visits`（＝`sc_visits_dp`）vs `erfc_form`（(10.83) 更正版）、`green_sc`（$\to\frac3{2\pi R}$）
- Redner–Reynolds：`torus_moments`／`torus_rho`（＝`torus_rho_fd`）／`torus_peak`／`torus_curve`、`mf_rho`／`mf_chi`（$\beta'=1$、$\gamma'=\alpha'=0$）、`rr_exponents`＝`table101`、`rr_laws`、`dc_mean_field`（$d_c=4$）

`python -m combkit.walk_cond` 是自我檢查：每個結果兩種以上獨立算法、11 處筆誤的更正與 walk017 第 10 章場景卡上的每個數字（約 4 秒，其中 5×5 環面的列舉約 1 秒）。

### 4.2q `combkit.walk_rg`（Python＋numpy／scipy）

標度、碎形與重整化，《Elements of the Random Walk》第 11 章用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_rg.py`），用 `tools/sync_from_walk017.py` 複製過來。用相對 import 引用同套件的 `walk_saw`（四方向步伐、全部 $4^N$ 條漫步、$2\nu_{\rm eff}$）與 `walk_crit`（$\nu$ 的已知值、指數表、Landau 的根、Widom 指數）。
約定：Landau 用 $\frac r2m^2+\frac u4m^4-hm$（$u=1$）；一維 Gaussian 鏈 $\sum[\frac r2S_i^2+\frac a2(S_i-S_{i+1})^2]$、$a=1$；抽格 $b=2$、連續 $b=e^\ell$；指數用 Fraction。
- 碎形：`koch_points`／`koch_counts`、`cantor_intervals`／`cantor_counts`、`mass_table`／`mass_radius_dim`（$n(R/3)=n(R)/4$ 精確 → $\ln4/\ln3$）＝`box_dim`（盒子計數）＝`replica_dim`（(11.6)）；
  `rw_rg2`＝`rg2_exact`（$\langle R_g^2\rangle=\frac{N(N+2)}{6(N+1)}$）、`rw_df`（→ 2）、`saw_df`（$1/\nu$）、`saw_df_eff`
- 被動的重整化：`landau_m`（＝二分法＝`walk_crit.landau_min`）、`landau_scaled`／`collapse_curve`（$y+y^3=x$ 的塌縮）、`mf_rescale_check`、`gauss_rescale_check`、`gauss_exponents`（Ex. 11.4）
- 一維 Gaussian 鏈：`square_check`（完成平方＝數值積分）、`decimate_gauss`、`rr_gauss`（$r''=r(r+4a)/a$）、`hh_gauss`、`ring_decimation_check`（環的行列式；(11.31) 應為 $N/4$ 次方）、`cobweb_points`、`small_r_table`、`chain_lambdas`
- 一維 Ising：`ising_A`／`ising_Kp`、`ex118`（加外場）、`ring_Z_brute`＝`ring_Z_transfer`＝`ring_Z_decimation`（(11.63) 的指數應為 $N/2^j$）、`x_map`（$x'=\frac{2x}{1+x^2}$）、`xi_ratio`（$\to\frac12$）、`free_energy_ratio`（$F\to-NJ$）
- 動量殼與流：`alpha_q`（$a(1-\cos q)$；書上 $-2a(1-\cos q)$）、`fourier_check`、`shell_rg`＝`shell_rg_steps`、`flow`＝`flow_numeric`、`pde_residual`（(11.24) 滿足 (11.87)）、`char_check`
- 不動點與指數：`toy_eigen`／`toy_flow`（＝矩陣指數）、`rg_exponents`（(11.105)–(11.108)）＝`rg_exponents_widom`＝`walk_crit.EXPONENTS`（Gaussian $d=4$、二維 Ising、二維自迴避）、`xi_from_flow`＝`xi_from_flow_numeric`

`python -m combkit.walk_rg` 是自我檢查：每個結果兩種以上獨立算法、8 處筆誤的更正與 walk017 第 11 章場景卡上的每個數字（約 1 秒）。

### 4.2r `combkit.walk_eps`（Python＋numpy／scipy）

重整化群的 ε 展開，《Elements of the Random Walk》第 12 章（全書最後一章）用的數學；正本在 walk017 專案
（`C:\mypython\manim\projects\walk017\codes\utils\walk_eps.py`），用 `tools/sync_from_walk017.py` 複製過來。用相對 import 引用同套件的 `walk_saw`（$c_N$、`count2`）、
`walk_crit`（$\nu$ 的已知值、$z_c$、Flory 差距）與 `walk_rg`（$\ln4/\ln3$、複本公式）。
約定：遞迴方程把書上 (12.37) 的 $d$ 併進 $u$：$\frac{dr}{d\ell}=2r+\frac{Au}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{Bu^2}{(1+r)^2}$，$A=4(n+2)$、$B=4(n+8)$（Ising 12、36）；分數用 Fraction。
- 積掉自由度：`toy_Heff`（兩變數玩具，閉式＝數值積分）、`toy_coeffs`（$r'=r+\frac{12u}{1+r}$、$u'=u-\frac{36u^2}{(1+r)^2}$）＝`toy_coeffs_numeric`、`toy_Z`（$Z=A\,Z'$）；
  `shell_fraction`（$1-e^{-d\delta\ell}$）＝`shell_fraction_lattice`（數格點）、`rescale_factors`（(12.33) 的筆誤）
- 配對與連結叢集：`matchings`、`type_census`／`diagram_census`（列舉）＝`pair_factor`（(12:S-4)）＝`general_factor`（(12:S-6)，`ring_aut` 的 $f_3=6$）、`ring_weight`、`egf_check`；
  `lc_average`（數值積分＝所有圖＝連通圖的指數＝閉式）、`lc_order`、`lc_cumulants_rings`＝`lc_cumulants_moments`、`lc_multi`、`on_counts`（einsum：$2n+4$、$4(n+8)$）
- 遞迴與不動點：`beta`／`beta_book`、`fixed_point`（Fraction 精確，Ising $(-\frac\epsilon{6+\epsilon},\frac\epsilon{(6+\epsilon)^2})$）＝`fixed_point_numeric`、`jacobian`＝`jacobian_fd`、`eigen`、`lambda_r_closed`、
  `flow_line`／`flow_starts`、`r_critical`（二分法；Ex. 12.4）、`u_analytic`＝`u_numeric`（(12.41)）、`int_1248`（(12.48) 的筆誤）、`t_closed`＝`t_numeric`
- $\nu$：`nu_first`（$\frac12+\frac{(n+2)\epsilon}{4(n+8)}$）、`nu_unexpanded`、`nu_eigen`＝`nu_from_flow`、`flory_series`（$\frac3{6-\epsilon}$ 的展開）、`NU_3D_CITED`（只引用）
- 圖解法：`g_numeric`（(12.86)＝(12.39)）、`pc_exponent`（$u^2r^{d-3}$；書上 $2d-3$）＝`pc_numeric`（實空間 $\int G^3$ 的數值冪次）
- 全書回顧：`BOOK_RECAP`、`recap_numbers`、`BOOK_SCENES`（共 87 幕）、`scene_range`

`python -m combkit.walk_eps` 是自我檢查：每個結果兩種以上獨立算法、11 處筆誤的數值與 walk017 第 12 章場景卡上的每個數字（約 1.2 秒）。

### 4.2s `combkit.young_extras`（Python＋numpy）

youngart（Fulton《Young Tableaux》全書 60 張 4K 藝術桌布）裡可重複使用的計算：`patience(w)`（patience sorting 的堆號與一組最長遞增子序列）、`slide_path(T)`（拿掉左上角後一次 jeu de taquin 的洞的路徑與新表；evacuation 每一步）、`bump_routes(w)`（列插入的擠壓路徑與倒序行插入的路徑，得到同一個 P）、`shadow_lines(w)`／`four_corners(w)`（Viennot 影子線；從四個角各做一次，條數是 P 的第一列長／第一行長）、`rsk_rows(w)`（bisect 版 R-S-K，n 上萬也快）、`spectral_xy`（圖的譜嵌入）、`pipe_route(D, m, i, cs)`（pipe dream 的水管折線）、`ruling`／`line_meets_hyperboloid`（單葉雙曲面的兩族直線、四條線問題的兩條橫截線）、`ellipse_curve`／`line_curve`／`cubic_curve`＋`meet`（平面曲線的實交點，Bézout 的 d·e）。自我測試：`python -m combkit.young_extras`。

### 4.2t `combkit.sym_mat`（Python＋numpy／scipy／sympy）

Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 1 章（矩陣代數與量子力學的複習）用的數學；正本在 harter_sym 專案（`C:\mypython\manim\projects\harter_sym\codes\utils\sym_mat.py`），用 `tools/sync_from_harter_sym.py` 複製過來。精確值用 sympy（Rational、sqrt、I）。
約定：分析器轉 $\theta$ 的轉換矩陣是書上 (1.1.3) 的 $U(\theta)=\begin{pmatrix}\cos\frac\theta2&-i\sin\frac\theta2\\-i\sin\frac\theta2&\cos\frac\theta2\end{pmatrix}$；久期多項式取首一的 $\det(\varepsilon\mathbf1-M)$。
- 分析器與 Axiom 4：`U`、`T`（$T(\beta\leftarrow\alpha)$）、`T_real`、`pass_fraction`、`analyzer_chain`（密度矩陣；只留一通道／相干合併／量了走哪條）＝`axiom4`（振幅求和）、`analyzer_chain_num`、
  `series_table`（$0\to90\to180^\circ$：0、0、$\frac12$、$\frac14$）、`problem_111`、`zeno`（$\cos^{2n}\frac\pi{4n}$）
- bra／ket 與換基底：`elementary`、`ket_primed`、`change_basis`（$T\mathcal HT^\dagger$）、`outer`、`dag`、`is_unitary`、`to_np`
- 對角化：`secular`、`eigen_mult`、`minimal_poly`、`poly_at`（Hamilton–Cayley）、`hc_factors`（$\mathfrak p_j$）、`projectors`（$\mathcal P_j$，精確）、`mat_func`（$f(M)=\sum f(\varepsilon_j)\mathcal P_j$）、
  `eig_book`（左右本徵向量、對偶）、`adj_vectors`（伴隨矩陣法）、`joint_projectors`（對易矩陣）、`lagrange_basis`、`singular_values`、`iterate_angles`、`sqrt_all`、`secular_coeffs`（三種算法）
- 散射、密度矩陣、兩態拍頻：`S_EX`、`flux_in_out`、`rho_pure`、`purity`、`H2`、`populations`＝`populations_numeric`＝`psi_t`
- 耦合振子（Fig. 1.4.1）：`force_matrix`、`mass_matrix`、`accel_matrix`（$a=m^{-1}\mathcal F$）、`omega2`、`omega2_generalized`（`eigh(F, m)`）、`osc_analytic`＝`osc_numeric`、`beat_times`、`modes_symbolic`、`spring_force_matrix`
- 連續空間：`d1_matrix`（反對稱）、`fourier_matrix`、`delta_gauss`、`delta_integral`、`delta_deriv_integrals`（(1.5.10) 的符號）、`adjoint_fd`、`resolvent`、`resolvent_adj`、`contour`（圍道積分＝投影算符）、`contour_PN`
- 習題 1.1.1–1.2.8 的資料：`P121`、`P122`、`P124`、`problem_125`、`P126_M`／`P126_N`、`P127_*`、`shift_matrix`

`python -m combkit.sym_mat` 是自我檢查：每個結果兩種以上獨立算法、可數值確認的書上筆誤與 harter_sym 第 1 章場景卡上的每個數字（約 2 秒）。

### 4.2u `combkit.sym_abel`（Python＋numpy／scipy／sympy）

Harter PSDS 第 2 章（Abel 對稱群與對稱表示）用的數學；正本在 harter_sym 的 `codes\utils\sym_abel.py`，用 `tools/sync_from_harter_sym.py` 複製過來；最小方程的投影算符借用同一個套件的 `sym_mat`（相對 import）。
約定：運動方程 $\ddot x=-\mathcal Ax$；單擺環 $\omega^2(k)=b+2a(1-\cos k)$；$\hbar=1$；書上的筆誤（$a\equiv k/m\ell^2$、(2.5.3) 的 $2A+b$、(2.7.9) 少的 $b$、(2.8.1a) 的負號、(2.12.38b) 的 $D$ 符號等）都已更正，原印的寫法另有 `book=True` 或 `*_BOOK` 可比對。
- 群：`ROT3`、`D2_TABLE_BOOK`、`table_from_matrices`／`table_from_rows`、`axioms`、`assoc_failures`（「不是群」的 (2.2.11)：100 組不結合）、`cyclic_table`、`order_counts`、`GROUP_COUNTS`（引用值，64 階 267）、`abelian_count`／`abelian_types`（分割數）
- 投影算符與特徵標：`c2_projectors`（最小方程）＝`c2_projectors_char`、`cn_projectors_minimal`＝`cn_projectors_char`（精確，$n\le4$）、`cn_projectors_np`、`cn_char_table`、`rep_projectors`、`multiplicities`、`d2_projectors_product`＝`d2_projectors_char`、`outer_char_table`（Kronecker）、`abelian_char_table`、`is_cyclic`、點群 `POINT3`／`is_group`／`point_table`
- 耦合擺：`accel2`、`Modal`（$\ddot x=-\mathcal Ax$ 的解析解，給 updater）＝`ode_solve`、`beat_times2`、`detuned_max_transfer`、`lissajous`＝`lissajous_numeric`、`closure`（$t_{\rm LJ}$、拍數）、`ellipse_axes`、`phasors`、`accel3`／`reduce3`（部分約化）、`ring_matrix`／`ring_omega2`、`wave_vector`、`dispersion`、`brillouin`、`alias`、`broken_ring_matrix`＝`broken_bands`（聲學／光學帶）、`broken_standing`、`d2_matrix`／`d2_omega2`
- 波：`swr`、`gallop_psi`、`gallop_envelope`、`gallop_zeros`、`zero_velocity`（在 $c\Delta$ 與 $c/\Delta$ 之間）
- 量子：`rotor_matrix`／`rotor_block`／`rotor_bands`（$C_{12}$ 轉子）、`pair_exact`／`pair_pert`、`tight_fit`、`tb_matrix`／`tb_levels`、`mathieu_edges`＝`mathieu_edges_scipy`、`rotor_edges`、`square_well_G`＝FFT、Kronig–Penney `kp_half_trace`＝`kp_half_trace_ode`、`kp_edges`、`kp_gap_near`、`two_level`、`tunnel_populations`、`lz_formula`＝`lz_numeric`、方井 `s_matrix`＝`s_matrix_direct`、`phase_shifts`、`transmission`、`resonances`、`crossing_step`／`crossing_well`、`s_from_c`／`c_from_s`、`bound_geometric`＝`bound_fd`、`c_barrier`／`c_valley`、`wkb_cluster`、`wkb_vs_exact`
- 習題 2.3.3–2.12.3：`problem_233` … `problem_2123`（2.12.3 用節點計數 `line_levels`＝有限差分 `line4_fd`／`ring4_fd`）

`python -m combkit.sym_abel` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 2 章場景卡與導讀上的每個數字（約 3 秒）。

### 4.2v `combkit.sym_nonabel`（Python＋numpy／scipy／sympy）

Harter PSDS 第 3 章（非 Abel 對稱群）用的數學；正本在 harter_sym 的 `codes\utils\sym_nonabel.py`，用 `tools/sync_from_harter_sym.py` 複製過來；借用同一個套件的 `sym_abel`（平面波轉子、單擺環、外積）（相對 import）。
約定：群乘積 $ab$ ＝ 先 $b$ 後 $a$；$C_n$ 的不可約表示一律 $\mathcal D^{k_n}(r)=e^{+2\pi ik/n}$；$C_{3v}$ 的鏡線 $\sigma_3$ 沿 $x$ 軸、$\sigma_2$ 60°、$\sigma_1$ 120°；書上的筆誤（(3.3.18) 的 $\varepsilon^*r^2$、Fig. 3.4.4 型 2 的向量、(3.4.24a) 的 $\varepsilon$ 號、(3.5.13) 的 ${}^\circ c_g$、NH₃ 的 (3.7.15)–(3.7.24) 與 p. 216 的擬合）都已更正，原印的寫法另有 `book=True` 或 `*_BOOK` 可比對。
- 有限群：`FiniteGroup(mats)`、`FiniteGroup.from_generators(gens)`（任何精確矩陣群）：乘法表、共軛類、正規化子、陪集與雙陪集、換位子群、類代數的結構常數與正則表示、類算子的最小方程、全對易冪等元（逐一分裂）、群代數（`ga_mul`、`ga_dag`…）、正則表示、子群分裂與基本算子 $P^\alpha_{ij}$、由基本算子求不可約表示、特徵標表（類代數＝跡）、頻率公式、相關（subduction）表、同構；`book_order`、`outer_table`（Kronecker）、`cn_chars`、`ga_tex`（群代數元素 → TeX）
- 具體的群：`C3V`、`D3`、`D4`、`C4V`、`D2D`、`D6`、`D3H`、`D3D`、`D2H`、`D4H`、`QUAT`（四元數群）與書上的表（`*_TABLE_BOOK`、`*_CHARS_BOOK`）
- Hamilton 轉：`quat`、`qmul`、`quat_to_rot`、`su2`、`hamilton_arc`、`hamilton_product`（半角弧頭尾相接＝四元數乘積＝兩次反射）、`d3_arcs`、`conjugate_rotation`、`body_product`（本體座標的反序）
- $C_{3v}$：`c3v_split_cv`、`c3v_elementary`、`c3v_irreps`（駐波＝幾何矩陣；圓形 $\mathcal D^{3^\circ}$）、`c3v_regular_T`（正則表示的塊對角化）、三擺的 `pendulum_frequencies`、`pendulum_trace_eigs`；六自由度彈簧質量模型 `six_dof_*`（第一列公式＝彈簧幾何、對稱約化、Fig. 3.4.4 的駐波、圓周運動）
- $D_n$ 位能裡的轉子：`dn_levels`（本徵向量的對稱定標記＝$m\bmod n$）、`dn_level_paths`、`dn_gap_2x2`、`dn_density`、`dn_zeeman`
- 氨 NH₃：`nh3_geometry`、`nh3_force`（12×12）、`nh3_rep`、`nh3_projector`、`nh3_cvectors`、`nh3_reduced`＝`nh3_reduced_formulas`、`nh3_secular`、`nh3_fit`、`nh3_ratio_fit`、`nh3_crossing`、`nh3_scan`、`nh3_modes`
- 習題 3.1.1–3.5.6：`problem_311` … `problem_356`、`d4_square_matrix`

`python -m combkit.sym_nonabel` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 3 章場景卡與導讀上的每個數字（約 3 秒）。

### 4.2w `combkit.sym_oct`（Python＋numpy／scipy／sympy）

Harter PSDS 第 4 章（八面體對稱與誘導表示）用的數學；正本在 harter_sym 的 `codes\utils\sym_oct.py`，用 `tools/sync_from_harter_sym.py` 複製過來；有限群工具全部借用同一個套件的 `sym_nonabel`（相對 import）。
約定：群乘積 $ab$ ＝ 先 $b$ 後 $a$；元素名照書（`r1`…`r4^2`、`R1^2`、`R1`、`R1^3`、`i1`…`i6`；$O_h$ 的 `I.r1`…），$r_k$ 的軸是一個正四面體的四個頂點；**$C_4$ 照本章書上 $\mathcal D^{1_4}(R_3)=-i$**（與 `sym_nonabel` 的 $+i$ 相反，相關表不受影響）；
書上的筆誤（(4.1.7) 的 $c_i^5$、(4.2.18) 的前置因子、(4.2.28) 三個元素、(4.2.36) 的對角、(4.3.46)ₓ、(4.3.63)、(4.3.66)、(4.3.71)、(4.4.22) 兩個力常數、(4.4.33) 的 $j-b$ 等）都已更正，原印的寫法另有 `book=True` 或 `*_BOOK` 可比對。
- 群：`O`、`OH`（＝ $O\times C_i$）、`TD`（$\cong O$）、`TGRP`（$T$）、`TH`；Hamilton 弧 `nomogram_product`、四條體對角線的置換 `diag_perm`、`cycle_str`（$O\cong S_4$）
- 類代數 (4.1.6)、$c_i$ 的冪與最小方程 $c^5-40c^3+144c$、全對易冪等元、特徵標表 `O_CHARS`、`OH_CHARS`、`T_CHARS`（複數）
- 不可約表示 `irrep(α, g, basis)`：四條子群鏈 `"tetragonal"`（$O\supset D_4\supset D_2$）、`"trigonal"`（$D_3\supset C_2$）、`"c4"`、`"c3"`（行進波）；鏈的冪等元分裂 `chain_split`、`d4_idempotents`、`d3_idempotents`；
  子群 `SUB` 與相關表 `correlation(key)`（(4.2.42)–(4.2.48)、$D_{4h}$、$D_{3d}$）、`conjugates`
- 誘導表示 `induced`（(4.3.34)）、`induced_char`、`frobenius(key)`（矩陣的跡＝誘導特徵標公式＝相關表的行）；六點穿隧模型 `six_*`（$H+T+4S$、$H-T$、$H+T-2S$；投影＝陪集分解）；
  外場 `vq_*`（$D_{4h}$）、`vd_*`（$C_{4v}$）、`vb_*`（磁場 $C_4$，Peierls 相位模型）、`v111_*`（(111) 電場）與能階掃描
- 六氟化物 XF$_6$：21 個座標 `XF6_NAMES`、彈簧幾何的力矩陣 `xf6_force`（＝(4.4.22) 三列＋對稱 `xf6_force_from_rows`）、`xf6_rep`、`xf6_classify`、`xf6_projected`、`xf6_numeric`＝`xf6_formulas`、
  $T_{1u}$ 的 2×2 `xf6_t1u_matrix`／`xf6_secular`、`xf6_scan`（Fig. 4.4.5–4.4.6）、`uf6_fit`、`xf6_isotope`（SF$_6$ 1.0172、1.0065）、`xf6_t1u_vectors`、`xf6_modes`
- 量子振動：截斷的升降算子 `ladder`、`ladder_commutator`、`osc_levels`、簡併度 `degeneracy`＝`degeneracy_count`（Pascal 三角形）、`harmonic_states`、`state_count`、`turning_point`
- 習題 4.1.1–4.4.2：`problem_411` … `problem_442`、`xy4_fit`（$XY_4$ 的兩參數擬合）

`python -m combkit.sym_oct` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 4 章場景卡與導讀上的每個數字（約 11 秒）。

### 4.2x `combkit.sym_rot`（Python＋numpy／scipy／sympy）

Harter PSDS 第 5 章（連續轉動群的表示與應用）用的數學；正本在 harter_sym 的 `codes\utils\sym_rot.py`，用 `tools/sync_from_harter_sym.py` 複製過來；譜分解、Hamilton 轉與 $O$ 的特徵標借用同一個套件的
`sym_mat`、`sym_nonabel`、`sym_oct`（相對 import）。約定：角度一律用度；$R(\alpha\beta\gamma)=R_z(\alpha)R_y(\beta)R_z(\gamma)$（實驗室的軸，最右先做）；
$\mathcal D^j_{m'm}=\langle jm'|e^{-i\alpha J_z}e^{-i\beta J_y}e^{-i\gamma J_z}|jm\rangle$（標準 Wigner 相位，基底 $m=j\ldots-j$）；書上的筆誤（(5.4.25)ₓ 的 $\frac{\sqrt6}2$、(5.4.35) 四格、
(5.5.17) 的因子 2、(5.6.26)–(5.6.27)、(5.7.14)、(5.7.26b)、(5.8.6)–(5.8.10) 等）都已更正，照印的另有 `*_BOOK`。
- $O_3$、$O_2$、$D_{\infty h}$：`o2_rot`、`o2_refl`、`o2_standing`、`dinfh_chars`、`dinfh_to_dnh`（＝ $m$ 規則）
- 轉動的參數：`euler_matrix`、`euler_from_matrix`、`euler_alt`、`darboux_matrix`（＝ 五個曲柄）、`euler_to_darboux`（(5.5.2)：$(50,60,70)\to[80.00,33.69,128.68]$）、`darboux_sequence`（Fig. 5.3.5 的 720° 序列）、
  `su2_euler`、`su2_darboux`、`su2_product`（(5.5.10)）、`lab_frame_product`／`body_frame_product`（本體操作的反序）
- 生成元與表示：`GX`、`GY`、`GZ`、`gz_projectors`（譜分解）、`small_commutator`（Fig. 5.4.2）、`jmat(j)`、`wigner_d`（Wigner 公式）＝ `D_expm`（矩陣指數）＝ `D_schwinger`（兩個振子）、
  `wigner_D`、`D_axis`、`d_sym`、`D_from_u`、`d_symmetry_checks`
- Stern–Gerlach：`sg_probs`、`sg_amps`、`mean_jz`、`cone_angle`、`cone_theta`；量子轉子：`rotor_ops`（本體 ＝ $J^*$）、`rotor_levels`（數值）＝ `rotor_formula`、`asym_matrix_formula`
- 球諧：`ylm` ＝ `ylm_from_D`、`addition_theorem`、`rotated_ylm`；分裂：`chi`、`o3_to_o`（＝ 24 個矩陣的跡）、`o3_to_d6`、`o3_to_y`、`o3_to_td`、`invariant_count`、`d_orbital_states`、
  `multipole_coeffs`、`v4_identities`、`sf6_levels`（SF$_6$ 型離心扭曲，含半整數 $J$）、`sf6_clusters`、助記輪 `WHEEL_*`、`wheel_read`、`wheel_o`（讀出的次序 ＝ 能階次序）
- 射線表示：`SpinAlgebra`（帶號乘法表、類算子、冪等元）、`D6S`（(5.7.4) 121 格）、`D2S`、`OS`（八面體）、`half_correlation`、`cn_spin_chars`、`dn_ray_irreps`
- $R_4$ 與 $U_3$：`r4_generators`、`r4_regular`、`r4_ladders`、`coulomb_shell`（$M\otimes1$、$1\otimes N$ 的 $n^2$ 維表示）、`coulomb_energy`、`kepler_check`、`runge_lenz_poisson`、
  `E3`、`su3_cartan`（六邊形的根）、`su3_weights`（半標準 Young 表）、`osc_shell`、殼層閉合 `COULOMB_CLOSURES`、`OSC_CLOSURES`
- 習題 5.1.1–5.6.2：`problem_511` … `problem_561`、`card_values()`（harter_sym 第 5 章場景卡的數字）

`python -m combkit.sym_rot` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 5 章場景卡與導讀上的每個數字（約 5 秒）。

### 4.2y `combkit.symosc`（Python＋numpy）與 `combkit.draw.osc_draw`（ManimGL）

對稱群與耦合振盪子的小演示（symosc 專案，正本 `C:\mypython\manim\projects\symosc\codes\utils\`，用 `tools/sync_from_symosc.py` 複製過來）。
- D3：`d3_matrix`、`K_d3`、`d3_chars_of`、`decompose`、`projector_d3`、`class_sum_values`（類和＝不可約表示的標籤）、`d3_freqs`、`beat`
- 破缺：`P_CYC`、`v_m`、`H_flux`／`flux_levels`（帶相位的耦合）、`gyro_freqs`／`gyro_freqs_numeric`（陀螺項，裂距 √3 g）、`precession`、`K_c2`／`c2_levels`
- Q8：`q8_mul`、`q8_inv`、`Q8_REP2`、`L`／`R`（左右乘法）、`K_q8`、`Q8_CHARS`、`q8_levels`、`q8_irrep_vector`、`q8_E_basis`（標籤 (L(i), R(i))）、`frobenius_schur`
- 畫圖 `osc_draw`：`osc_triangle`、`spring`、`masses_updater`、`radial_bars`、`phasor(s)`、`level_diagram`、`q8_positions`、`q8_graph`

`python -m combkit.symosc` 是自我檢查（暴力積分＝公式、陀螺項的解析式＝一階化的本徵值、左右乘法互換、四元數型的 K 在 E 上是常數）。

### 4.2z `combkit.sym_prod`（Python＋numpy／scipy／sympy）

Harter PSDS 第 6 章（表示的乘積：耦合係數、張量、Wigner–Eckart、光譜選擇規則、Jahn–Teller）用的數學；正本在 harter_sym 的 `codes\utils\sym_prod.py`，用 `tools/sync_from_harter_sym.py` 複製過來；
有限群、八面體 irreps 與 XF$_6$ 模型借用同一個套件的 `sym_nonabel`、`sym_oct`、`sym_rot`（相對 import）。約定：數值的 irrep 組 `IrrepSet`（`oset(basis)`、`ohset(basis)`、`d4set`、`d2set`、`d3set`、`c3set`、`tset`、`qset`、`c3vset`）；
**耦合係數 `cg(S, a, b, gammas)` 用書上 (6.1.18b) 的相位約定**（每一份取 $N_{nkl}\ne0$ 的第一組，$\mathcal C_{kl,n}=\sqrt N>0$），書上的七張表 (6.1.19c)、(6.1.23)、(6.1.26)、(6.1.30)、(6.1.35a,b)、(6.1.36) 逐格相同（`cg_book`、`compare_table`）；
書上的筆誤（(6.1.29)–(6.1.30) 的 $d=-\sqrt2c$、(6.1.36) 的表頭、(6.3.28)ₓ、(6.3.47)、(6.4.8)、(6.4.15b) 的 $1/\ell^\gamma$、(6.5.1)、(6.5.21) 的 $2\Gamma$、(6.7.15) 等）都已更正。
- 兩粒子與外積群：`site_perm`、`outer_op`、`outer_frequencies`、`kron_example`、`two_particle_levels`（＝ 36×36 對角化）、`pair_V`、`t1u_split`（耦合係數 ＝ 直接對角化）
- 耦合係數：`cg`、`cg_norms`、`cg_reduces`、`factorization_lemma`、`unitarity_check`、`solve_6129`（(6.1.27b) 的未知數）、`d2_table`；irrep 三型 `fs_indicator`（Frobenius–Schur）＝ `irrep_type`（不變雙線性型）、`sym_squares` ＝ `cg_symmetry`
- 張量與固體：`unit_tensors`、`tensor_components`、`poly_from_tensors`、`stress_examples`、`elastic_counts`（暴力法）＝ `elastic_counts_chars`、`elastic_tensor`、`moduli` ＝ `moduli_numeric`、
  `spring_model` ＝ `spring_model_bonds`（鈣鈦礦彈簧逐根加總）、`anisotropy`（Kittel 的常數）；$C_{3v}$ 的場 `c3v_field`、`p3_12`、`striction_counts`
- 算子：`six_invariants`、`six_from_levels`、`param_count` ＝ `commutant_count`、`we_check`（六井的偶極矩陣元 ＝ 耦合係數 × 約化矩陣元）、`we_complex_check`（複係數時共軛的位置）
- 光：`lorentz_steady`、`lorentz_transient`、`lorentz_closed` ＝ `lorentz_numeric`（scipy）、`lorentz_trace`、`smith`、`epsilon`、`book_examples`、`rayleigh`、`phasor_points`；
  XF$_6$ 的 `xf6_selection`（特徵標）＝ `xf6_selection_numeric`（偶極、鍵極化率的導數）
- 泛音與 Jahn–Teller：`overtone`（Newton 遞迴 ＝ 單項式基底）、`pascal`、`combination`；`jt_matrix`、`jt_sheets`、`renner_extrema`、`jt_eigvec_track`（繞一圈 −1）、
  `jt_vibronic_levels`（二維振子基底）＝ `jt_radial_levels`（更正的 (6.7.15)）、`jt_effective`
- 習題 6.1–6.6：`problem_61` … `problem_66`（含 $D_{6h}$ 的 `d6h_group`、`d6h_chi`）、`card_values()`（harter_sym 第 6 章場景卡的數字）

`python -m combkit.sym_prod` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 6 章場景卡與導讀上的每個數字（約 6 秒）。

### 4.2z2 `combkit.sym_coup`（Python＋numpy／scipy／sympy）

Harter PSDS 第 7 章（連續轉動群的表示乘積：$R_3$ 耦合係數、張量算子、轉子的 RE 曲面、旋量與偏振、分子軌域）用的數學；正本在 harter_sym 的 `codes\utils\sym_coup.py`，用 `tools/sync_from_harter_sym.py` 複製過來；
$J$ 矩陣、$\mathcal D^j$、球諧與 SF$_6$ 能階借用 `sym_rot`，$O$ 群與誘導表示借用 `sym_oct`（相對 import）。約定：Condon–Shortley 相位；乘積基底 $(m_1,m_2)$ 由大到小；約化矩陣元照書上的 $C^{k\,j_2\,j_1}$ 次序。
書上的筆誤（(7.1.6a)、(7.1.15)、(7.1.17)、(7.1.26)、(7.2.7)、(7.2.9)、(7.2.54)、(7.3.20)、$T_{1u}$／$T_{2u}$ 的標記、(7.3.37b)、Fig. 7.4.2、(7.4.20)、(7.4.24)、(7.5.5e)、(7.5.13)、(7.5.20)、(7.5.23)–(7.5.24)、(7.6.13)、(7.6.32a)、(7.6.40)、(7.6.47)、苯與丁二烯等 79 列）都已更正。
- 耦合係數：`cg`（sympy）＝ `cg_racah`（(7.2.57)）＝ `cg_lowering`（降算子）＝ `cg_float`（大 $J$）＝ `boson_cg`（Schwinger 玻色子 (7.2.54)）；`cg_table`、`stretched_cg`、`threej`、`threej_symmetries`、`sixj` ＝ `sixj_from_3j` ＝ `sixj_from_cg`、`factorization_check`（Euler 角數值積分）
- 原子：`hyperfine_matrix`、`hyperfine_numbers`（21 cm）、`two_electron_terms` ＝ `terms_by_symmetry`（Pauli）、`spin_orbit`、`dirac_level`、`fine_structure`、`p2_levels`（LS→jj）、`cg_curve`、`classical_range`（向量錐）
- 張量：`table_727`、`rr_tensor`、`tensor_T`、`unit_v`、`v_table`（Table 7.1–7.2）、`mixed_table`（Table 7.4）、`elementary_expansion`、`we_v_check`、`we_operator_check`、`crystal_matrix`、`crystal_levels`（O 標記）＝ `crystal_direct`（球諧積分）、`reduced_X`、`r4_hydrogen`、`d2_factor`（$6j$）＝ `d2_factor_numeric`
- 轉子：`re_rigid`、`rigid_levels`（D2 標記）、`fig742`、`rigid_correlation`、`d2_induced`、`octa_re`、`octa_clusters`（$C_4$／$C_3$ 叢與分裂樣式的擬合）、`t46_operator`、`t46_clusters`、`re46`、`sphere_contours`（球面上的 marching squares）、`rigid_paths`、`octa_paths`、`re46_paths`
- 旋量與振子：`spinor_euler`、`spin_vector`、`h_to_omega`、`evolve` ＝ `crank_evolve`、`lissajous`（解析）＝ `lissajous_numeric`（scipy）、`stokes_of`、`stokes_circular`、`stokes_linear`、`ellipse_card`、`beat`、`closure` ＝ `closure_numeric`、`anharmonic_fixed_points`、`osc_re`、`local_mode_levels`
- 分子：`h2p_integrals` ＝ `h2p_quadrature`（球體座標）、`h2p_V` ＝ `h2p_V_generalized`、`two_center`（兩中心 Coulomb 積分的多極展開）、`h2_integrals`、`hl_energy`、`h2_generalized` ＝ `vb_matrix` ＝ `mo_matrix`（更正後同本徵值）、`fig763_levels`、`hybrid_angle`、`xy2_hybrids`、`benzene`、`butadiene`、`octa_sigma`、`ligand_states`
- `card_values()`（harter_sym 第 7 章場景卡的數字）

`python -m combkit.sym_coup` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 7 章場景卡與導讀上的每個數字（約 25 秒）。

### 4.2z3 `combkit.sym_dyn`（Python＋numpy／scipy／sympy）

Harter PSDS 第 8 章（半古典與量子動力學：接觸轉換、作用量與波前、相干態、波包與自相關、半古典輻射、兩能階、量子化的光、Jaynes–Cummings）用的數學；正本在 harter_sym 的 `codes\utils\sym_dyn.py`，
用 `tools/sync_from_harter_sym.py` 複製過來；$J$ 矩陣與 $\mathcal D^j$ 借用 `sym_rot`，自旋向量、曲柄向量、轉子能階借用 `sym_coup`（相對 import）。約定：諧振子 $\hbar=\mu=\omega=1$；兩能階照書上 (8.5.5) 的矩陣、$\Delta=\varepsilon-\Omega$；
JC 的下能階 1 = ↓、上能階 2 = ↑。書上的筆誤（(8.1.27)、(8.1.30b)、(8.2.1)、(8.2.15a)、(8.2.21)、(8.3.3)、(8.3.20)–(8.3.22)、(8.3.32)、(8.4.14)、(8.4.18a)、(8.4.37)、TRK、(8.5.10)、(8.5.21)、(8.5.32)、Fig. 8.5.3 的 $2r$、(8.6.3)、(8.6.20d)、(8.6.31) 的 $1/2\pi$、(8.7.8)、CSRS 的 $\delta/2$ 等 71 列）都已更正。
- 接觸轉換與作用量：`legendre_numeric`、`line_family`、`line_envelope`、`circle_envelope`、`SP_closed` ＝ `SP_initial` ＝ `SP_numeric`、`SH_closed` ＝ `SH_numeric`、`fall_checks`、`trajectory`、`caustic`、`wavefront`、`fig816_curves`、`orthogonality`、`semiclassical_check`、`wkb_ratio`、`huygens_card`
- 相干態：`ladder`、`xp_ops`、`disp_D`、`boost_B`、`Q_op`、`bch_check`（照印的 $+\alpha^*a$ 不么正）、`coherent_cn` ＝ `coherent_numeric`、`coherent_stats`、`ho_eigfun`、`packet_closed` ＝ `packet_sum` ＝ `packet_split_operator`、`SH_ho`
- 波包與譜：`snapshot_sum`、`snapshot_error`（$c_{15}/c_5=0.0299$）、`ac_ho`、`ftac`、`ho_maslov`（$E_{FT}=E_{SC}-2\pi\alpha/\tau$）、`rotor_exact`、`rotor_localized`、`rotor_sc`（自旋相干態 TDVP：重現書上 $E_{SC}=53.134$ 與 (8.3.31) 的係數）、`rotor_ebk`
- 輻射：`gauge_check`（sympy：$A\cdot p\leftrightarrow E\cdot r$）、`hamiltonian_check`、`S_amp`、`spectral_card`、`first_order_check`、`golden_rule_check`、`ho_strengths`、`h_strengths` ＝ `h_box_strengths`（氫 1s 的 TRK）
- 兩能階：`crank`、`h_lab`、`evolve_lab`、`spin_precession`、`rotating_check`、`rabi` ＝ `rabi_numeric`、`floquet_population`、`bloch_siegert`（Floquet：$r^2/4\varepsilon+r^4/64\varepsilon^3$）、`bloch_steady` ＝ `bloch_numeric`、`thermal_sz`、`dressed`、`ac_stark`、`circle_geometry`、`adiabatic`（＝ Landau–Zener）
- 量子化的光與 JC：`photon_elements`、`mode_energy_check`、`angular_sum` ＝ `angular_sum_mc`、`einstein_A` ＝ `einstein_A_gauss`、`h2p_A`（$6.2683\times10^8$ s⁻¹，約化質量 $6.2649\times10^8$）、`golden_A_numeric`、`number_vs_coherent`、`jc_matrix`、`jc_block`、`jc_levels`、`jc_splitting`、`mollow`、`mollow_jc`、`jc_ladder`
- `RECAP`（Harter 全書八章的回顧表）、`card_values()`（harter_sym 第 8 章場景卡的數字）

`python -m combkit.sym_dyn` 是自我檢查：每個結果兩種以上獨立算法、harter_sym 第 8 章場景卡與導讀上的每個數字（約 15 秒）。

### 4.3 `combkit.tikz`（純 Python）

用 Python 組 TikZ：先在原點附近建立圖元，算外框後平移到指定位置（`place`），輸出成 TikZ 指令。
- 基本：`pt`、`emit`、`place`、`node`、`picture`、`colored_math`
- 卡塔蘭物件：`tri_prims`（三角剖分，可依樹上色）、`tree_prims`（二元樹）、`dyck_prims`、`lattice_prims`、`matching_prims`、`paren_math`、`word_math`
- 顏色名稱 `cA`、`cTop`、`cLeft`、`cRight`、`cBad` 要在 LaTeX 導言區用 `\colorlet` 定義（見 `examples/ex03_tikz_figure.py`）。

### 4.3b `combkit.tikz_trees`（純 Python）

平面樹（`counting.binary_to_plane` 的輸出）的 TikZ 畫法：`plane_tree_layout(pt, dx, dy)` → (前序, 座標)，葉子等距、父節點在子節點正中上方、根在 x=0；`plane_tree_prims(pt, dx, dy, labels)`（`labels=True` 標前序編號，給字串如 `"ijk"` 就當頂點名稱），圖元可直接交給 `tikz.place`。（`tikz.py` 由 catalan930 同步產生，所以新的 TikZ 畫法放在這個獨立模組。）

### 4.4 `combkit.draw`（ManimGL）

| 模組 | 內容 |
|---|---|
| `styles` | 顏色、字型、尺寸常數；`zh()` 中文字、`mixed()` 中文與 LaTeX 混排、`boxed()` 加框、`block_color(size)` |
| `scene_base` | `ChuScene`：`show_header(text)` 左上標題、`say(*parts)` 底部字幕（以 `$` 開頭的段是 LaTeX）、`clear_stage()`、`fade_all()` |
| `circle_points` | 圓周等距點（從正上方順時針、0 起算）、弦、點標籤、兩弦交點 |
| `polygon_groups` | 依點編號畫群組（2 點＝弦、3 點以上＝多邊形）、`Configuration`、`make_configuration`、根群組切出的弓形區域 |
| `root_figure` | `RootFigure`：Chu 證明用的「根分解」示意圖；`block_icon` |
| `lattice_paths` | Raney 證明用的路徑：`PathMob`、`word_tiles`、`path_thumbnail`、`level_line` |
| `catalan_draw` | 卡塔蘭物件：`TriangulationMob`、`BinaryTreeMob`、`StepPathMob`（`dyck_path_mob`、`lattice_path_mob`）、`paren_mob`（含括號由內而外長出的動畫）、`matching_item`；統一入口 `catalan_item(family, tree, scale, colored)` |
| `mit_draw` | `ArcDiagramMob`（直線上的弧圖）、`YoungDiagramMob`（英式 Young 圖：`.cell(i,j)`、填數、鉤）、`PlaneTreeMob`（平面樹）、`box_row`、`perm_row` |
| `rep_draw` | Fulton 第 7 章用：`TabloidMob`（列 tabloid）、`ColumnTabloidMob`（帶方向的行 tabloid）、`signed_sum`（帶號線性組合）、`MatrixMob`（`.col(j)`、`.entry(i,j)`）、`perm_tex`、顏色常數 |
| `schubert_draw` | Fulton 第 9–10 章：`CellMatrixMob`（列梯形矩陣，`.stars()`、`.ones()`）、`BoxDiagramMob`（方框裡的 Young 圖）、`HasseMob`（Bruhat 序的 Hasse 圖）、`PipeDreamMob`、`poly_tex`、`word_tex`。注意：`Tex(color=)` 在這版無效，同一個 play 不要同時動畫父物件與子物件 |
| `raney_draw` | 整數序列的 Raney 引理用：`IntPathMob`（每步改變高度 $x_j$ 的路徑，`.vertex(j)`、`.y_of(h)`）、`level_line`、`step_tiles`（數字方塊列）、`readings_row`（各種循環讀法排成一列、基線對齊，標 ✓/✗ 或正點數）、`CardMob`（摸牌遊戲的牌）、`draws_badge`；Raney (1960) 用：`symbol_color`、`word_tiles`（字的方塊列，可切開字列）、`WordTreeMob`（前序的字畫成平面樹，節點依元數上色，可寫標籤）、`forest_mob`、`rank_walk`（秩路徑）、`chu_walk`（倒過來讀的 Chu 路徑）；顏色規則同 Chu 講義（符號 $j$ ↔ $j$ 點群組） |
| `fulton_draw` | Fulton 第 2–4 章用：`tableau_mob`（填好數的 Young 圖）、`SkewMob`（斜表／有洞的表）、`WordMob`（字母方塊，`.apply_move` 播放一步 Knuth 變換）、`word_tex`、`TwoRowArray`（雙列陣列）、`BallMatrixMob`（矩陣球）、版面檢查 `assert_in_frame`、`assert_apart` |
| `young_art` | Young 表的藝術桌布（youngart）：匯入 `art_kit` 與 `young_extras` 的一切，再加 `YPlate`（右上角 PLATE NN · CHAPTER I…X／APPENDIX A、B，改 `PARTS` 可換書）、`cell_center`／`cell_squares`／`tab_squares`（英式方格，逐格上色）、`outline_pts`（分割外框；`[1:-1]` 只剩外側階梯，可當等高線）、`perm_pts`、`rot3`／`proj`（3D 正交投影與深度）、`turn3d(items, view, angle)`（動畫：登記好的 3D 折線每格重新投影，1.26.0）、`shape_of`。注意：`ramp()` 回傳顏色物件，給 `lerp_hex` 前要 `rgb_to_hex(color_to_rgb(c))`；空的點串列不要交給 `dots`／`glow_dots` |
| `formulas` | 常用公式的 Tex：卡塔蘭數、廣義公式、生成函數方程、遞迴式 |
| `walk_draw` | 一維隨機漫步（walk017）：`coin`／`coin_row`（正反硬幣）、`number_line`＋`line_point`、`walker_dot`、`walk_axes`＋`path_point`＋`walk_path`（走法圖）、`fan`（很多條路線）、`sqrt_envelope`（$\pm\sqrt N$）、`bar_chart`＋`bar_of`、`curve`、`plot_axes`（Axes＋Tex 刻度）、`pascal_triangle`＋`pascal_entry`、`delta_grid`＋`grid_cell`、`tex_table`、`erratum_box`（「書上筆誤」紅框）、`contour_circle`、`pole_marker`；元件用原點＋單位長排好，不要事後 `.scale()` |
| `lattice_draw` | 格子上的漫步（walk017 第 2 章）：`lat_point`、`lattice_lines`、`lattice_dots`、`count_grid`＋`grid_number`（格點上的數字）、`neighbor_star`（Fig. 2.3 的四個鄰居）、`trail`（足跡）、`heat_color`／`heatmap`（布里淵區上的色塊）、`heat_cells`（機率擴散）、`bz_frame`、`roots_of_unity`／`root_chain`（單位根頭接尾）、`iso_point`／`lattice3d`（sc／bcc／fcc 的最近鄰斜投影）、`branch_cut`；不 import `walk_draw` |
| `art_kit` | 藝術桌布（walkart）：調色盤 `PAL`（10 組）、`background()`（4K 漸層＋暗角＋顆粒，快取在目前資料夾的 `_artcache/`）、`path`／`gradient_path`／`glow`（三層光暈，保留每條子線的粗細與透明度）、`dots`／`glow_dots`（半徑、顏色、透明度可逐點）、`contours`＋`grid`（contourpy 等高線）、`lattice_walk`／`gauss_walk`／`walk1d`／`pivot_saw`（樞軸法自迴避漫步）、`Plate`（只寫 `art()`；背景、右上角 PLATE 編號、右下角英文標題＋說明或公式自動加；`panel=True` 墊柔邊底色）、`render_png(cls, res, out)`（直接讀相機存 PNG，不用 `-s`；import manimlib 前要把 `sys.argv` 設成 `["manimgl", "-w"]`）、`contact_sheet`；動畫（1.24.0 起）：`animate_plate(Plate 子類別, cfg)` 照原本的 `art()` 建好物件（最後一格＝桌布），0–3 秒依種類「長出來」（線條沿路徑畫出、點雲從 1/4 大小炸開、填色放大淡入）、3–5 秒是 `cfg["B"]` 的動態段，整段鏡頭緩慢推近、文字固定在畫面上；`cfg`：`together`、`A(scene)`、`B(scene)`、`drift`；輔助 `win`（時間窗 rate function）、`travelers`（光點沿折線走，可和畫線同步）、`grow_points`（點雲依 key 逐顆亮起）、`keep`（每格重設點）、`spin`；`render_video(cls, res, fps, outdir, stem)` 直接輸出 mp4；1.26.0 起（youngart 動畫版）：`cfg["prep"](scene)` 在 `art()` 之後執行、`keyed_reveal(m, keys)`（依 key 一格格出現，一個 updater，2500 格也快）、`mirror(items, axis)`（沿平面上的軸 3D 翻轉再翻回＝鏡射）、`twinkle`（點雲逐點閃爍）、`pulse`（點雲逐點呼吸）。注意：這版 wgpu 渲染器的 DotCloud 緊接在 VMobject 後面畫會出錯，`Plate.add()` 會自動在中間夾透明小圖 |
| `sym_draw` | Harter PSDS（harter_sym，第 1–8 章共用）：`mat_tex`（圓括號矩陣，格子可依 (i,j)／列／行上色；`mat_entry`、`mat_row`、`mat_col`、`labeled_mat`）、`analyzer`（Stern–Gerlach 分析器盒，盒內箭頭是分析器軸；`recombine=True` 是合併盒）＋`port`、`beam`（線寬隨強度）、`split_beams`、`vec_plane`、`vec_arrow`、`eigen_line`、`circle_image`（單位圓 → 橢圓的變形）、`contour_around`（圍道）、`SpringChain`（牆—彈簧—質量鏈，`set_displacements` 給 updater 用）；第 2 章（1.28.0）加 `PendulumRow`（一排單擺＋彈簧）、`PendulumRing`（俯視的單擺環，切向／徑向位移，彈簧可逐根上色、改粗細）、`PhasorClock`（相量鐘面，`set_value(z)`）、`spring`／`spring_points`（任意方向的彈簧）、`band_dots`、`band_curve`、`level_lines`（能帶、色散、能階）、`quad_contours`（二次位能的等位橢圓）；第 3 章（1.29.0）加 `mirror_lines`、`BobTriangle`（固定在紙上的鏡子與三個擺錘）、`SpringTriangle`（六自由度彈簧質量模型，`set_displacements`）、`sym_table`（群表／特徵標表）、`perm_grid`（正則表示的格子）、`View3D`＋`SphereView`（投影的 Hamilton 球與半角弧）、`NH3Model`（投影的球棍氨分子）、`polar_curve`（環上的波）；第 4 章（1.30.0）加 `CubeView`（投影的立方體＋上色的體對角線＋內接八面體，`set_rotation`）、`OctaSites`（八面體六頂點，振幅的大小與相位上色）、`level_diagram`（能階分裂／相關圖）、`XF6Model`（投影的球棍 XF6，`set_displacements`、`arrows`）、`rot_matrix`、`irrep_color`；第 5 章（1.31.0）加 `TriadBall`（投影的球＋跟著轉的赤道／子午線＋本體三腳架＋實驗室軸，`set_rotation`（可非正常）、`add_crank`、`trail`）、`VectorCone`（角動量向量錐，`set_beta` 傾斜）、`cone_family`、`ProbBars`（機率／有號振幅長條，`set_values`）、`ladder_column`（能階梯，簡併度寫成文字）、`weight_diagram`（根圖／權圖，重數畫成同心圈）、`mnemonic_wheel`＋`wheel_marker`（助記輪）；第 6 章（1.33.0）加 `PairOcta`（八面體上的兩個粒子，`set_pair`）、`PairGrid`（6×6 的兩粒子格子）、`StrainCube`（投影立方體的應變 `set_strain` 與面上的應力箭頭 `face_arrows`）、`field_glyph`（⊙ ⊗ ±）、`current_loop`、`trace_curve`（時間曲線）、`SmithChart`（Γ 圓、Δ 圓、相量）、`hat_contours`（Jahn–Teller 墨西哥帽的等高線）、`hat_profile`（截面）；第 7 章（1.34.0）加 `SpinTop`、`SpinPair`（兩個自旋小陀螺）、`vector_triangle`（角動量相加的三角形）、`CouplingCone`（$j_1$ 繞 $J$ 的錐與古典範圍，`set_phase`）、`RESurface`（RE 曲面的投影等高線＋淡色網，`set_view`；`mesh_lines`）、`PoincareSphere`（$A,B,C$ 三軸的 Stokes 球，`add_stokes`、`set_stokes`）、`polarization_ellipse`、`ellipse_from_stokes`、`projected_curve`（3D 折線，振子環面）、`box_frame`、`potential_curves`、`min_marker`；第 8 章（1.35.0）加 `polyline`、`line_family`（直線族）、`trajectory_fronts`（軌跡族＋等作用量波前）、`PhasePhasor`（相空間相量＋Gauss 圓斑，`set_alpha`）、`snapshot_stack`（波包快照列）、`CrankSphere`（曲柄與自旋向量，`add_vector`、`set_vector`、`add_trail`、`cone_rim`）、`dressed_ladder`（JC 穿衣能階梯）——3D 一律投影成 2D（這台 wgpu 渲染器的 Surface／DotCloud 緊接 VMobject 畫會出錯）；也匯出 walk_draw 的 `tex_table`、`erratum_box`、`plot_axes`、`curve`、`pole_marker` |

---

## 5. 範例（`examples/`）

| 檔案 | 內容 | 執行（在 `examples\` 裡） |
|---|---|---|
| `ex01_dyck_paths.py` | 算出 n=3 的 5 條 Dyck 路徑並和暴力列舉比對，再上色畫出 | `manimgl ex01_dyck_paths.py DyckPathsExample` |
| `ex02_row_insertion.py` | Fulton p.7 的列插入：算出新表與擠壓路徑，畫出前後兩個 Young 表並標出路徑 | `manimgl ex02_row_insertion.py RowInsertionExample` |
| `ex03_tikz_figure.py` | 不用 manim：產生五邊形三角剖分＋Dyck 路徑的 TikZ 圖 | `python ex03_tikz_figure.py`，再 `xelatex ex03_tikz_figure.tex` |
| `ex04_rsk.py` | R-S-K：算出排列 4 1 3 5 2 的 (P, Q) 並確認能還原，再逐步畫出 P 長大、Q 記錄步數 | `manimgl ex04_rsk.py RSKExample` |
| `ex06_specht.py` | Specht 模組：$\lambda=(2,1)$ 的 $v_T=\{T\}-\{(13)T\}$ 與 $\rho(12)$ 的矩陣 | `manimgl ex06_specht.py SpechtExample` |
| `ex07_schur_module.py` | Schur 模（Fulton 第 8 章）：$\lambda=(2,1)$、$m=3$ 的 8 個半標準表按權重分 7 堆，合起來是 $s_{(2,1)}(x_1,x_2,x_3)$ | `manimgl ex07_schur_module.py SchurModuleExample` |
| `ex08_raney_cycle.py` | Raney 的循環引理：$\sigma=0020030$（秩 $-2$）的 7 種循環讀法中恰有 2 種是 2 字列，各畫成 Chu 講義的路徑；雙重計數得 $L(2;0,1,1)=12$ | `manimgl ex08_raney_cycle.py RaneyCycleExample` |
| `ex11_random_trees.py` | 用循環引理產生均勻隨機的樹：組成 $(1;0,2,1)$ 的 21 個字中抽 4 次，畫出打亂的字串 → 字 → 樹 → 運算式 | `manimgl ex11_random_trees.py RandomTreesExample` |
| `ex12_schubert.py` | Schubert 胞腔與 pipe dream（Fulton 第 9–10 章）：$w=21543$ 的列梯形矩陣（星號數 = 逆序數）、$S_3$ 的 Bruhat 序、$\mathfrak S_{21543}$ 裡係數 2 的兩張 pipe dream | `manimgl ex12_schubert.py SchubertExample` |
| `ex13_evacuation_keys.py` | evacuation 與鑰匙（Fulton 附錄 A）：書上 p.184 的 $T\mapsto T^*$、p.210 的 $K_-(T)$、$T$、$K_+(T)$ | `manimgl ex13_evacuation_keys.py EvacuationKeysExample` |
| `ex10_raney_lemma.py` | 整數序列的 Raney 引理：$x=(2,-1,-3,1,2)$ 的 5 種循環讀法中只有從最後一個最低點開始的是好讀法；加強版的正點數 $3,1,2,5,4$ | `manimgl ex10_raney_lemma.py RaneyLemmaExample` |
| `ex09_wick_moments.py` | 不用 manim：6 點的 15 種配對裡 5 種不交叉；GOE 的 $\frac1N\mathbb{E}\,\mathrm{Tr}\, J^4=2+5/N+5/N^2$；$k=3$ 的樹形閉路徑 ↔ Dyck ↔ 平面樹 ↔ 配對，畫成 TikZ | `python ex09_wick_moments.py`，再 `xelatex ex09_wick_moments.tex` |
| `ex05_lr_rule.py` | LR 規則：(3,2,1)/(2,1) 上內容 (2,1) 的 3 個斜表裡有 2 個 LR 斜表 = $s_{(2,1)}^2$ 裡 $s_{(3,2,1)}$ 的係數 | `manimgl ex05_lr_rule.py LRRuleExample` |
| `ex14_random_walk.py` | 一維隨機漫步（Rudnick & Gaspari 第 1 章）：十次擲幣 HTHHTHTTHH 的走法圖、40 條路線與 $\pm\sqrt N$ 包絡、$N=10$ 的二項分布（$C(10,0)=252$）與高斯曲線 | `manimgl ex14_random_walk.py RandomWalkExample` |
| `ex15_lattice_gf.py` | 格子上的漫步（Rudnick & Gaspari 第 2 章）：正方格子 $N=4$ 的走法數格子（回到原點 36＝遞迴＝兩枚硬幣＝布里淵區積分）、$(\chi/4)^4$ 的熱圖、sc／bcc／fcc 的最近鄰 | `manimgl ex15_lattice_gf.py LatticeGFExample` |
| `ex16_recurrence.py` | 回家與維度（Rudnick & Gaspari 第 3 章）：一維從不回家的走法 1…252（生成函數＝暴力列舉）、二維 400 步的足跡（159 個點）、一二三維 1000 步還沒回家的機率 0.0252／0.3218／0.6685、sc／bcc／fcc 的逃脫機率（Watson 積分） | `manimgl ex16_recurrence.py RecurrenceExample` |
| `ex17_image_walker.py` | 鏡像走路人與吸收牆（Rudnick & Gaspari 第 4 章）：碰過牆的路線 $1\to0\to1\to2\to3$ 與它的鏡像路線、$x_0=3$ 的存活機率 0.6562／0.2356／0.0756 vs $x_0\sqrt{2/\pi N}$、一維穩定態的帳篷 1.2、2.4、1.6、0.8 與流進兩牆的 3/5、2/5 | `manimgl ex17_image_walker.py ImageWalkerExample` |
| `ex18_persistent_walk.py` | 有慣性的漫步與電報方程（Rudnick & Gaspari 第 5 章）：$Q=\frac34$ 與 $\frac14$ 走 20 步的精確分布（變異數 56.00／7.11 vs $NQ/(1-Q)$ 60／6.67）與高斯、持續長度 3 與 1/3、電報方程 $t=1,3,8$ 的分布與尖峰重量 0.1839／0.02489／0.0001677 | `manimgl ex18_persistent_walk.py PersistentWalkExample` |
| `ex19_walk_shape.py` | 隨機漫步的形狀（Rudnick & Gaspari 第 6 章）：「東東北西」的 $T=\frac1{25}\begin{pmatrix}14&3\\3&6\end{pmatrix}$（定義式＝Kramers），二維 1000 步的重心與主軸（$R_1^2=258.8$、$R_2^2=31.6$，旋轉不變）、三維 2 萬條 100 步的個別非球度直方圖、$A_3=\frac{10}{19}$ vs 模擬 0.5231、$A'_3=0.394274$ | `manimgl ex19_walk_shape.py WalkShapeExample` |
| `ex20_saw_flory.py` | 自迴避漫步與 Flory 指數（Rudnick & Gaspari 第 7 章）：二維精確列舉 $c_N$（DFS＝文獻＝暴力）、$\langle R^2\rangle$ 對 $N$ 的點 vs 普通漫步的直線 $N$ 與 $N^{3/2}$（$\langle R^2\rangle_{14}=42.79$、$2\nu_{\rm eff}=1.4556$）、固定種子的 10 步自迴避／普通漫步足跡、$1-f$ 部分和 $256\to28\to124\to96\to100=c_4$、Flory $\nu=\frac3{d+2}$ | `manimgl ex20_saw_flory.py SawFloryExample` |
| `ex21_scaling_collapse.py` | 標度入門（Rudnick & Gaspari 第 8 章）：一維格子 $g(z;q)$ 在 $z_c-z=0.1,0.05,0.01$ 的三條曲線（$g(z;0)=5,10,50$）塌縮成 $(z_c-z)g=\frac2{4+x^2}$、長尾巴步伐的寬度 $\propto N$（1822 vs 28.06）、最陡下降更正版 $e^{-x^2/42}$ vs 書上 $e^{-x^2/21}$（$N=20$、$x=8$：0.2098 精確）、上臨界維度 $N^{(4-d)/2}$ | `manimgl ex21_scaling_collapse.py ScalingCollapseExample` |
| `ex22_ising_scaling.py` | 臨界現象與 $n\to0$（Rudnick & Gaspari 第 9 章）：平均場 vs 二維 Ising（Onsager）的 $M_s(T)$（log–log 斜率 $\beta=0.5000$、0.1250）、Landau 雙井 $r=1,0,-1$（谷底 $m^*=\pm1$）、一維 Ising $\chi=e^{2K}=7.3891$、三組指數的標度律、Metropolis 翻轉等待 133、608、2151 sweep、$2\times3$ 梯子上 $n\to0$ 只剩自迴避路徑（$K+K^3+K^5=0.65625$） | `manimgl ex22_ising_scaling.py IsingScalingExample` |
| `ex23_torus_saw.py` | 關在盒子裡的自迴避漫步（Rudnick & Gaspari 第 10 章）：4×4 環面的 $c_N$（$N<4$ 與無限平面相同，$N=4$ 少 4 條，走滿 16 格 5712 條，$16\times5712=91392$＝四維 Gray 碼）、一條走滿的漫步（越界的一步拆成兩個半段）、3×3、4×4、5×5 的步數密度 $\rho(z)=\langle N\rangle/V$（漲落峰 $z^*=0.4443$、0.4380、0.4272 往 $z_c$ 移） | `manimgl ex23_torus_saw.py TorusSawExample` |
| `ex24_decimation_flow.py` | 重整化的流（Rudnick & Gaspari 第 11 章）：一維 Ising 環的 $Z$（$N=4,8$；暴力＝轉移矩陣＝抽格連乘 121.2329、9160.4387）、$K_n$：1、0.6625、0.3501、0.1137 → 0，$x=e^{-2K}$ 的蛛網 $x'=\frac{2x}{1+x^2}$ 流向 1；一維 Gaussian 鏈 $r''=r(r+4)$ 的蛛網 0.05 → 0.2025 → 0.8510 → 4.1282（小 $r$ 時每次 ×4 ⇒ $\Lambda_r=2$、$\nu=\frac12$） | `manimgl ex24_decimation_flow.py DecimationFlowExample` |
| `ex25_wilson_fisher.py` | Wilson–Fisher 不動點（Rudnick & Gaspari 第 12 章）：一迴圈遞迴 $\frac{dr}{d\ell}=2r+\frac{12u}{1+r}$、$\frac{du}{d\ell}=\epsilon u-\frac{36u^2}{(1+r)^2}$ 在 $\epsilon=1$ 的 $(r,u)$ 流線（臨界線用二分法求出、流進 $(-\frac17,\frac1{49})$，特徵值 1.8968、−1.2301）、$\nu$ 的一階公式 $n=0..3$：0.5625、0.5833、0.6000、0.6136 對已知值（引用）、5 個片段的配對 $\frac{5!}{3!2!}=10$（列舉 945 種） | `manimgl ex25_wilson_fisher.py WilsonFisherExample` |
| `ex26_art_wallpaper.py` | 藝術桌布（walkart 的做法）：`Plate` 子類別只寫 `art()`，`render_png` 直接存 PNG；畫 400 條一維漫步的扇形（±√N 包絡）與第 4 代 Koch 雪花（D = ln4/ln3 ≈ 1.2619），再拼成接觸表 | `python ex26_art_wallpaper.py` |
| `ex27_art_animation.py` | 桌布變 5 秒動畫（walkart 動畫版的做法）：`animate_plate` 讓 400 條漫步一起畫出、6 顆領頭光點同步跑、鏡頭推近、文字最後淡入；`render_video` 輸出 640×360、15 fps 的 mp4 | `python ex27_art_animation.py` |
| `ex28_young_wallpaper.py` | Young 表的桌布（youngart 的做法）：`YPlate` 的角落標籤；(3,3,3) 的 42 個標準表畫成 3×3 方框 Young 格裡兩端收攏的 42 條路（$f^{(3,3,3)}=42=\deg Gr(3,6)$）、900 個數的隨機排列的 Q 表依長出來的時間上色＋等高線（`rsk_rows`，第一列長＝最長遞增子序列）；加 `--anim` 另出 5 秒動畫（`animate_plate`＋`prep`＋`keyed_reveal`：格子依長出來的時間一格格出現） | `python ex28_young_wallpaper.py`（或加 `--anim`） |
| `ex29_eigen_projectors.py` | 本徵投影算符（Harter PSDS 第 1 章 §1.2）：$M=\begin{pmatrix}4&1\\3&2\end{pmatrix}$ 的 $\mathcal P_1=\frac14\begin{pmatrix}1&-1\\-3&3\end{pmatrix}$、$\mathcal P_5=\frac14\begin{pmatrix}3&1\\3&1\end{pmatrix}$（冪等、正交、完備；$M^{50}=\mathcal P_1+5^{50}\mathcal P_5$＝直接乘；預解式的圍道積分重現 $\mathcal P_j$），畫出投影矩陣（行＝ket、列＝bra 上色）、單位圓→橢圓與本徵方向、$v=\mathcal P_1v+\mathcal P_5v$、60° 分析器 | `manimgl ex29_eigen_projectors.py EigenProjectorsExample` |
| `ex30_cn_ring.py` | C_n 對稱的單擺環（Harter PSDS 第 2 章 §2.6–2.7）：$C_3$ 投影算符「最小方程連乘」＝「特徵標」＝ `sym_mat.projectors`；6 擺環 $\omega^2=0,1,3,4,3,1$（＝數值）；24 擺環每隔一條彈簧變弱（$\bar a=\frac32$、$\underline a=\frac12$）裂成聲學帶 $[0,1]$ 與光學帶 $[3,4]$（公式＝24×24 數值），畫出 $k=1$ 行波（六個相量鐘依序差 60°）、摺疊的 $C_{24}$ 頻譜 $2\mp2\cos\frac k2$ 與兩條帶、$C_3$ 特徵標表 | `manimgl ex30_cn_ring.py CnRingExample` |
| `ex31_character_table.py` | 由生成元算特徵標表（Harter PSDS 第 3 章 §3.2、§3.5）：只給兩個 2×2 生成元，`FiniteGroup.from_generators` 封閉出 $C_{3v}$（6 個）與 $D_4$（8 個）；共軛類、類代數、最小方程（$c^3-9c$；$D_4$ 要兩個 $c^3-4c$）、全對易冪等元 ⇒ $\ell=({}^\circ G\,p_1)^{1/2}$、$\chi^*=p\,{}^\circ G/\ell$；＝不可約表示的跡＝書上的 (3.5.8)、(3.6.3)（差一個重新命名），列正交、$\sum\ell^2={}^\circ G$；畫出兩張表 | `manimgl ex31_character_table.py CharacterTableExample` |
| `ex32_octahedral_group.py` | 由轉動矩陣產生八面體群（Harter PSDS 第 4 章 §4.1–§4.2）：兩個 3×3 生成元（繞 $x$、$y$ 各 90°）封閉出 24 個轉動；類大小 1, 8, 3, 6, 6、$c_i$ 的最小方程 $c^5-40c^3+144c$、類代數的特徵標表＝書上 (4.1.11)＝`sym_oct.O_CHARS`、四條體對角線的置換 ⇒ $O\cong S_4$；`sym_oct.correlation("D4")` ＝ (4.2.42a)，並用 Frobenius（誘導表示的跡）核對每一行；畫兩張表與轉動中的立方體（`SD.CubeView`） | `manimgl ex32_octahedral_group.py OctahedralGroupExample` |
| `ex33_wigner_d.py` | Wigner 的轉動矩陣（Harter PSDS 第 5 章 §5.4D–§5.5B）：$d^j_{mn}(\beta)$ 用 Wigner 的階乘公式與 $e^{-i\beta J_y}$ 各算一次（$j=\frac12\ldots4$ 一致、$j=20$ 誤差 $<10^{-9}$），第三種算法是 Schwinger 的兩個振子；$\mathcal D^{1/2}(360^\circ)=-\mathbf 1$；書上 (5.4.35) 照印的 $\mathcal D^2$ 在 $\beta=0.7$ 時 $(\mathcal D\mathcal D^T)_{11}=1.3174$（不正交）；Stern–Gerlach：轉過 $\beta$ 的 $\vert2,1\rangle$ 五個通道的機率 $\vert\mathcal D^2_{m'1}\vert^2$（$\beta=90^\circ$：0.25, 0.25, 0, 0.25, 0.25；$\langle J_z\rangle=\cos\beta$）；畫隨 $\beta$ 傾斜的向量錐（`SD.VectorCone`）與機率長條（`SD.ProbBars`） | `manimgl ex33_wigner_d.py WignerDExample` |

| `ex34_symmetric_oscillators.py` | 對稱群與耦合振盪子（symosc）：三個擺的 D3 約化 $A_1\oplus E$ 與 $E$ 模式、$C_3$ 的相量 $(1,\omega,\omega^2)$、Q8 的 Cayley 圖與 $E$ 的相量（$E$ 出現兩次、四重簡併） | `manimgl ex34_symmetric_oscillators.py SymmetricOscillatorsExample` |
| `ex35_clebsch_gordan.py` | 投影算子算 Clebsch–Gordan 係數（Harter PSDS 第 6 章 §6.1B）：$O_h$ 的特徵標給 $T_{1u}\otimes T_{1u}=A_{1g}\oplus E_g\oplus T_{1g}\oplus T_{2g}$；投影公式 (6.1.18b)（書上的相位約定）算出的 9 × 9 表與書上 (6.1.19c) 逐格相同，$\mathcal C^\dagger(\mathcal D\otimes\mathcal D)(R_1)\mathcal C$ 是 1、2、3、3 的方塊（48 個元素都對）；$T_{1g}$ 的 $N_{123}=\frac12$；畫直積矩陣與方塊對角的格子（`SD.perm_grid`）與表的前三行（`SD.sym_table`） | `manimgl ex35_clebsch_gordan.py ClebschGordanExample` |
| `ex36_re_surface.py` | 不對稱陀螺的轉動能量（RE）曲面與量子化路徑（Harter PSDS 第 7 章 §7.4A）：$J=10$、$A,B,C=0.2,0.4,0.6$ 的 21 個能階（對角化 ＝ 三對角公式），分隔線 $BJ(J+1)=44$，上下兩邊的能階成對（$C_2(x)$ 叢 $(A_1B_1)$、$C_2(z)$ 叢 $(B_2A_1)$），最低一對只差 $8.6\times10^{-7}$ cm⁻¹；在每個本徵能量畫 RE 曲面的等高線（`SC.rigid_paths`），用 `SD.RESurface` 投影、鏡頭繞轉 | `manimgl ex36_re_surface.py RESurfaceExample` |
| `ex37_coherent_packet.py` | 相干態波包在諧振子裡來回不散（Harter PSDS 第 8 章 §8.2–8.3）：$|\alpha|^2=5$ 的 Poisson 振幅（$Q|0\rangle$ 的截斷矩陣 ＝ 閉式，$c_4=c_5=0.4189$、$\Delta n=\sqrt5$），波包閉式 (8.2.28) ＝ 本徵函數和 ＝ FFT 分裂算子；`SD.PhasePhasor` 順時針轉一圈、波包跟著來回；十張快照乘 $e^{i(5+1/2)t_p}$ 相加 ＝ $\langle x|5\rangle$（誤差 $c_{15}/c_5=0.0299$） | `manimgl ex37_coherent_packet.py CoherentPacketExample` |

寫新場景時照範例的結構：檔頭加 `sys.path`（或已安裝）→ 從數學層取資料並 `assert` → 用 `combkit.draw` 畫 → 繼承 `ChuScene`。

---

## 6. 工具與測試

| 指令 | 用途 |
|---|---|
| `python tests/run_selftests.py` | 純 Python 部分的自我測試；加 `--draw` 再檢查 `combkit.draw` 每個模組都能載入 |
| `python tools/render_check.py <場景檔> <類別>` | 低畫質試渲染＋總覽圖（`runs/<類別>/sheet.png`）；`--frames 3,8` 另抽原尺寸畫面、`--hd` 輸出 1080p。每次在獨立資料夾、用獨立的公式快取與嚴格 LaTeX 模式，可以同時跑好幾個 |
| `python tools/make_api.py` | 重新產生 `API.md` |
| `python tools/sync_from_walk017.py` | walk017 的 `codes/utils/walk1d.py`、`walk_draw.py`、`walk_gf.py`、`lattice_draw.py` 改過之後執行：複製成 `combkit/walk1d.py`、`combkit/walk_gf.py`、`combkit/draw/walk_draw.py`、`combkit/draw/lattice_draw.py`（加一行來源註記），再更新 `API.md` |
| `python tools/sync_from_harter_sym.py` | harter_sym 的 `codes/utils/sym_mat.py`、`sym_abel.py`、`sym_nonabel.py`、`sym_oct.py`、`sym_rot.py`、`sym_prod.py`、`sym_coup.py`、`sym_dyn.py`、`sym_draw.py` 改過之後執行：複製成 `combkit/sym_mat.py`、`combkit/sym_abel.py`、`combkit/sym_nonabel.py`、`combkit/sym_oct.py`、`combkit/sym_rot.py`、`combkit/sym_prod.py`、`combkit/sym_coup.py`、`combkit/sym_dyn.py`、`combkit/draw/sym_draw.py`（加來源註記），再重做 `API.md` |
| `python tools/sync_from_catalan930.py` | catalan930 的 `codes/utils` 改過之後執行：只覆寫來自 catalan930 的那幾個檔案（它們請改 catalan930 那份再同步，不要在這裡直接改），其他計畫加進來的模組與 `__init__.py` 不會被動到；最後更新 `API.md` |
| `python tools/make_zip.py` | 先跑自我測試，通過後打包成 `C:\mypython\manim\projects\combkit_v<版本>.zip`；同版本的 zip 已存在時不覆寫（要覆寫加 `--force`），所以舊版 zip 都會留著 |
| `python tools/make_bundle.py` | 連同 ManimGL 一起打包（給新機器用）：`..\combkit_v<版本>_with_manimgl.zip`，內含 `manimgl\`（`C:\mypython\manim` 裡 git 追蹤的原始碼，commit fafa083）、`combkit\`、`docs\`（環境與交接、使用說明、新計畫起步指南）、`先讀我_README.md`、`install.bat`（裝固定版本的套件、`pip install --no-deps -e` 兩個套件、跑自我測試）、`requirements_pinned.txt`；先跑自我測試並確認 ManimGL 工作樹沒有修改，同名 zip 不覆寫 |

`runs/` 是試渲染的工作資料夾，可以隨時刪掉。

---

## 7. 在 catalan930 裡，所有程式怎麼分工

catalan930（`C:\mypython\manim\projects\catalan930`）用的是同一批程式（放在 `codes/utils/`，import 名稱是 `utils.X`；
本工具包是它的複本，import 名稱改成 `combkit.X`）。

| catalan930 的程式 | 做什麼 | 用到的模組 |
|---|---|---|
| `codes/scenes/scene00a…00m`（卡塔蘭基礎模組） | 三角剖分、括號、二元樹、格子／Dyck 路徑、投票、配對、生成函數 | counting 第 5 節、catalan_draw、polygon_groups |
| `codes/scenes/scene01…09`（Chu 1987） | 不交叉群組、根分解、遞迴、生成函數、Raney 證明 | counting 第 1–4 節、circle_points、polygon_groups、root_figure、lattice_paths、formulas |
| `codes/scenes/sceneM01…M14`（MIT 18.212） | Cauchy 乘積、醉漢走路、平面樹、弧圖、排序、模式、Young 表、鉤長公式 | counting 第 6 節、mit_draw、catalan_draw、lattice_paths |
| `codes/scenes/sceneF01…F08`（Fulton 第 1 章） | Young 表、Schur 多項式、列插入、乘積、jeu de taquin | tableaux、mit_draw |
| `codes/scenes/sceneF09…F19`（Fulton 第 2–4 章） | Knuth 等價、主定理、Pieri、Kostka、行字、L(w,k)、唯一性、R-S-K、矩陣球、應用 | tableaux 第 7 節、fulton_draw、mit_draw |
| `codes/scenes/sceneF28…F41`（Fulton 第 7 章） | 對稱群的表示：tabloid、Specht 模、標準基底、矩陣與特徵標、Young's rule、誘導積與 $\Lambda\cong R$、Frobenius 公式、分支律、straightening、Young 對稱化子、模 3 | rep_data、rep_draw、fulton_draw、tableaux |
| `codes/scenes/sceneF20…F27`（Fulton 第 5–6 章） | 表的因數分解與 LR 數、對稱性、反格字、LR 規則、Remmel–Whitney、對稱多項式的基底、Jacobi–Trudi、對稱函數環 | tableaux 第 8 節、fulton_draw |
| `codes/scenes/sceneF42…F45`（Fulton 第 8 章） | Schur 模 $E^\lambda$、權重與 $\operatorname{Char}=s_\lambda$、Schur–Weyl $E^{\otimes n}\cong\bigoplus E^\lambda\otimes S^\lambda$、Plücker 二次關係（小工具寫在場景檔裡） | tableaux、fulton_draw |
| `codes/scenes/sceneF46…F55`（Fulton 第 9–10 章） | Plücker 座標、旗流形的座標環、Grassmannian 的胞腔、四條線問題（Schubert 計算）、旗流形的胞腔、差商、Schubert 多項式、pipe dream、Bruhat 序、Grassmann 置換 = Schur 多項式 | schubert、schubert_draw、tableaux、fulton_draw |
| `codes/scenes/sceneF56…F65`（Fulton 附錄 A–B 與全書回顧） | evacuation、行插入、Burge、四個角、0-1 矩陣與 Gale–Ryser、鑰匙；胞腔與同調、子簇的類與次數、Chern 類；全書地圖 | fulton_appendix、schubert、tableaux、fulton_draw |
| `codes/scenes/sceneR01…R06`（Raney 引理筆記） | 整數序列的 Raney 引理、最後一個最低點的證明、加強版、出棧序列與 $C_n$、摸牌遊戲 ＝ Raney 的字、勘誤小結 | raney 第 9 節、raney_draw、fulton_draw（版面檢查） |
| `codes/scenes/sceneN01…N12`（Raney (1960) 導讀） | 前序記法與字、還欠幾個引數、三種讀法、Raney 的循環引理證明、雙重計數與 Chu Fig. 1、隨機例子產生器、切開與第一個符號、L*、Lagrange 反演、三項式與 Problem 4983、卷積、總結 | raney、raney_examples、raney_draw、polygon_groups |
| `codes/main.py` | 每個場景的同名子類別＋十四種連續播放（含 `FultonDemo2`、`FultonDemo3`、F28–F41 的 `FultonDemo4`、F42–F45 的 `FultonDemo5`、F46–F55 的 `FultonDemo6`、F56–F65 的 `FultonDemo7`、R01–R06 的 `RaneyLemmaDemo`、N01–N12 的 `RaneyPaperDemo`）（`play_sequence` 依序呼叫各場景的 `construct`） | — |
| `codes/render_all.py` | 一次輸出很多支影片（`--intro`、`--chu`、`--mit`、`--fulton`、`--fulton2`、`--fulton3`、`--fulton4`、各種 demo） | — |
| `codes/check_env.py` | 環境檢查（Python、ManimGL、LaTeX、字型、兩個數學模組的自我測試、試渲染） | counting、tableaux |
| `notes/raney/make_figs.py`（Raney 導讀） | 先檢查導讀寫到的每一個數字，再產生 TikZ 圖與表（沒有對應的動畫） | raney、counting |
| `notes/make_catalan_figs.py`、`notes/make_mit_figs.py`、`notes/mit/make_figs_partK.py` | 講義第 1 節與附錄 A 的 TikZ 圖；先檢查講義寫到的數字再產生 | counting、tikz 的那些函式（本工具包的 `tikz.py` 就是從這裡抽出來的） |

各場景的內容與長度見 catalan930 的 `codes/README.md` 與 `交接說明_HANDOFF.md`。

---

## 8. 注意事項（在 catalan930 踩過的坑）

1. 這是 **ManimGL**（`manimgl`、`ShowCreation`、`from manimlib import *`），不是 Manim Community（`manim`、`Create`）。
2. `custom_config.yml` 只能寫英文（繁中 Windows 用 cp950 讀它）。
3. 中文不要放進 `Tex`，用 `styles.zh()`／`mixed()`；`t2c` 上色要整個因子一起（例如 `x_1^{n_1}`）。
4. Tex 字串開頭／結尾不要有 `\;`、`\,`；用 `t2c` 時不要用 `\left( … \right)`（改 `\Big(`）。否則「帶顏色標記」的編譯會失敗，
   ManimGL 卻不報錯、沿用上一個公式。`tools/render_check.py` 的嚴格模式會把這種錯誤直接報出來。
5. `Tex("a", "b")[1]` 是第 2 個字形，不是第 2 段。
6. 要連續播放多個場景（像 catalan930 的 `play_sequence`）時，輔助函式寫成**模組層級的函式**，不要寫成場景的方法。
7. 不要在同一個資料夾同時跑多個 `manimgl`（`latex_cache/working.tex` 會互相覆蓋）；ManimGL 預設的公式快取
   （`%LOCALAPPDATA%\manim\manim\Cache`，SQLite）在很多程式同時寫入時曾經損毀（錯誤訊息 `database disk image is malformed`）。
   平行渲染請用 `tools/render_check.py`（各自的資料夾與快取）。快取損毀時把那個資料夾刪掉或改名即可，ManimGL 會重建。
8. 自訂 Mobject 的屬性不要和 ManimGL 撞名（`.parents`、`.color`、`.stroke_width`）；放在 dict／list 裡的子物件不會跟著 `.copy()` 複製。

---

## 9. 版本、來源與授權

- 版本號只寫在 `combkit/__init__.py` 的 `__version__`（`pyproject.toml` 會讀它）。
- 1.0.0（2026-10-06）：由 catalan930 的 `codes/utils/`（2026-10-05 版）與 `notes/make_catalan_figs.py` 產生。
- 1.1.0（2026-10-10）：加入 Fulton 第 2–4 章（`tableaux` 第 7 節、`draw/fulton_draw.py`、`examples/ex04_rsk.py`），附上環境說明。
- 1.5.0（2026-10-14）：加入 `raney`（Raney 1960 論文導讀用的數學）與 `examples/ex08_raney_cycle.py`；環境說明的副本更新到 2026-10-13 版；
  `tools/make_api.py` 改成自動找出所有模組（之前漏列 `rep_data`、`fulton_draw`、`rep_draw`）。
- 1.17.0（2026-10-18）：加入 walk017 第 7 章（Path integrals and self-avoidance，動畫 W701–W706）的 `walk_saw`（路徑和的切段、二維自迴避漫步的精確列舉、1−f 展開與弧線圖、Dyson 方程、冪次數算、平均場、Flory）與 `examples/ex20_saw_flory.py`；`tools/sync_from_walk017.py` 多同步 `walk_saw.py`，`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.18.0（2026-10-18）：加入 walk017 第 8 章（Properties of the random walk: introduction to scaling，動畫 W801–W804）的 `walk_scale`（長尾巴步伐與非高斯、生成函數的標度塌縮、最陡下降與 Fisher 指數、γ = (2−η)ν 與比值法、重新標度與上臨界維度、Flory）與 `examples/ex21_scaling_collapse.py`；`tools/sync_from_walk017.py` 多同步 `walk_scale.py`（第一個用相對 import 引用其他 walk 模組的同步檔），`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.19.0（2026-10-18）：加入 walk017 第 9 章（Scaling of walks and critical phenomena，動畫 W901–W906）的 `walk_crit`（平均場與 Onsager、小 Ising 的精確列舉與 Metropolis 翻轉時間、臨界指數與標度律、Widom 形式、Landau 與 (9.34)、係數漸近與 ⟨R²⟩ ∝ N^(2ν)、輔助場、高斯積分與 Wick、Trotter、n → 0 的梯子）與 `examples/ex22_ising_scaling.py`；`tools/sync_from_walk017.py` 多同步 `walk_crit.py`（用相對 import 引用 `walk_saw`、`walk_scale`），`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.20.0（2026-10-19）：加入 walk017 第 10 章（Walks and the O(n) model: mean field theory and spin waves，動畫 W1001–W1005）的 `walk_cond`（環面上的自迴避漫步精確列舉、平均場積分與 n → 0、註 3 的對照、自旋波矩陣的特徵值與逆矩陣、牆邊的 tanh、相關洞的 Coulomb 形與 Ewald 常數、erfc 形的步數、Redner–Reynolds 的 ρ 與指數）與 `examples/ex23_torus_saw.py`；`tools/sync_from_walk017.py` 多同步 `walk_cond.py`（用相對 import 引用 `walk_saw`、`walk_crit`、`walk_rec`、`walk_gf`），`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.21.0（2026-10-19）：加入 walk017 第 11 章（Scaling, fractals, and renormalization，動畫 W1101–W1106）的 `walk_rg`（Koch 曲線與 Cantor 集的維度、Landau／Gaussian 的標度塌縮、一維 Gaussian 鏈與 Ising 鏈的抽格、動量殼與連續重整化的流、線性化與臨界指數）與 `examples/ex24_decimation_flow.py`；`tools/sync_from_walk017.py` 多同步 `walk_rg.py`（用相對 import 引用 `walk_saw`、`walk_crit`），`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.35.2（2026-10-26）：只修錯。從 symosc 同步 `symosc.py`：`precession()` 的進動方向原本反了（$\cos(2\pi j/3+\delta t/2)$ → $\cos(2\pi j/3-\delta t/2)$，圖樣往 1→2→3 轉），
  `check()` 加 RK4 積分陀螺方程核對方向，檔尾改成 `main()`；`osc_draw.py` 不變。`tools/sync_from_symosc.py` 可以用參數指定 symosc 的 `codes/utils` 位置（預設仍是原機器路徑）。
- 1.35.1（2026-10-26）：加 `tools/make_bundle.py`：combkit 連同 ManimGL 1.7.2 原始碼、環境文件、固定版本清單與 `install.bat` 打包成一個 zip（程式與 1.35.0 相同）。
- 1.35.0（2026-10-23）：加入 harter_sym 第 8 章（半古典與量子動力學，全書最後一章，動畫 H801–H813）的 `sym_dyn`（直線族與圓族的包絡、Legendre 轉換；拋射體的 $S_P$、$S_H$（閉式 ＝ 數值積分）、焦散與正交的等作用量波前；$e^{iS/\hbar}$ 的半古典極限（sympy）；截斷的升降算子、平移／推動與 BCH（照印 $+\alpha^*a$ 不么正）、相干態的 Poisson 振幅、波包閉式 ＝ 本徵函數和 ＝ FFT 分裂算子；快照合成 $\langle x|5\rangle$ 與誤差；自相關與 FTAC、Maslov 與 EBK、$J=10$ 轉子的自旋相干態半古典傳播（重現書上 53.134 與 (8.3.31)）；$A\cdot p\leftrightarrow E\cdot r$ 的 sympy 證明、譜函數、黃金律、振子強度與氫的 TRK（閉式 ＝ 盒中對角化）；兩能階的曲柄與自旋進動、轉動座標、Rabi（閉式 ＝ ODE）、Bloch–Siegert（Floquet 數值）、Bloch 方程穩態、熱平衡、穿衣態與 AC Stark、Fig. 8.5.4 的圓、絕熱跟隨 ＝ Landau–Zener；光子矩陣元、模密度、Einstein $A$（SI ＝ Gauss ＝ 原子單位；氫 2p）、Jaynes–Cummings 的 2×2 塊、$g\sqrt N$ 分裂與 Mollow 三重線；全書回顧表 `RECAP`）；`draw/sym_draw.py` 加 `polyline`、`line_family`、`trajectory_fronts`、`PhasePhasor`、`snapshot_stack`、`CrankSphere`、`dressed_ladder`（第 1–7 章的元件不變）；`examples/ex37_coherent_packet.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_dyn.py`；`tests/run_selftests.py` 多檢查它的成功訊息。
- 1.34.0（2026-10-23）：加入 harter_sym 第 7 章（連續轉動群的表示乘積，動畫 H701–H714）的 `sym_coup`（$R_3$ 耦合係數的四種算法（sympy、Racah、降算子、Schwinger 玻色子）、$3j$ 對稱與 $6j$；氫的超精細與 21 cm、$(2p)^2$ 的 Pauli 項、自旋–軌道與 LS／jj 相關、向量錐的古典範圍；向量的純量／向量／張量積、張量算子與單位張量表（Table 7.1–7.4）、$R_3$ Wigner–Eckart、$j=2,3$ 的晶體場（O 標記；WE ＝ 球諧積分）、氫的 $\langle r^4\rangle$、$d^2$ 的 $6j$；剛體轉子的 RE 曲面、$J=10$ 的叢與超精細、$D_2\supset C_2$ 誘導、相關圖；八面體轉子 $J=30$ 的 $C_4$／$C_3$ 叢與分裂樣式、四階／六階混合的 12 態 $C_2$ 叢；旋量、軸角演化、Coriolis 振子、Stokes 向量與偏振橢圓、Faraday／雙折射、Rabi 拍、Lissajous 閉合、非諧局域模；H$_2^+$ 的 LCAO（閉式＝球體座標積分）、H$_2$ 的 Heitler–London、VB／MO／廣義本徵值三種 4×4、兩中心 Coulomb 積分的多極展開、混成角、Hückel 苯與丁二烯、八面體 σ／π 鍵的誘導表示）；`draw/sym_draw.py` 加 `SpinTop`、`SpinPair`、`vector_triangle`、`CouplingCone`、`RESurface`、`PoincareSphere`、`polarization_ellipse`、`ellipse_from_stokes`、`projected_curve`、`box_frame`、`potential_curves`、`min_marker`（第 1–6 章的元件不變）；`examples/ex36_re_surface.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_coup.py`；`tests/run_selftests.py` 多檢查它的成功訊息。
- 1.33.0（2026-10-23）：加入 harter_sym 第 6 章（表示的乘積，動畫 H601–H614）的 `sym_prod`（八面體上兩粒子的外積群、對角子群與相互作用的分裂；投影法的 Clebsch–Gordan 係數（書上的相位約定，$O_h$ 四方、$O\supset D_4\supset D_2$、$O\supset D_3\supset C_2$ 七張表逐格相同）與 (6.1.27)–(6.1.29) 的未知數；Frobenius–Schur 三型、對稱／反對稱平方、$T$ 群的重複 irrep、正交完備與因式分解引理；單位張量、對稱化多項式、應力應變、彈性常數的計數（81→36→21→3，暴力法＝張量定理）、模數、鈣鈦礦彈簧模型、Kittel 常數的各向異性；$C_{3v}$ 的電場、磁場與伸縮常數；參數定理與 Wigner–Eckart（六井的偶極矩陣元、複係數的共軛）；Lorentz 振子（解析＝scipy）、Smith 圓、介電函數、輻射衰減與 Rayleigh；XF$_6$ 的紅外／Raman 選擇規則；泛音的對稱冪；$E\otimes e$ Jahn–Teller 的位能面、Renner 極值與振電能階（二維振子基底＝徑向方程）；習題 6.1–6.6）；`draw/sym_draw.py` 加 `PairOcta`、`PairGrid`、`StrainCube`、`field_glyph`、`current_loop`、`trace_curve`、`SmithChart`、`hat_contours`、`hat_profile`（第 1–5 章的元件不變）；`examples/ex35_clebsch_gordan.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_prod.py`；`tests/run_selftests.py` 多檢查它的成功訊息。
- 1.32.0（2026-10-23）：加入 symosc 專案（D3／C3／C2／Q8 與耦合振盪子，動畫 S01–S03）的 `symosc`、`draw/osc_draw.py`、`examples/ex34_symmetric_oscillators.py` 與 `tools/sync_from_symosc.py`。
- 1.31.0（2026-10-23）：加入 harter_sym 第 5 章（連續轉動群的表示與應用，動畫 H501–H516）的 `sym_rot`（$O_3$、$O_2$、$D_{\infty h}$；Euler ↔ Darboux ↔ 四元數／SU(2) 的互換與 Fig. 5.3.5 的 720° 序列、
  實驗室與本體的反序；生成元的譜分解與小轉動的交換子；任意 $j$ 的 $J$ 與 Wigner $\mathcal D^j$（階乘公式＝矩陣指數＝Schwinger 振子，含 (5.4.35) 的更正）；Stern–Gerlach 機率與古典錐；
  球形／對稱／不對稱轉子與 Zeeman、Stark；球諧與加法定理；$O_3\downarrow O_h$、$D_6$、$T_d$、$I_h$ 相關與助記輪、$d$ 軌域、多極展開、SF$_6$ 型能階叢；$D_6$、$D_2$、$O$ 的射線代數；
  $R_4$ 的 Coulomb 對稱（Runge–Lenz、$n^2$ 簡併）與 $U_3$ 振子殼層、$SU_3$ 根與權；習題 5.1.1–5.6.2）；`draw/sym_draw.py` 加 `TriadBall`、`VectorCone`、`cone_family`、`ProbBars`、`ladder_column`、
  `weight_diagram`、`mnemonic_wheel`、`wheel_marker`（第 1–4 章的元件不變）；`examples/ex33_wigner_d.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_rot.py`；`tests/run_selftests.py` 多檢查它的成功訊息。
- 1.30.0（2026-10-22）：加入 harter_sym 第 4 章（八面體對稱與誘導表示，動畫 H401–H413）的 `sym_oct`（$O$、$O_h$、$T_d$、$T$、$T_h$ 由轉動矩陣產生、Hamilton 弧與 $S_4$ 置換；類代數與最小方程、特徵標表；四條子群鏈（四方／三方駐波、$C_4$／$C_3$ 行進波）的冪等元分裂與不可約表示；相關表、共軛與正規子群；誘導表示與 Frobenius 互反；六點穿隧模型與外場 $V_Q$、$V_D$、磁場、(111) 電場的能階掃描；XF$_6$ 的 21 座標簡正模＝對稱公式、$T_{1u}$ 的 2×2、UF$_6$／SF$_6$ 的 $k/j$ 曲線與同位素位移；升降算子與簡併度；習題 4.1.1–4.4.2）；`draw/sym_draw.py` 加 `CubeView`、`OctaSites`、`level_diagram`、`XF6Model`、`rot_matrix`、`irrep_color`（第 1–3 章的元件不變）；`examples/ex32_octahedral_group.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_oct.py`。
- 1.29.0（2026-10-22）：加入 harter_sym 第 3 章（非 Abel 對稱群，動畫 H301–H312）的 `sym_nonabel`（通用的有限群工具 `FiniteGroup`：乘法表、共軛類、類代數與最小方程、全對易冪等元、子群分裂的基本算子與不可約表示、正則表示的約化、特徵標表（類代數＝跡＝外積）、頻率公式、相關表；$C_{3v}$、$D_3$、$D_4$、$C_{4v}$、$D_6$、$D_{3h}$、$D_{3d}$、$D_{2h}$、$D_{4h}$、四元數群；Hamilton 轉（半角弧＝四元數＝SU(2)）；六自由度彈簧質量模型與圓周運動；$D_4$／$D_6$ 轉子能帶的對稱標記與 Zeeman；NH₃ 的 12×12 簡正模、2×2 久期方程與擬合；習題 3.1.1–3.5.6）；`draw/sym_draw.py` 加 `mirror_lines`、`BobTriangle`、`SpringTriangle`、`sym_table`、`perm_grid`、`View3D`、`SphereView`、`NH3Model`、`polar_curve`（第 1–2 章的元件不變）；`examples/ex31_character_table.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_nonabel.py`。
- 1.28.0（2026-10-22）：加入 harter_sym 第 2 章（Abel 對稱群，動畫 H201–H212）的 `sym_abel`（群表與公理、群個數、$C_2$／$C_3$／$C_n$／$D_2$ 投影算符（最小方程＝特徵標）、耦合擺的拍與 Lissajous、相量、部分約化、單擺環色散與混疊、對稱破缺的能帶、奔馳波、8 階 Abel 群與點群、$C_{12}$ 轉子與 Bloch 波、緊束縛、Mathieu、Kronig–Penney、氨兩能階與 Landau–Zener、方井 $S$ 矩陣、穿越矩陣與 WKB、習題 2.3.3–2.12.3）；`draw/sym_draw.py` 加 `PendulumRow`、`PendulumRing`、`PhasorClock`、`spring`、`band_dots`、`band_curve`、`level_lines`、`quad_contours`（第 1 章的元件不變）；`examples/ex30_cn_ring.py`；`tools/sync_from_harter_sym.py` 多同步 `sym_abel.py`。
- 1.27.0（2026-10-22）：加入 harter_sym 專案（Harter《Principles of Symmetry, Dynamics, and Spectroscopy》第 1 章，動畫 H101–H109）的 `sym_mat`（分析器與 Axiom 4、么正性、換基底、兩態拍頻、久期方程、Hamilton–Cayley 投影算符與譜分解、最小方程、對易矩陣、伴隨矩陣、散射與密度矩陣、耦合振子、d/dx 與 δ 序列、預解式的圍道積分、習題 1.1.1–1.2.8；第一個用 sympy 的模組）、`draw/sym_draw.py`、`examples/ex29_eigen_projectors.py` 與 `tools/sync_from_harter_sym.py`；`combkit/__init__.py` 的說明（之前存檔時被轉碼弄壞）重寫。
- 1.26.0（2026-10-20）：youngart 動畫版（60 支 5 秒 4K 動畫＋合成）的共用工具：`art_kit` 的 `animate_plate` 多 `prep`，加 `keyed_reveal`、`mirror`、`twinkle`、`pulse`；`young_art` 加 `turn3d`；`examples/ex28_young_wallpaper.py` 加 `--anim`。
- 1.25.0（2026-10-20）：加入 youngart（Fulton《Young Tableaux》全書 60 張 4K 藝術桌布）的 `young_extras`（計算）、`draw/young_art.py`（`YPlate` 與方格／外框／投影工具）與 `examples/ex28_young_wallpaper.py`；`tests/run_selftests.py` 多檢查 `young_extras` 的成功訊息。youngart 自己保留一份 `codes/art/ykit.py`，讓 60 張桌布可以原樣重現。
- 1.24.0（2026-10-19）：`draw/art_kit.py` 加入動畫（`animate_plate`、`render_video`、`travelers`、`grow_points`、`keep`、`spin`、`win`）與 `examples/ex27_art_animation.py`，來自 walkart 動畫版（60 張桌布各做成 5 秒動畫＋合成）。
- 1.23.0（2026-10-19）：加入 walkart（《Elements of the Random Walk》全書 60 張 4K 藝術桌布）的 `draw/art_kit.py` 與 `examples/ex26_art_wallpaper.py`；這是 combkit 第一個不靠 `manimgl` 指令、用 `render_png()` 直接輸出靜態 PNG 的工具。
- 1.22.0（2026-10-19）：加入 walk017 第 12 章（More on the renormalization group，全書最後一章，動畫 W1201–W1206）的 `walk_eps`（兩變數玩具的積掉自由度、殼層體積、配對列舉與連結叢集、O(n) 縮併、一迴圈遞迴的不動點與線性化、ν 的 ε 展開與 Flory、(12.86) 與冪次數算、全書回顧的數字）與 `examples/ex25_wilson_fisher.py`；`tools/sync_from_walk017.py` 多同步 `walk_eps.py`（用相對 import 引用 `walk_saw`、`walk_crit`、`walk_rg`），`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.16.0（2026-10-18）：加入 walk017 第 6 章（The shape of a random walk，動畫 W601–W607）的 `walk_shape`（迴轉半徑張量、Kramers 恆等式、非球度、1/d 展開、Fixman）與 `examples/ex19_walk_shape.py`；`tools/sync_from_walk017.py` 多同步 `walk_shape.py`，`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.15.0（2026-10-18）：加入 walk017 第 5 章（Variations on the random walk，動畫 W501–W507）的 `walk_var`（第一個用到 numpy／scipy.special 的數學模組）與 `examples/ex18_persistent_walk.py`；`tools/sync_from_walk017.py` 多同步 `walk_var.py`，`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.14.0（2026-10-18）：加入 walk017 第 4 章（Boundary conditions, steady state，動畫 W401–W408）的 `walk_bc` 與 `examples/ex17_image_walker.py`；`tools/sync_from_walk017.py` 多同步 `walk_bc.py`，`tests/run_selftests.py` 多檢查它的成功訊息（沒有新的畫圖模組）。
- 1.13.0（2026-10-17）：加入 walk017 第 3 章（Generating functions II，動畫 W301–W308）的 `walk_rec` 與 `examples/ex16_recurrence.py`；`tools/sync_from_walk017.py` 多同步 `walk_rec.py`（沒有新的畫圖模組）。
- 1.12.0（2026-10-17）：加入 walk017 第 2 章（Generating functions I，動畫 W201–W212）的 `walk_gf`、`draw/lattice_draw.py` 與 `examples/ex15_lattice_gf.py`；`tools/sync_from_walk017.py` 多同步這兩個模組。
- 1.11.0（2026-10-17）：加入 walk017 專案（Rudnick & Gaspari《Elements of the Random Walk》第 1 章，動畫 W101–W112）的 `walk1d`、`draw/walk_draw.py`、`examples/ex14_random_walk.py` 與 `tools/sync_from_walk017.py`；`combkit/__init__.py` 的說明補上 `schubert`、`fulton_appendix`。
- 1.10.1（2026-10-16）：更新附的環境文件 `環境與交接_ENVIRONMENT.md`（Fulton 全書 F01–F65、可攜版 106）；程式與 1.10.0 相同。
- 1.10.0（2026-10-16）：加入 Fulton 附錄 A–B（`fulton_appendix`、`examples/ex13_evacuation_keys.py`）；來自 catalan930 的動畫 F56–F65。
- 1.9.0（2026-10-15）：加入 Fulton 第 9–10 章（`schubert`、`draw/schubert_draw.py`、`examples/ex12_schubert.py`）；`tools/render_check.py` 每次試渲染改用自己的 TEMP 資料夾，兩個 render_check 可以同時跑。
- 1.8.0（2026-10-15）：加入 `raney_examples`（例題產生器：均勻隨機的字／樹／不交叉分法、Markdown 例題集）、`draw/raney_draw.py` 的字／樹／秩路徑工具與 `examples/ex11_random_trees.py`；來自 catalan930 的動畫 N01–N12（Raney 導讀）。
- 1.7.0（2026-10-15）：加入整數序列的 Raney 引理（`raney` 第 9 節、`draw/raney_draw.py`、`examples/ex10_raney_lemma.py`），來自 catalan930 的動畫 R01–R06（`notes/ran_note015_改寫.pdf`）。
- 1.6.0（2026-10-14）：加入 `rmt_moments`（Wick 配對、GOE/GUE 精確動差、樹形閉路徑）、`tikz_trees`（平面樹的 TikZ）與 `examples/ex09_wick_moments.py`；來自 spglass 專案的 `codes/make_solution_tikz.py`（先修講義練習解答的圖），該檔現在改用這兩個模組。
- 1.4.0（2026-10-12）：加入 Fulton 第 8 章的範例 `examples/ex07_schur_module.py`（第 8 章沒有新增共用模組）。
- 1.3.0（2026-10-12）：加入 Fulton 第 7 章（`rep_data` 與資料檔、`draw/rep_draw.py`、`examples/ex06_specht.py`）；`render_check.py` 加 `--tile`。
- 1.2.0（2026-10-10）：加入 Fulton 第 5–6 章（`tableaux` 第 8 節：LR 規則、對稱多項式；`examples/ex05_lr_rule.py`）。
  catalan930 仍使用自己的那一份；兩邊若要一致，改 catalan930 後執行 `tools/sync_from_catalan930.py`。
- 程式由 Claude（Anthropic）在 catalan930 專案中撰寫。數學內容依據：W. Chu (1987)、A. Selimi & M. Saračević (2019)、
  MIT 18.212 Lecture 1–4（使用者提供的整理）、W. Fulton《Young Tableaux》、G. N. Raney (1960)、I. M. Gessel (2016)；`rmt_moments` 依據標準的隨機矩陣動差法（Anderson–Guionnet–Zeitouni 第 2.1 節）與 Harer–Zagier (1986)。這些參考文獻本身不包含在工具包裡。
- ManimGL 為 MIT 授權（3Blue1Brown 與貢獻者），不包含在本工具包裡。

---

## 10. 加入新的演示（今後的擴充方式）

以後做的新演示（不論在哪個計畫），可重複使用的程式都收進 combkit。做法：

1. **數學**放 `combkit/<主題>.py`（純 Python，不 import manim），檔尾加 `if __name__ == "__main__":` 自我測試，
   用暴力列舉交叉驗證；最後印一句「…全部驗證通過」。
2. **畫圖**放 `combkit/draw/<主題>_draw.py`；import 一律寫 `from combkit.X import …`／`from combkit.draw.X import …`。
3. **範例**：在 `examples/` 加一個小場景 `exNN_<主題>.py`，示範「先算、assert、再畫」，並在第 5 節的表加一行。
4. **說明**：在第 4 節加該模組的介紹、第 3 節加新的資料表示法（如果有）；在 `combkit/__init__.py` 的說明列出新模組。
5. **版本**：`combkit/__init__.py` 的 `__version__` 加一版（新功能 1.x→1.(x+1)，只修錯 1.0.0→1.0.1），在第 9 節記一行。
6. **檢查與打包**（都在本資料夾執行）：
   ```bash
   python tests/run_selftests.py --draw
   ```
   ```bash
   python tools/make_api.py
   ```
   ```bash
   python tools/make_zip.py
   ```
   `run_selftests.py` 會自動找到新模組，不用改它。

catalan930 裡新增的演示則照舊：程式寫在 catalan930 的 `codes/utils/`，再執行 `tools/sync_from_catalan930.py`
（新模組會自動帶進來），然後做上面的 3–6。

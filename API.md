# combkit API 索引

由 `tools/make_api.py` 用 ast 從原始碼自動產生；說明取自各函式 docstring 的第一行。手寫的介紹請看 README.md。

## `combkit.counting`

counting.py —— 純 Python 的計數工具（不依賴 manim）

常數：`LEAF`、`PATTERNS_3`、`MIT_EXAMPLE_TREE`

- `catalan(n: int)` — 經典卡塔蘭數 C_n = binom(2n, n) / (n + 1)。
- `inner(n_vec, v_vec)` — 內積 <n, v> = sum n_i v_i，也就是圓周上的總點數 N。
- `closed_form(n_vec, v_vec)` — 論文式 (2)：
- `blocks_cross(a, b)` — 兩個群組 a, b（點的編號集合）是否交叉。
- `is_noncrossing(blocks)` — 
- `enumerate_noncrossing(n_vec, v_vec)` — 列出所有不交叉分組。回傳 list，每個元素是 list[tuple]（各群組的點編號，0 起算）。
- `brute_force(n_vec, v_vec)` — 
- `recurrence(n_vec, v_vec)` — d(n) = sum_i sum_{r_1+...+r_{v_i} = n - e_i} d(r_1) ... d(r_{v_i})，d(0) = 1。
- `step_weight(letter: int, v_vec)` — 
- `path_heights(word, v_vec, start: int=0)` — 回傳 [S_0, S_1, ..., S_M]，S_0 = start。
- `is_good_path(word, v_vec)` — 
- `first_bad_index(word, v_vec)` — 第一個高度 <= 0 的點（不含起點）；好路徑回傳 None。
- `multinomial_vec(M: int, n_vec)` — binom(M; n) = M! / (n_1! ... n_k! (M - |n|)!)；不合法時回傳 0。
- `raney_A(t: int, n_vec, v_vec)` — A_t(n) = t / (N + t) * binom(N + t; n)。A_1 就是廣義卡塔蘭數。
- `word_content(word, k: int)` — 路徑用了幾個第 1, ..., k 型步。
- `path_words(t: int, n_vec, v_vec)` — 所有排法（依字典序排列，不含重複）。
- `good_paths(t: int, n_vec, v_vec)` — 
- `rotate(word, j: int)` — 從第 j 步開始讀（循環）。
- `good_rotations(word, v_vec)` — M 種讀法中，哪幾種是好路徑（回傳起點 j）。
- `cut_index(word, s: int, v_vec)` — 最後一個高度 <= s 的點；好路徑在這裡切成「終點 s」與「終點 t」兩段。
- `sub_vectors(n_vec)` — 所有 0 <= r <= n（逐分量）的向量 r，依字典序。
- `binary_trees(n: int)` — n 個內部節點的所有滿二元樹，依「左子樹大小」由小到大排列。
- `tree_size(t)` — 內部節點的個數。
- `tree_to_paren_tokens(t, factors='abcdefghij')` — 完全括號化（最外層不加括號），用 · 表示每一次乘法。
- `parenthesizations(n: int, factors='abcdefghij')` — n+1 個因子的所有完全括號化（字串，最外層不加括號）。
- `parenthesizations_bruteforce(seq: str)` — 獨立的列舉：把因子串切成左右兩段、各自遞迴（字串形式與上面相同）。
- `tree_to_triangulation(t, i: int, j: int)` — 凸多邊形頂點 i, i+1, ..., j（底邊為 i–j）的三角剖分，回傳對角線串列。
- `triangulation_apex(t, m: int)` — 凸 m 邊形（頂點 0..m-1，底邊 0–(m-1)）中，含底邊的那個三角形的第三個頂點。
- `triangulations(m: int)` — 凸 m 邊形的所有三角剖分（對角線的 tuple），依「含底邊三角形的頂點」由左到右排列。
- `diagonals_cross(a, b)` — 
- `triangulations_bruteforce(m: int)` — 獨立的列舉：從所有對角線裡挑 m-3 條互不相交的。
- `tree_to_dyck(t)` — Dyck 路徑（U 上坡、D 下坡）：t = (L, R) 對應 U·dyck(L)·D·dyck(R)。
- `dyck_parts(t)` — Dyck 路徑每一步屬於 'root'（第一個上坡與它對應的下坡）、'left' 或 'right'。
- `is_dyck(word: str)` — 
- `dyck_words(n: int)` — 獨立的列舉：所有 n 個 U、n 個 D、從不低於 0 的字（字典序）。
- `dyck_to_lattice(word: str)` — Dyck 上坡 -> 格子路徑往右 R；下坡 -> 往上 U。（不越過 y=x 的路徑。）
- `lattice_words(n: int)` — (0,0) 到 (n,n) 的所有單調格子路徑（n 個 R、n 個 U），字典序。
- `lattice_points(word: str, start=(0, 0))` — 
- `lattice_ok(word: str)` — 從不走到對角線 y = x 的上方（可以碰到對角線）。
- `first_bad_step(word: str)` — 第一次走到 y = x + 1 的那一步之後的位置（點的編號）；合法路徑回傳 None。
- `reflect_tail(word: str)` — 反射原理：第一次碰到 y = x + 1 之後的部分，沿著 y = x + 1 翻轉（R 與 U 互換）。
- `dyck_to_brackets(word: str)` — 
- `dyck_to_ballot(word: str)` — 
- `ballot_ok(word: str)` — 卡塔蘭特例：開票過程中 A 的票數始終 >= B 的票數，而且最後同票。
- `bertrand_strict_count(a: int, b: int)` — 暴力列舉：A 得 a 票、B 得 b 票，開票過程中 A 始終「嚴格」領先的順序數。
- `ballot_more_than_k_times(a: int, b: int, k: int)` — 暴力列舉：開票過程中 A 的票數始終「大於」B 的 k 倍的順序數（推廣的投票問題）。
- `dyck_to_matching(word: str)` — 每個上坡與它對應的下坡配成一對：圓周 2n 點的不交叉完美配對。
- `matching_parts(t)` — 配對裡每條弦屬於 'root'（含第 0 點的弦）、'left'（弦內側）或 'right'（外側）。
- `series_square_shift(coeffs: list)` — 回傳 1 + x * C(x)^2 的前 len(coeffs) 個係數。
- `series_mul(a, b, terms: int | None=None)` — 兩個冪級數（係數串列）的 Cauchy 乘積前 terms 項：[x^n](ab) = sum_{k=0}^{n} a_k b_{n-k}。
- `binom_general(r, m: int)` — 廣義二項係數 binom(r, m) = r (r-1) ... (r-m+1) / m!（r 可以是分數，例如 1/2）。
- `sqrt_1m4x_coeffs(terms: int)` — sqrt(1-4x) = sum_m binom(1/2, m) (-4)^m x^m 的前 terms 個係數：1, -2, -2, -4, -10, ...（都是整數）。
- `catalan_from_sqrt(terms: int)` — C(x) = (1 - sqrt(1-4x)) / (2x)：C_n = -[x^{n+1}] sqrt(1-4x) / 2。
- `drunkard_fall_words(n: int)` — 暴力列舉：長度 2n+1 的 U/D 字，從 0 出發、在最後一步「第一次」走到 -1（掉下懸崖）。字典序。
- `drunkard_fall_prob(n: int)` — 第 2n+1 步第一次掉下去的機率 = C_n / 2^(2n+1)（Dyck 路徑接一步往下）。
- `drunkard_partial_sums(N: int)` — [S_1, ..., S_N]，S_K = sum_{n<K} C_n / 2^(2n+1) = 1/2, 5/8, 11/16, ...（精確分數，趨近 1）。
- `cycle_words(n: int)` — n+1 個 +1（字母 0 = LEAF）與 n 個 -1（字母 1）的所有排法，長度 2n+1，共 binom(2n+1, n) 個。
- `cycle_classes(n: int)` — 依循環移位把 cycle_words(n) 分組：[(好的字, [rotate(好的字, j) for j in 0..2n]), ...]，
- `cycle_good_to_dyck(word)` — 好的字 = 一個 +1 接一條 Dyck 路徑：去掉第一個字母，0 -> 'U'、1 -> 'D'。
- `binary_to_plane(t)` — 滿二元樹（n 個內部節點）-> 平面樹（n+1 個頂點）：左小孩 = 第一個小孩、右小孩 = 下一個兄弟。
- `plane_to_binary(pt)` — binary_to_plane 的反函數。
- `plane_tree_edges(pt)` — 邊數（= 頂點數 - 1）。
- `plane_trees(n: int)` — n 條邊的所有平面樹，順序與 binary_trees(n) 相同（所以和 00k 的欄位一致）。
- `plane_trees_bruteforce(n: int)` — 獨立的列舉：根有 k 個小孩，剩下 n-k 條邊任意分給 k 棵小孩子樹。
- `plane_tree_dfs(pt)` — 深度優先走訪：走向小孩記 'U'、回到父親記 'D'。plane_tree_dfs(binary_to_plane(t)) == tree_to_dyck(t)。
- `plane_tree_preorder(pt)` — 頂點依前序（DFS 第一次到達的順序）編號 0, 1, ...（0 是根）；回傳 [(父親編號或 None, 深度), ...]。
- `plane_tree_walk(pt)` — 沿輪廓走一圈依序經過的頂點（前序編號），長度 2n+1；第 k 步從 walk[k] 走到 walk[k+1]。
- `triangulation_dual(t)` — tree_to_triangulation(t, 0, m-1)（m = n+2，底邊 0–(m-1)）的對偶樹。每個三角形一筆，依「中序」排列
- `triangulation_to_tree(m: int, diagonals)` — 獨立的反方向：只看對角線重建對偶二元樹。從底邊 (0, m-1) 出發，找含這條邊的三角形 (a, k, b)，
- `perfect_matchings(n_points: int)` — n_points 個點（0 起算）的所有完美配對，共 (n_points-1)!! 個；每個是排序好的 ((i, j), ...)，i < j。
- `arcs_cross(p, q)` — 直線上兩條弧交叉 <=> 端點交錯 i < j < k < l：(i, k)、(j, l)。
- `arcs_nest(p, q)` — 直線上兩條弧巢狀 <=> 一條在另一條裡面 i < j < k < l：(i, l)、(j, k)。
- `matching_crossings(pairs)` — 交叉的弧：[(第 x 條, 第 y 條), ...]（pairs 裡的位置）。
- `matching_nestings(pairs)` — 巢狀的弧：[(外面那條, 裡面那條), ...]（pairs 裡的位置）。
- `matching_is_noncrossing(pairs)` — 
- `matching_is_nonnesting(pairs)` — 
- `dyck_to_nonnesting(word: str)` — Dyck 路徑 -> 不巢狀配對：上坡依序排隊，每個下坡配給最早、還沒配的上坡（先進先出）。
- `matching_to_dyck(pairs)` — 配對 -> 字：每條弧的左端記 'U'、右端記 'D'（不交叉或不巢狀配對都適用）。
- `stack_sort(perm)` — 用一個堆疊排序，目標輸出 1, 2, ..., n（mit_plan.md 第 4 節 M09 的規則），每一步：
- `queue_bypass_sort(perm)` — 用「佇列＋插隊線」排序（mit_plan.md 第 1 節第 5 點），每一步：
- `sort_states(perm, ops)` — 重播 ops：回傳每一步之後的 (輸入, 容器, 輸出)，第 0 筆是開始時。
- `standardize(seq)` — 相對大小：(4, 5, 3) -> (2, 3, 1)。
- `contains_pattern(w, p)` — w 含模式 p（有子序列和 p 的相對大小相同）時，回傳第一個出現的位置（0 起算的 index tuple，
- `avoiders(n: int, p)` — 1..n 的排列中避開 p 的（依字典序）。
- `integer_partitions(n: int)` — n 的所有整數分割（由大到小的 tuple），反字典序：(5,), (4, 1), (3, 2), (3, 1, 1), ...。
- `partition_count(n: int)` — 獨立的算法：Euler 五角數定理的遞迴 p(m) = sum_k (-1)^(k+1) [p(m - k(3k-1)/2) + p(m - k(3k+1)/2)]。
- `integer_compositions(n: int)` — n 的所有合成（有順序的正整數和），共 2^(n-1) 個（n >= 1）。字典序。
- `conjugate(lam)` — 共軛分割 λ'（Young 圖轉置）：λ'_j = 第 j 行的長度。
- `young_cells(lam)` — 所有格子 (i, j)（第 i 列、第 j 行，0 起算），逐列由左到右。
- `hook_cells(lam, i: int, j: int)` — 鉤上的格子：[(i, j) 自己, 同列右邊（arm）..., 同行下面（leg）...]。
- `hook_length(lam, i: int, j: int)` — 鉤長 h(i, j) = arm + leg + 1（0 起算）。
- `hook_lengths(lam)` — 每一列的鉤長，例如 (4, 3, 1) -> [(6, 4, 3, 1), (4, 2, 1), (1,)]。
- `hook_product(lam)` — H(λ) = 所有格子鉤長的乘積。
- `hook_formula(lam)` — 鉤長公式 f^λ = n! / H(λ)。
- `corners(lam)` — 角落（拿掉後仍是 Young 圖的格子）：[(i, λ_i - 1), ...]，由上到下。
- `remove_corner(lam, c)` — λ - c：拿掉角落 c（列長變 0 就去掉那一列）。
- `syt_count(lam)` — 角落遞迴：最大的數 n 一定在角落 -> f^λ = sum_{角落 c} f^{λ-c}，f^() = 1。（λ 可以是 list 或 tuple）
- `syt(lam)` — 形狀 λ 的所有標準 Young 表（每列 →、每行 ↓ 都遞增），依「列字」的字典序。
- `syt_bruteforce(lam)` — 獨立的列舉：把 1..n 的每一種排列逐列填入，留下列、行都遞增的。
- `syt_to_rowword(T)` — 列字：第 k 個字母 = 數字 k+1 在第幾列（0 起算）。
- `rowword_to_syt(word)` — syt_to_rowword 的反函數：依序把 1, 2, ... 放到列字指定那一列的尾端。
- `rowword_counts(word, k: int | None=None)` — 累計次數表（k 維格子路徑）：第 t 筆 = 讀完前 t 個字母後各列的個數 (x_0, ..., x_{k-1})，第 0 筆全為 0。
- `is_lattice_rowword(word)` — 
- `syt_to_dyck(T)` — 形狀 (n, n)：數字 i 在第一列 -> 第 i 步 'U'、在第二列 -> 'D'。
- `dyck_to_syt(word: str)` — syt_to_dyck 的反函數：'UUDUDD' -> ((1, 2, 4), (3, 5, 6))。
- `corner_prob_formula(lam, c)` — (1/n) H(λ) / H(λ - c)。對所有角落加總 = 1 <=> 鉤長公式（對 n 歸納）。
- `hook_walks(lam, start)` — 從格子 start 出發的所有鉤上漫步：每步跳到目前格子鉤上的另一格（機率 1/(h-1)），停在角落（h = 1）。
- `hook_walk_probs(lam)` — {角落: 停在這個角落的機率}（精確分數）：起點均勻 1/n，列舉所有漫步加總。等於 corner_prob_formula。
- `walk_rows_cols(path)` — 漫步的投影：(經過的列集合 A, 經過的行集合 B)，都是排序好的 tuple。
- `hook_identity(lam, cell, c)` — 關鍵恆等式（c = (α, β) 是角落、cell = (a, b) 在 c 的左上方）：
- `hook_lemma_rhs(lam, c, A, B)` — 引理的右邊：從 (min A, min B) 出發、經過的列集合 A、行集合 B、停在角落 c = (α, β) 的所有漫步權重和
- `corner_prob_product(lam, c)` — M14 第 5 段：P(c) = (1/n) prod_{i<α} (1 + 1/(h(i,β)-1)) · prod_{j<β} (1 + 1/(h(α,j)-1))。

## `combkit.tableaux`

tableaux.py —— 純 Python 的 Young 表組合演算（不依賴 manim）

常數：`FULTON_TABLEAU`、`FULTON_STANDARD`、`FULTON_INSERT_T`、`FULTON_PRODUCT_U`、`FULTON_SKEW`、`KNUTH_RULES`、`CH2_WORD`、`CH2_SLIDE_BEFORE`、`CH2_SLIDE_AFTER`、`CH2_SLIDE_CHAIN`、`CH3_WORD`、`CH3_WORD2`、`CH3_TABLEAU`、`CH4_WORD`、`CH4_ARRAY`、`CH4_MATRIX`、`CH4_P`、`CH4_Q`、`CH4_HOOK_SHAPE`、`CH5_RL_GOOD`、`CH5_RL_BAD`、`CH5_DECOMP`、`CH5_S21_SQUARED`、`CH5_ZERO`

- `partitions(n: int, max_parts: int | None=None)` — n 的所有分割（由大到小的 tuple），反字典序；max_parts 限制列數。例：partitions(3) = [(3,), (2, 1), (1, 1, 1)]。
- `is_partition(lam)` — λ 是否為由大到小（弱遞減）的正整數序列。
- `conjugate(lam)` — 共軛分割 λ̃（沿主對角線翻轉 Young 圖）：λ̃_j = 第 j 行的長度。
- `young_boxes(lam)` — Young 圖的所有格子 (i, j)，逐列由左到右、由上到下。
- `corners(lam)` — 角（右邊與下面都沒有格子的格子）[(i, λ_i - 1), ...]，由上到下。對 λ 就是外角，對 μ 就是內角。
- `addable_boxes(lam)` — 加上去之後仍是 Young 圖的格子，由上到下（最後一個是新的一列 (len(λ), 0)）。
- `contains(lam, mu)` — μ ⊂ λ（μ 的 Young 圖包含在 λ 的裡面），也就是對所有 i 有 μ_i <= λ_i。
- `subpartitions(lam)` — 所有 μ ⊂ λ（含 λ 自己與 ()）。
- `skew_boxes(lam, mu)` — 斜圖形 λ/μ 的格子，逐列由左到右。
- `is_horizontal_strip(lam, mu)` — λ/μ 是水平帶（μ ⊂ λ，且沒有兩格在同一行）<=> λ_{i+1} <= μ_i。
- `is_vertical_strip(lam, mu)` — λ/μ 是垂直帶（μ ⊂ λ，且沒有兩格在同一列）<=> λ_i - μ_i <= 1。
- `dominates(lam, mu)` — 支配序 λ ⊵ μ：|λ| = |μ| 且每個前綴和 λ_1+...+λ_k >= μ_1+...+μ_k。
- `shape(T)` — 表（或斜表）的外形 λ = 各列長度。
- `num_boxes(S)` — 有數字的格子數（斜表不算 None）。
- `is_tableau(T)` — （半標準）Young 表：外形是分割、數字是正整數、每列弱遞增、每行嚴格遞增。
- `is_standard(T)` — 標準表：是表，且數字恰好是 1..n 各一次。
- `transpose(T)` — 轉置 T^τ（共軛圖上的編號）：T^τ 的 (j, i) = T 的 (i, j)。標準表轉置後仍是標準表，一般的表不一定。
- `content(T, m: int | None=None)` — 內容（= 單項式 x^T 的指數）：(1 出現幾次, 2 出現幾次, ..., m 出現幾次)。m 預設為最大的數。
- `ssyt(lam, m: int | None=None, content: Sequence[int] | None=None)` — 形狀 λ、數字在 [m] 的所有表，依逐列讀法的字典序。
- `standard_tableaux(lam)` — 形狀 λ 的所有標準表（= 內容為 (1, 1, ..., 1) 的表）。
- `tableaux_of_size(n: int, m: int)` — 恰有 n 格、數字在 [m] 的所有表（所有形狀）。
- `bump_steps(T, x: int)` — 列插入 T ← x 的每一步：[(格子, 放進去的數, 被擠出來的數 或 None), ...]，由上到下每列一筆。
- `row_insert(T, x: int)` — 列插入 T ← x。回傳 (新表, 擠壓路徑, 新格子)。
- `reverse_bump(T, box: Box)` — 反向擠壓（書 p.8）：從外角 box 的數 y 開始，在上一列找「最右邊、嚴格小於 y」的數，換上去並把它往上擠，
- `tableau_word(S)` — 表（或斜表）的字 w(T)：由最下面一列讀到最上面一列，每列由左到右（斜表略過 None）。
- `insert_word_steps(word, T=())` — 把字母依序列插入 T，每一步記下 (字母, 插入後的表, 擠壓路徑, 新格子)。
- `insert_word(word, T=())` — ((T ← w_1) ← w_2) ← ... ← w_k；T 預設為空表，此時就是 P(w)。
- `product(T, U)` — 乘積 T·U（§1.1）：把 U 的字 w(U)（由下而上、由左到右）依序列插入 T。
- `skew_shape(S)` — 斜表的形狀 (λ, μ)：λ_i = 第 i 列長度，μ_i = 第 i 列開頭 None 的個數。
- `skew_cells(S)` — 斜表裡有數字的格子 {(i, j): 數字}。
- `skew_from_cells(lam, mu, cells: Dict[Box, int])` — skew_cells 的反函數：給 λ、μ 與 λ/μ 每一格的數字，組回斜表。
- `is_skew_tableau(S)` — 斜表：每列是「若干 None + 正整數」，λ、μ 都是分割且 μ ⊂ λ，數字每列弱遞增、每行嚴格遞增。
- `inside_corners(S)` — 內角：μ 裡面、下面與右邊都不在 μ 的格子，由上到下。沒有內角 <=> S 已是一般的表。
- `outside_corners(S)` — 外角：λ 裡面、下面與右邊都不在 λ 的格子，由上到下（可能同時是內角，見書 p.13）。
- `forward_slide(S, corner: Box)` — 從內角 corner 做一次滑動（挖洞）。回傳 (新斜表, 洞的路徑)。
- `reverse_slide(S, box: Box)` — 反向滑動（書 p.14）：在 λ 外加上格子 box 當洞，洞往上或往左跟較大的鄰居交換（一樣大選上面），
- `jdt_steps(S, choose: Callable[[List[Box]], Box] | None=None)` — jeu de taquin：一直選內角滑動，直到沒有內角。回傳每一步 (內角, 洞的路徑, 滑動後的斜表)。
- `rectify(S, choose: Callable[[List[Box]], Box] | None=None)` — Rect(S)：jeu de taquin 做到底得到的一般表（書的 Claim 2：與內角的選法無關）。
- `star(*tabs)` — T*U（§1.2）：μ 是「T 的行數 × U 的列數」的空白矩形，T 放在矩形下面、U 放在矩形右邊（T 左下、U 右上）。
- `skew_tableaux(lam, mu, m: int)` — 形狀 λ/μ、數字在 [m] 的所有斜表，依逐列讀法的字典序。
- `schur_poly(lam, m: int)` — s_λ(x_1, ..., x_m) = Σ_T x^T（T 跑遍形狀 λ、數字在 [m] 的表）。回傳 {指數: 係數}；λ 列數 > m 時是 {}（零）。
- `complete_h(n: int, m: int)` — 完全齊次對稱多項式 h_n(x_1..x_m)：所有 n 次單項式，係數都是 1。
- `elementary_e(n: int, m: int)` — 基本對稱多項式 e_n(x_1..x_m) = Σ_{i_1 < ... < i_n} x_{i_1} ... x_{i_n}。
- `kostka(lam, mu)` — Kostka 數 K_{λμ}：形狀 λ、內容為 μ（數字 k 出現 μ_k 次）的表的個數。μ 可以是任意（弱）合成。
- `word_runs(w)` — 把字在「後一個字母嚴格變小」的地方切開，得到極大的弱遞增片段（依原來順序）。書 p.17。
- `word_to_tableau(w)` — w 若是某個表的字，傳回那個表（片段由下往上當各列），否則 None。書 p.17。
- `skew_from_word(w)` — 字為 w 的一個斜表：每個弱遞增片段當一列，每列放在下一列的整個右上方。書 p.17–18。
- `column_word(S)` — 行字 w_col(S)：由左到右逐行，每行由下往上讀（斜表略過 None）。書 p.27。
- `knuth_moves(w)` — 一步基本 Knuth 變換能到的字：[(i, 規則, 新字)]，i = 被改的三個字母的起點（0 起算）。書 p.19：
- `knuth_class(w)` — 與 w Knuth 等價的所有字（廣度優先搜尋；只適合短字）。
- `knuth_path(w1, w2)` — 從 w1 到 w2 的最短基本變換序列 [(i, 規則, 變換後的字)]；不等價時傳回 None。
- `insertion_knuth_chain(T, x)` — Knuth 對列插入的拆解（書 p.18–19）：從 w(T)·x 一步一步變到 w(T ← x)，每步一個基本變換 [(i, 規則, 新字)]。
- `max_increasing(w, k: int=1, decreasing: bool=False)` — L(w, k)（書 p.30）：從 w 取出 k 個互不相交的弱遞增子序列，長度總和的最大值，以及一組達到它的子序列
- `greene_numbers(w, decreasing: bool=False)` — L(w,1), L(w,2)-L(w,1), ...：依 Greene 定理就是 P(w) 的列長（decreasing=True 時是行長）。
- `remove_largest(w, p: int=1)` — 去掉 p 個最大的字母（相同時右邊的算比較大，先去掉）。書 p.33 Lemma 3。
- `tableau_from_greene(w)` — 第 3 章的唯一性證明當演算法（書 p.34）：只用 L(w,k) 決定形狀，一次去掉一個最大字母，
- `is_lex_array(top, bottom)` — 雙列陣列是字典序：上列弱遞增，上列相同時下列弱遞增（書 p.39 (1)(2)）。
- `rsk_steps(top, bottom=None)` — R-S-K（書 p.36–39）：陣列 (top 在上、bottom 在下，須是字典序) → 每一步 (u, v, P_k, Q_k, 新格子)：
- `rsk(top, bottom=None)` — R-S-K 對應的 (P, Q)；用法同 rsk_steps。
- `inverse_rsk(P, Q)` — 由同形狀的 (P, Q) 還原字典序陣列 (top, bottom)：每次取 Q 最大的數（相同取最右邊）那一格，對 P 反向擠出。
- `matrix_to_array(A)` — 矩陣 A（a(i,j) = (i 在上、j 在下) 出現的次數，i、j 從 1 起算）→ 字典序陣列。書 p.41。
- `array_to_matrix(top, bottom, m: int | None=None, n: int | None=None)` — 
- `rsk_matrix(A)` — 矩陣的 R-S-K：P 的數來自行號、Q 的數來自列號。
- `transpose_matrix(A)` — 
- `ball_numbering(A)` — 矩陣球的編號（書 p.42）：(i, j) 放 a(i,j) 個球（0 起算）；球的號碼 = 比它西北方所有球的號碼都大的最小整數
- `derived_matrix(A)` — 導出矩陣 A♭（書 p.43–44）：同號的 ℓ 個球由西南排到東北，除了最後一個，每個球的右邊、
- `matrix_ball(A)` — 矩陣球構造（書 §4.2）：每一輪把球編號，P 的這一列 = 各號碼最左邊的行號、Q 的這一列 = 各號碼最上面的列號
- `pieri_shapes(lam, p: int, vertical: bool=False)` — λ 加 p 格、任兩格不同行（水平條）或不同列（vertical=True，垂直條）得到的所有 μ。書 p.24 (4)(5)。
- `poly_mul(f: Poly, g: Poly)` — 
- `poly_sum(polys)` — 
- `lex_greater(lam, mu)` — 字典序 λ > μ：第一個不同的位置 λ_i > μ_i。書 p.26。
- `kostka_matrix(n: int)` — n 的所有分割（字典序由小到大，同書 p.26）與 Kostka 矩陣 K[λ][μ]（下三角、對角線為 1）。
- `hook_content_count(lam, m: int)` — Stanley 的鉤–內容公式（書 p.55）：數字在 [m] 的 λ 形表的個數 d_λ(m) = Π (m + j − i) / h(i,j)。
- `count_involutions(n: int)` — S_n 裡對合（σ² = 1）的個數 = Σ_k n! / ((n−2k)! 2^k k!)（書 p.52 Exercise 6）。
- `factorizations(V, lam, mu)` — 書 p.61 的 𝒯(λ, μ, V)：所有 [T, U]，T 形狀 λ、U 形狀 μ、T·U = V（暴力列舉）。
- `rectifying_skews(nu, lam, U0)` — 書 p.60 的 𝒮(ν/λ, U0)：形狀 ν/λ、rectification 等於 U0 的斜表。
- `is_reverse_lattice(w)` — 反格字（Yamanouchi 字，書 p.63）：從字尾往前讀到任何位置，1 的個數 ≥ 2 的個數 ≥ 3 的個數 ≥ …。
- `lattice_tableau(mu)` — U(μ)：第 i 列全是 i 的表（LR 斜表 rectify 之後的樣子，書 p.65）。
- `lr_tableaux(nu, lam, mu=None)` — ν/λ 上的 Littlewood–Richardson 斜表：字是反格字（μ 給定時只取內容為 μ 的）。
- `lr_coefficient(lam, mu, nu)` — Littlewood–Richardson 數 c^ν_{λμ} = ν/λ 上內容為 μ 的 LR 斜表個數（書 p.64 命題 3）。
- `lr_product(lam, mu)` — s_λ s_μ = Σ c^ν_{λμ} s_ν（書 p.66 (4)）：{ν: c}。
- `skew_decomposition(nu, lam)` — s_{ν/λ} = Σ_μ c^ν_{λμ} s_μ（書 p.65）：{μ: c}。
- `reverse_numbering(nu, lam)` — ν/λ 的反向編號（書 p.68）：由上到下逐列，每列由右往左編 1, 2, …。
- `remmel_whitney(nu, lam, mu)` — Remmel–Whitney（書 p.68–69 命題 4）：μ 上滿足 (i)(ii) 的標準表；個數 = c^ν_{λμ}。
- `poly_scale(f: Poly, c)` — 
- `poly_prod(polys, m: int)` — 
- `monomial_sym(lam, m: int)` — m_λ：x^λ 的所有不同排列之和（λ 的列數 > m 時是 0）。
- `power_sum(r: int, m: int)` — 
- `h_lambda(lam, m)` — 
- `e_lambda(lam, m)` — 
- `p_lambda(lam, m)` — 
- `z_lambda(lam)` — z(λ) = Π r^{m_r} m_r!（書 p.74 (3)）；n!/z(λ) 是 S_n 中循環型為 λ 的排列個數。
- `schur_expand(f: Poly, m: int)` — 把對稱多項式寫成 Σ c_λ s_λ（每次扣掉字典序最大的單項式；需要 m ≥ 次數才完整）。
- `poly_det(M, m: int)` — 多項式矩陣的行列式（Leibniz 展開；只適合小矩陣）。
- `jacobi_trudi(lam, m: int, dual: bool=False)` — (5) s_λ = det(h_{λ_i+j-i})；dual=True 時是 (6) s_λ = det(e_{λ̃_i+j-i})（書 p.75）。
- `alternant(exps, m: int)` — a(ℓ_1..ℓ_m) = det(x_j^{ℓ_i})；書 p.75 (7)：s_λ = a(λ+δ)/a(δ)，δ = (m-1, …, 1, 0)。
- `principal_specialization(lam, m: int)` — s_λ(1, q, …, q^{m-1}) 當成 q 的多項式 {次數: 係數}（書 p.76 練習 5）。
- `character_table(n: int)` — χ^λ_μ（書 p.78 (11)：p_μ = Σ_λ χ^λ_μ s_λ）。回傳 (n 的分割, 表[λ][μ])，都用 partitions(n) 的順序。
- `schur_inner(f: Poly, g: Poly, m: int)` — 對稱函數的內積：以 s_λ 為正交規範基底（書 p.78；m 要 ≥ 次數）。
- `omega(f: Poly, m: int)` — 對合 ω：s_λ ↦ s_λ̃（書 p.78）。

## `combkit.rep_data`

rep_data.py —— Fulton 第 7 章（場景 F28–F41）的固定算例與小工具（純 Python）

常數：`_HERE`、`_DATA`、`_CACHE`

- `load(key: str)` — 
- `compose(p, q)` — (p∘q)(i) = p(q(i))：先 q 後 p。
- `sign(p)` — 
- `cycles(p, keep_fixed=True)` — 
- `cycle_type(p)` — 
- `transposition(a, b, n)` — 
- `act(p, rows)` — 
- `tabloid_key(rows)` — 
- `tabloid_str(key)` — 
- `columns(rows)` — 
- `row_group(rows)` — 
- `col_group(rows)` — 
- `polytabloid(rows)` — v_T = Σ_{q∈C(T)} sgn(q) {qT}，回傳 {tabloid_key: 係數}（略去 0）。
- `act_vec(p, vec)` — 
- `t21(i)` — 
- `specht21(vec)` — 零和向量 x1 t1 + x2 t2 + x3 t3 = a v_a + b v_b：a = x3、b = x2（且 x1 = -a - b）。
- `rep_matrix21(p)` — 
- `gr_mul(a, b)` — 
- `symmetrizers(rows)` — 
- `check()` — 

## `combkit.raney`

raney.py —— G. N. Raney, *Functional composition patterns and power series reversion*,

常數：`NOTE015_SEQ`

- `rank(s)` — r(σ) = Σ (σ(x) − 1)。
- `head_ranks(s)` — [r(空), r(σ 的前 1 項), …, r(σ)]，也就是 Raney 的「秩路徑」。
- `is_list(s, n: int)` — s 是否為 n 字列（Rosenbloom 判別法的推廣）：r = −n，真開頭的秩都 > −n。
- `is_word(s)` — 
- `is_list_recursive(s, n: int)` — 直接照定義：空字串是 0 字列；字 = 0，或 k α_1 … α_k（α_j 是字）；n 字列 = 字接 (n−1) 字列。
- `split_words(s)` — 把字列切成字：每次取最短的、秩為 −1 的開頭（Raney 定理 2.1 證明裡的做法）。
- `word_to_tree(s)` — 字 → 樹：(arity, (子樹, …))；葉子是 (0, ())。
- `tree_to_word(tree)` — 
- `composition(s, K: int | None=None)` — (n; a_1, …, a_K)：n = −r(σ)，a_i = 符號 i 的個數。
- `a0_of(n: int, a)` — 引數（符號 0）的個數 a_0 = n + Σ (i − 1) a_i。
- `length_of(n: int, a)` — 字串長度 m = n + Σ i a_i。
- `multinomial(*parts)` — M(a_0, a_1, …) = (Σ a_i)! / Π a_i!。
- `strings_with(n: int, a)` — 所有組成為 (n; a) 的字串（字典序）。
- `lists_with(n: int, a)` — 
- `L_formula(n: int, a)` — Raney 定理 2.2：m·L(n; a) = n·M(a_0, a_1, …)。
- `L_brute(n: int, a)` — 
- `cyclic_variant(s, k: int)` — σ g_k：(σ g_k)(x) = σ(x + k mod m)，也就是從第 k 項開始讀一圈。
- `list_shifts(s)` — 所有使 σ g_k 成為 n 字列的 k（n = −r(σ)）。
- `raney_alpha(s)` — Raney 證明定理 2.1 的第一步：在所有循環變形的所有開頭裡，取秩最大、同秩中最長的 α。
- `raney_proof_trace(s)` — 照 Raney 的證明走一遍：α、σ_1、σ_1 = α_1 … α_n、n 個循環 σ_i，以及它們對應的 k。
- `split_head(s, m: int)` — (m+n) 字列唯一的「m 字列開頭」：回傳 (前 m 個字, 其餘)。
- `sub_sequences(a)` — 
- `thm23_rhs(m: int, n: int, a)` — 
- `thm24_rhs(k: int, a)` — 
- `to_chu_path(s, v_vec)` — Raney 字串 → Chu 講義（counting.py 第 4 節）的路徑字母：倒過來讀，0 → LEAF(=0)，符號 v_i → i。
- `from_chu_path(path, v_vec)` — 
- `chu_heights(s)` — Chu 講義路徑的高度（倒過來讀、每步 1 − j）。等於「後綴的秩」取負號。
- `partition_to_word(blocks, N: int)` — Chu 的不交叉分法 → Raney 的字（「固定一個點」的分解，Chu 講義第 4 節）。
- `word_to_partition(s, N: int)` — partition_to_word 的反函數。
- `gessel_suffix_code(trees)` — Gessel §4.5 的 suffix code：c(T) = c(T_1)…c(T_j) j（先寫子樹、最後寫根的小孩數）。
- `prefix_code(trees)` — Raney 的前序碼（n 字列）。
- `a_from_chu(n_vec, v_vec)` — Chu 的 (n⃗, v⃗) → Raney 的 (a_1, …, a_K)，K = max v_i。
- `ps_mul(f, g, N: int)` — 
- `ps_pow(f, m: int, N: int)` — 
- `ps_compose(phi, w, N: int)` — φ(w)，w 的常數項為 0。
- `solve_w(phi, N: int)` — Raney 式 (5.3)：w = x·Σ a_n w^n 的唯一解（逐次代入）。
- `lagrange_coeff(phi, n: int, t: int=1)` — Lagrange（–Bürmann）：[x^n] w^t = (t/n) [x^{n−t}] φ^n（n ≥ 1）。t = 1 就是 Raney 式 (5.5)。
- `exp_series(c, N: int)` — e^{c x} 的係數。
- `lstar_coeffs(K: int, max_size: int)` — L*(p) = L(1; p)，p = (p_1, …, p_K)，只取 Σ i p_i ≤ max_size。
- `mv_mul(f: dict, g: dict, max_size: int)` — 
- `mv_pow(f: dict, n: int, K: int, max_size: int)` — 
- `a_table(phi, M: int, N: int)` — a(m, n) = [x^n] φ^m，0 ≤ m ≤ M、0 ≤ n ≤ N。
- `b_from_a(phi, m: int, n: int)` — Raney 式 (7.2)：b(m, n) = m/(m+n)·a(m+n, n)，b(0, 0) = 1。
- `c_from_b(phi, alpha, k: int)` — Raney 式 (7.5)：c(α, k) = Σ_m α^m b(m, k−m)/m!。
- `gould_G(alpha, n: int, beta, gamma)` — Raney 式 (7.7)：G(α, 0; β, γ) = 1，G(α, n; β, γ) = α/n!·Π_{m=1}^{n−1} (α + βn − γm)。
- `gould_phi(beta, gamma, N: int)` — Raney 式 (7.9)：φ(x) = γ x e^{βx} / (e^{γx} − 1)（γ = 0 時是 e^{βx}）。
- `trinomial_coeff(n: int, k: int)` — Raney 式 (6.4)：w = 1 + x w^n 的 x^k 係數 = 1/(1+(n−1)k)·C(nk, k)。
- `kelly_P(n: int, r: int)` — Problem 4983（Klamkin–Shepp，Kelly 解）：n 個因子、每次 r+1 個相乘的非結合乘積個數。
- `kelly_P_brute(n: int, r: int)` — 直接數：因子 a_1 … a_n，每一對括號恰好包 r+1 個項。
- `polygon_dissections_brute(sides: int, cell: int)` — 凸 sides 邊形用不交叉的對角線切成全是 cell 邊形的方法數（暴力：選對角線集合）。
- `heights(seq)` — [S_0, S_1, …, S_M]，S_0 = 0。
- `reading(seq, i: int)` — 從第 i 點開始讀一圈：x_{i+1}, …, x_{i+M}。
- `positive_count(seq)` — 出發後高度 > 0 的點數。
- `is_good_reading(seq)` — 
- `good_starts(seq)` — 所有好讀法的起點 i。
- `last_minimum(seq)` — S_1..S_M 中最後一個最小值的位置 x（1 起算）。從第 x 點開始讀就是好讀法（和為 1 時）。
- `positive_count_order(seq)` — 加強版的排序法：1..M 依 (S_k, −k) 由小到大。排第 r 位的起點有 M − r + 1 個正的點。
- `periodic_heights(seq, periods: int)` — 重複 periods 個週期的高度 [S_0, …, S_{periods·M}]。
- `stack_ops_to_steps(ops: str)` — 進棧 'P' → +1、出棧 'O' → −1。
- `stack_sequences(n: int)` — 所有合法的進出棧操作（n 次進、n 次出、出棧時棧不空），字典序。
- `pop_order(ops: str)` — 操作序列 → 出棧順序（進棧順序 1, 2, …）。
- `draws_left(cards)` — 摸牌遊戲：一開始能摸 1 張，摸到 w 元的牌變化 w − 1。回傳 [1, 摸完第 1 張後, …]。
- `deck_ok(cards)` — 前 m 張牌（勝利牌在最底下，不列入）能否摸到勝利牌：每張之後都還能再摸。
- `card_game_brute(w, m: int)` — 所有合格的牌序（牌有標號）：特殊牌 S1..Sn（w_i），無用牌 U1..U(m−n)。回傳標號 tuple 的串列。
- `card_game_formula(w, m: int)` — m!/(m − n + 1)，n = w_i ≥ 1 的張數。

## `combkit.tikz`

tikz.py —— 用 Python 產生 TikZ 程式碼（純 Python，不需要 manim）

常數：`COLOR`、`FILL`

- `f(x: float)` — 
- `pt(x: float, y: float)` — 
- `bbox(prims)` — 
- `emit(prims, dx: float=0.0, dy: float=0.0)` — 
- `place(prims, cx: float, cy: float)` — 把圖元的外框中心移到 (cx, cy)。
- `node(p, text, opts='', size=(0.12, 0.12))` — 
- `colored_math(chars, parts, wrap=lambda c: c)` — 把連續屬於同一部分的字元包成一組上色（數學模式）。
- `polygon_vertices(m: int, r: float)` — 底邊（第 0 與第 m-1 點）水平放在最下面；第 0 點在左下，順時針到右下的第 m-1 點。
- `tri_prims(m: int, diags, r: float, tree=None, diag_w=1.1, dot_r=0.03)` — 凸 m 邊形的三角剖分；給 tree 時把根三角形、左右子多邊形上色。
- `tree_prims(tree, gap: float, level: float, colored=True, labels=None, node_r=0.05, leaf_s=0.04, edge_w=0.8)` — 滿二元樹，樹根在 (0, 0)。葉子依中序等距排列，內部節點放在它所有葉子的正中間
- `dyck_prims(word: str, unit: float, parts=None, w=1.1, base=True)` — 
- `lattice_prims(word: str, n: int, unit: float, color='cA', colors=None, w=1.1, ny=None)` — 
- `circle_point(i: int, npts: int, r: float)` — 
- `matching_prims(pairs, npts: int, r: float, parts=None, root=False, w=1.1, labels=False)` — 圓周 npts 點（第 0 點在正上方、順時針）的配對；root=True 時第 0 點畫成橘色。
- `paren_math(tree, factors='abcd', colored=True)` — 
- `word_math(word: str, parts, kind: str)` — Dyck 字轉成括號字串或開票順序，逐字依根／左／右上色。
- `picture(lines, opts='')` — 

## `combkit.fulton_appendix`

fulton_appendix.py —— Fulton《Young Tableaux》附錄 A（組合變形）與附錄 B（代數簇的拓撲）的計算（純 Python，場景 F56–F65）

- `dual_word(w, N)` — 
- `evacuation_steps(T, N)` — [(ΔT, 空出的格子, 放進去的 x*), …]：每步拿掉左上角 x、把剩下的斜表滑成一般表，空出的格子放 N+1−x。
- `evacuation(T, N)` — 
- `column_insert(T, x)` — 行插入 x → T：在第一行擠掉最上面（最小）的 ≥ x 的數，被擠的數去下一行；沒有就放在該行底部。
- `column_insert_word(w, T=())` — 
- `matrix_to_pairs(A)` — 
- `lex_sort(u, v)` — 
- `antilex_sort(u, v)` — 
- `rsk_array(u, v)` — 
- `burge(A)` — 作法 (1b)：反字典序陣列 (u_1..u_r; v_1..v_r)，依序行插入 v_r, v_{r-1}, …, v_1，並依序把 u_r, …, u_1 放進新格子。
- `corners(A)` — §A.4.2 定理 1 的四個（弱序）配對，由轉換過的陣列用 R-S-K 重算；自我測試核對它們等於 (P,Q)、(P*,Q*)、(R,S*)、(R*,S)。
- `dual_rsk(A)` — 0-1 矩陣（§A.4.3）：字典序陣列，依序行插入 v_1, …, v_r，並把 u_i 放在新格子的共軛（轉置）位置。回傳 (P, Q)。
- `zero_one_matrices(rows, cols)` — 
- `gale_ryser(lam, mu)` — 練習 20：列和 λ、行和 μ 的 0-1 矩陣存在 ⇔ μ 被 λ 的共軛支配（μ ⊴ λ~）。
- `compact_skew(col_lengths)` — 回傳每一行的 (頂端列, 長度)：相鄰兩行左邊較長則頂端對齊，右邊較長則底端對齊；最後整體平移使最上面是第 0 列。
- `frank_tableau(T, col_lengths)` — 形狀 compact_skew(col_lengths)、rectify 成 T 的斜表（唯一）；回傳各行的內容（由上而下）。
- `keys(T)` — 左鑰匙 K_-(T)：長度 c 的行換成「最左行長度為 c 的 frank 斜表」的最左行；右鑰匙 K_+(T) 同理取最右行。
- `betti_from_cells(counts)` — counts[k] = 複維數 k 的胞腔個數 ⇒ b_{2k} = counts[k]、奇數 Betti 數為 0；χ = Σ counts。
- `chern_projective(n)` — 
- `plucker_degree(r, m)` — 
- `check()` — 

## `combkit.raney_examples`

raney_examples.py —— 用程式「製造」Raney (1960) 的例子（純 Python，不需要 manim）

常數：`OPS`、`ARGS`、`SECTIONS`

- `multiset(n: int, a)` — 組成 (n; a) 的所有符號：a_0 個 0、a_1 個 1、……
- `random_string(n: int, a, rng: random.Random)` — 
- `random_list(n: int, a, rng: random.Random, trace: bool=False)` — 均勻隨機的 n 字列：隨機字串 → 隨機挑一個好的起點轉過去。trace=True 時一併回傳 (字串, 起點)。
- `random_word(a, rng: random.Random, trace: bool=False)` — 
- `random_composition(rng: random.Random, max_size: int=7, max_arity: int=3)` — 隨機的 (a_1, …, a_K)：Σ i·a_i ≤ max_size，至少一個非零。
- `random_sum_one(M: int, rng: random.Random, lo: int=-4, hi: int=3)` — 總和為 1 的隨機整數序列（Raney 引理用）。
- `random_noncrossing(n_vec, v_vec, rng: random.Random)` — 均勻隨機的不交叉分法（Chu 的 (n⃗, v⃗)）：均勻隨機的字 → 「固定一個點」的一一對應。
- `word_to_expression(word)` — 把字寫成運算式：(前序寫法, 一般寫法)。一元用 f,g,h、二元 F,G,H,…、三元 T,U,V,…；引數 a,b,c,…
- `plain(s)` — 
- `ex_expression(rng)` — 
- `ex_cycle(rng)` — 
- `ex_count(rng)` — 
- `ex_chu(rng)` — 
- `ex_lagrange(rng)` — 
- `ex_convolution(rng)` — 
- `ex_sum_one(rng)` — 
- `ex_cards(rng)` — 
- `examples_markdown(seed: int=7, count: int=3)` — 

## `combkit.rmt_moments`

rmt_moments.py —— 隨機矩陣的動差與卡塔蘭數（純 Python，不需要 manim）

- `double_factorial_odd(k: int)` — (2k-1)!! = 1·3·5···(2k-1) = 2k 個點的配對數 = 標準高斯的 E X^{2k}。
- `pairing_counts(k: int)` — (全部配對數, 不交叉配對數) = ((2k-1)!!, C_k)，用暴力列舉。
- `goe_moment_poly(k: int)` — GOE 的 E Tr J^{2k}，表示成 {N 的冪次: 係數}。k=2：2N + 5 + 5/N；k=3：5N + 22 + 52/N + 41/N^2。
- `gue_moment_poly(k: int)` — GUE 的 E Tr H^{2k}（Harer–Zagier 展開）：N^{1-2g} 的係數是黏成虧格 g 曲面的方法數。k=2：2N + 1/N。
- `normalized_moment(poly: dict[int, int], N)` — (1/N) E Tr M^{2k} 在給定 N 的精確值（Fraction）。
- `tree_walks(k: int)` — 長度 2k 的「樹形閉路徑」（深度優先走一棵 k 條邊的平面樹再回到根），頂點按第一次到訪編號。共 C_k 條。
- `walk_to_dyck(walk: list[int])` — 樹形閉路徑 → Dyck 字：走到新頂點記 U，回到已到訪的頂點記 D。
- `walk_pairing(walk: list[int])` — 樹形閉路徑裡同一條邊的兩次經過配成一對（第 t 步是邊 walk[t]–walk[t+1]）。

## `combkit.schubert`

schubert.py —— Fulton《Young Tableaux》第 9–10 章（場景 F46–F55）的計算（純 Python，不依賴 manim）

- `length(w)` — 
- `inversions(w)` — 
- `descents(w)` — 
- `act_right(w, i)` — w·s_i：交換第 i、i+1 個位置的值（1 起算）。
- `w0(m)` — 
- `perms(m)` — 
- `lehmer_code(w)` — 
- `perm_str(w)` — 
- `rank_matrix(w)` — 
- `bruhat_leq(u, v)` — 
- `bruhat_leq_tableau(u, v)` — 練習 8：每個 p，把 u(1..p)、v(1..p) 各自排序後，u 的每一項 ≤ v 的對應項。
- `bruhat_chain(u, v)` — Lemma 11 的演算法：從 v 一步步往下（每步交換一對、長度減 1）；到得了 u 就回傳整條鏈，否則 None。
- `bruhat_levels(m)` — 
- `bruhat_covers(m)` — 
- `flag_cell_pattern(w)` — 第 p 列在第 w(p) 行放 1，1 的右邊是 0，每個 1 的下方是 0，其餘是 *（Fulton p.157）。
- `grass_cell_pattern(lam, r, n)` — Gr_r(C^{r+n}) 的胞腔 Ω°_λ（Fulton p.147）：第 i 列的 1 在第 n+i−λ_i 行。
- `q_factorial(m)` — 
- `gaussian_binomial(m, r)` — [m choose r]_q 的係數 = 各維數的 Schubert 胞腔個數（方框 r×(m−r) 裡 |λ| = k 的 λ 個數）。
- `count_subspaces_F2(m, r)` — F_2^m 裡 r 維子空間的個數（暴力：所有 r 個向量張成的空間去重）。
- `count_flags_F2(m)` — F_2^m 裡完全旗 0 ⊂ V_1 ⊂ … ⊂ V_{m-1} ⊂ F_2^m 的個數（暴力）。
- `plucker(A)` — 
- `plucker_from_coords(x, d, m)` — Fulton p.133 式 (2)：取 x_I ≠ 0 的 I（正規化成 1），a_{s,t} = x_{I 的第 s 個換成 t}。
- `flag_ring_dim(a, b)` — S(3;1,2) = Sym(E ⊕ ∧²E)/(關聯關係) 在雙次數 (a,b) 的維數；關係是非零因子，所以扣掉 (a−1,b−1) 的維數。
- `box_product(lam, mu, r, n)` — 
- `rectangle_count(r, m)` — 練習 17：與 r(m−r) 個一般的 (m−r) 維子空間都相交的 r 維子空間個數。
- `swap_vars(f, i)` — 
- `divided_difference(f, i)` — ∂_i：x_i^a x_{i+1}^b ↦ Σ_{k=0}^{a−b−1} x_i^{a−1−k} x_{i+1}^{b+k}（a > b），a < b 時取負號，a = b 時為 0。
- `schubert_poly(w)` — 𝔖_w，用 len(w) 個變數（Lemma 10 與命題 5：𝔖_{w·s_i} = ∂_i 𝔖_w 當 w(i) > w(i+1)）。
- `schubert_path(m, steps)` — 從 w0 出發依序做 ∂_i（i 在 steps 裡）；每步都要求目前 w(i) > w(i+1)。回傳 [(w, 𝔖_w), …]。
- `pipe_dreams(w)` — 縮減 pipe dream：階梯 {(i,j): i+j ≤ m} 中的十字集合 D，逐列由上而下、每列由右而左讀，
- `pipe_exits(D, m)` — 每條 pipe 從第 i 列的左邊進來、往右或往上走（十字直走；彎管：從左進→往上出，從下進→往右出），
- `pipe_weight(D, m)` — 
- `grassmannian_partition(w)` — w 只在 r 有下降 ⇒ λ = (w(r)−r, …, w(1)−1)（命題 8）；回傳 (λ, r)。
- `grassmannian_perm(lam, r, m)` — 
- `staircase_monomials(m)` — 
- `poly_str(f, var='x')` — {(2,1,0): 1, (1,1,1): 2} → 'x_1^{2}x_2+2x_1x_2x_3'（字典序由大到小）。
- `check()` — 

## `combkit.sym_abel`

sym_abel.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 2 章

常數：`I`、`R_`、`HBAR`、`ROT3`、`D2_NAMES`、`D2_TABLE_BOOK`、`C6_NAMES`、`C6_BOOK`、`NOT_GROUP_NAMES`、`NOT_GROUP_BOOK`、`C2_TABLE`、`C3_TABLE`、`GROUP_COUNTS`、`GROUP_COUNT_64_BOOK`、`ABELIAN_COUNTS`、`R2`、`PEND2`、`DETUNED`、`LISSA`、`PHASOR_NU`、`PHASOR_TIMES`、`C2_CHARS`、`O_R3`、`O_REP3`、`T3`、`PEND3`、`RING3`、`C3_STANDING`、`C3_MOVING`、`ALIAS`、`BROKEN`、`GALLOP`、`SWR_LIST`、`D2_TABLE`、`D2_REP4`、`_D2_BOOK_ACTION`、`D2_REP4_BOOK`、`D2_CHARS`、`D2_ABC`、`ABELIAN8`、`POINT3`、`C2H`、`C2V`、`C2V_BOOK`、`ROTOR_N`、`KP_W`、`KP_B`、`KP_CROSSINGS`、`NH3`、`NH3_NU`、`H_PLANCK_EV`、`LZ_BETAS`、`WELL`

- `to_np(M)` — 
- `table_from_rows(names, rows)` — 照抄的表 → {(g, h): gh}。rows[i][j] 是第 i 列（g = names[i]）、第 j 行（h = names[j]）的格子。
- `table_from_matrices(mats: dict)` — 由表示矩陣算乘法表：gh = 與 mats[g]·mats[h] 相等的那個名字。
- `table_rows(table: dict, names)` — 
- `assoc_failures(table: dict)` — 所有 (x, y, z) 使 (xy)z ≠ x(yz)。
- `axioms(table: dict)` — Table 2.2.1 的四條公理（外加「每列每行每個元素恰好一次」）。
- `is_abelian(table: dict)` — 
- `cyclic_names(n: int)` — 
- `cyclic_table(n: int)` — C_n = {1, r, …, r^{n−1}}：r^i r^j = r^{(i+j) mod n}（(2.2.12)、(2.6.2)）。
- `element_order(table: dict, g)` — 
- `order_counts(table: dict)` — 
- `abelian_count(n: int)` — 有限 Abel 群的個數 = Π_p p(e_p)（n = Π p^{e_p}，p(e) 是 e 的分割數）——§2.10 的定理。
- `abelian_types(n: int)` — n 階 Abel 群的全部結構（質數冪循環群的外積，大到小）。
- `accel2(a=PEND2['a'], b=PEND2['b'])` — (2.2.1c)：[[a+b, −a], [−a, a+b]]。
- `commutator(A, B)` — 
- `c2_projectors(R=R2)` — 最小方程 R² − 1 = 0 的 Hamilton–Cayley 投影算符（第 1 章的同一招，sym_mat.projectors）：{1: 𝒫⁺, −1: 𝒫⁻}。
- `c2_projectors_char(R=R2)` — 由特徵標：𝒫^± = (1 ± R)/2（(2.3.2)_x）。
- `normal2(a=PEND2['a'], b=PEND2['b'])` — (α₊, α₋) = (b, 2a + b)（p. 71）。
- `beat_times2(a=PEND2['a'], b=PEND2['b'])` — 完全轉移 π/Δω、拍的週期 2π/Δω（精確 sympy）。
- **class `Modal`** — ẍ = −𝒜x 的解析解：x(t) = Σ_j e_j (c_j cos ω_j t + s_j sin ω_j t/ω_j)（ω_j = 0 時 c_j + s_j t）。給 updater 用（純 numpy）。
  - `.velocity(t)` — 
- `ode_solve(A, x0, v0, ts)` — 數值積分 ẍ = −𝒜x（solve_ivp, DOP853）→ 陣列 (len(ts), n)。
- `detuned_matrix(a=DETUNED[0], b1=DETUNED[1], b2=DETUNED[2])` — 
- `detuned_max_transfer(a=DETUNED[0], b1=DETUNED[1], b2=DETUNED[2])` — x(0) = (1, 0)：x₂(t) = (sin 2θ/2)(cos ω₁t − cos ω₂t)，最大振幅 sin 2θ，tan 2θ = 2a/(b₂ − b₁)。
- `a_for_ratio(rho, b=1)` — ω₋/ω₊ = ρ ⇔ 2a + b = ρ²b ⇒ a = (ρ² − 1)b/2。
- `lissajous(t, rho)` — 從 x = (1, 0) 靜止放手：x₁ = (cos t + cos ρt)/2、x₂ = (cos t − cos ρt)/2（ω₊ = 1）。
- `lissajous_numeric(ts, rho)` — 
- `closure(rho)` — (2.3.11)–(2.3.13)：ω₊/ω₋ = n₊/n₋（互質）⇒ t_LJ = 2πn₊/ω₊ = 2πn₋/ω₋ = (n₋ − n₊)·t_beat（拍數 n₋ − n₊）。
- `quad_potential(a=PEND2['a'], b=PEND2['b'])` — (2.3.9) 的力常數矩陣 𝓕（m = ℓ = 1 時 = 加速度矩陣）；V = ½ xᵀ𝓕x。
- `ellipse_axes(a, b=1, V0=R_(1, 2))` — V = V0 的等高線：沿 (1,1) 的半軸 √(2V0/α₊)、沿 (1,−1) 的半軸 √(2V0/α₋)。
- `phasors(t, nu=PHASOR_NU)` — 𝒞 = Ce^{−iωt}：z₁ = (e^{−iω₊t} + e^{−iω₋t})/2、z₂ = (e^{−iω₊t} − e^{−iω₋t})/2（t = 0：等量的 (+)、(−)，x = (1, 0)）。
- `rep_projectors(reps: dict, chars: dict)` — 𝒪(𝒫^α) = (1/°G) Σ_g D^α(g)* 𝒪(g)。
- `multiplicities(reps: dict, chars: dict)` — 重複頻率 f^α = (1/°G) Σ_g D^α(g)* χ(g)，χ(g) = tr 𝒪(g)。
- `accel3(m=1, M=2, k=1, l=1, g=1)` — (2.5.1c)：a = k/mℓ²、A = k/Mℓ²、b = g/ℓ。
- `reduce3(**kw)` — (2.5.3)：T3†𝒜T3（部分約化：2×2 ⊕ 1×1）。
- `omega2_3(**kw)` — 三擺的 ω²（精確，由小到大）。
- `eps(n: int, k: int=1)` — e^{2πik/n}（精確，展成 cos + i sin）。
- `cyclic_rep(n: int)` — 𝓡(r)：r|p⟩ = |p+1⟩（(2.6.7)：𝓡(r)x = (x_n, x_1, x_2, …)）。
- `cn_char_table(n: int)` — D^{k_n}(r^p) = e^{2πikp/n}：列 k = 0..n−1、行 p = 0..n−1（(2.6.5)、Fig. 2.7.2）。
- `cn_projectors_char(n: int)` — (2.6.6b)：𝒫^{k_n} = (1/n) Σ_m e^{−2πikm/n} r^m。
- `cn_projectors_minimal(n: int)` — (2.3.2)：𝒫^{ε_k} = Π_{l≠k}(r − ε_l)/Π_{l≠k}(ε_k − ε_l)，r^n − 1 = Π(r − ε_l)。
- `cn_projectors_np(n: int, method: str='char')` — 浮點版（大 n 用）：method = "char"（(2.6.6b)）或 "minimal"（(2.3.2) 的連乘）。
- `ring_matrix(n: int, a=1, b=0)` — (2.6.1)、(2.7.1)：對角 2a + b、最近鄰 −a（環）。n = 2 時兩根彈簧並聯：[[2a+b, −2a], [−2a, 2a+b]]。
- `ring_omega2(n: int, a=1, b=0)` — (2.7.4)：ω²(k_n) = b + 2a(1 − cos 2πk/n)，k = 0..n−1（精確）。
- `wave_vector(n: int, k: int)` — (2.7.3)：|k_n⟩ 的分量 e^{−ik_n p}/√n。
- `cos_wave(n: int, k: int)` — (2.7.5a)：√(2/n) cos(k_n p)（k = 0 或 n/2 時再除以 √2 才歸一）。
- `sin_wave(n: int, k: int)` — (2.7.5b)：√(2/n) sin(k_n p)。
- `dispersion(k, a=1, b=0)` — (2.7.6)：ω(k) = √(2a(1 − cos k) + b)。
- `dispersion_series(k, a=1)` — 
- `brillouin(n: int)` — 第一 Brillouin 區的 n 個整數 k：−n/2 < k ≤ n/2（n 偶）、−(n−1)/2 ≤ k ≤ (n−1)/2（n 奇）。
- `alias(n: int, k: int)` — 把 k 換成第一 Brillouin 區裡等價的整數（(2.7.7)）。
- `alias_values(n: int=5, k: int=2)` — 
- `broken_ring_matrix(n=24, abar=R_(3, 2), alow=R_(1, 2), b=0)` — (2.7.8)：彈簧 p ↔ p+1 的常數 a̲（p 偶）或 ā（p 奇）。
- `broken_block(k, abar=R_(3, 2), alow=R_(1, 2), b=0)` — (2.7.9)（對角加回 b）：基底 |k⟩（偶數格點）、|k'⟩（奇數格點）。
- `broken_bands(k, abar=R_(3, 2), alow=R_(1, 2), b=0)` — (2.7.10)（加回 b）：ω² = b + ā + a̲ ∓ √(ā² + a̲² + 2āa̲ cos k)。回傳（聲學, 光學）。
- `broken_standing(n: int=24)` — k = π（C_{n/2} 的 6₁₂）的兩個駐波：(1, −1, −1, 1, …) 只拉軟彈簧（ω² = 2a̲ + b）、(1, 1, −1, −1, …) 只拉硬彈簧（2ā + b）；
- `swr(I_, R_amp)` — (2.7.13)：Δ = (I − R)/(I + R) = |ψ| 的極小/極大。
- `gallop_psi(x, t, I_, R_amp, k=1.0, w=1.0)` — (2.7.12)：ψ = (Ie^{ikx} + Re^{−ikx})e^{−iωt}。
- `gallop_envelope(x, I_, R_amp, k=1.0)` — (2.7.14)：|ψ| = [I² + R² + 2IR cos 2kx]^{1/2}。
- `gallop_zeros(t, delta, k=1.0, w=1.0, span=2 * math.pi)` — Re ψ = 0 的位置（I = (1+Δ)/2、R = (1−Δ)/2）：tan kx = −cot(ωt)/Δ；在 [0, span) 內。
- `zero_velocity(t, delta, c=1.0, w=1.0)` — 零點的速度 v = cΔ/(Δ² sin²ωt + cos²ωt)：在 cΔ（ωt = 0）與 c/Δ（ωt = π/2）之間（習題 2.7.2）。
- `regular_rep(table: dict, names)` — 正則表示：g|h⟩ = |gh⟩（基底順序 names）。
- `d2_projectors_product()` — (2.8.3)–(2.8.5b)：𝒫¹ = 𝒫^{x+}𝒫^{y+}、𝒫² = 𝒫^{x−}𝒫^{y+}、𝒫³ = 𝒫^{x+}𝒫^{y−}、𝒫⁴ = 𝒫^{x−}𝒫^{y−}（以 A1、B2、B1、A2 命名）。
- `d2_projectors_char()` — 
- `d2_matrix(A=D2_ABC[0], a=D2_ABC[1], b=D2_ABC[2], c=D2_ABC[3])` — (2.8.1a)（ẍ = −𝒜x）：[[A, a, b, c], [a, A, c, b], [b, c, A, a], [c, b, a, A]]。
- `d2_omega2(A=D2_ABC[0], a=D2_ABC[1], b=D2_ABC[2], c=D2_ABC[3])` — Fig. 2.8.2：ω²(α) = Σ_g D^α(g)·(第一列的係數)。
- `abelian_group(orders)` — C_{n1} × C_{n2} × …：元素是指數組 (e1, e2, …)，乘法＝逐分量加（mod n_i）。回傳 (elements, table)。
- `abelian_char_table(orders)` — 外積的特徵標表：χ_k(g) = Π_i ε_{n_i}^{k_i g_i}。回傳 (elements, irreps, sympy Matrix[列 irrep, 行 元素])。
- `outer_char_table(T1: sp.Matrix, T2: sp.Matrix)` — (2.10)：H × K 的特徵標表 = 兩張表的 Kronecker 積（列 (α, β)、行 (h, k)）。
- `is_cyclic(orders)` — 有一個元素的階 = 群的階 ⇔ 循環（⇔ 各 n_i 兩兩互質）。
- `is_group(names)` — 
- `generated_order(M: sp.Matrix)` — 
- `point_table(names)` — 點群的乘法表（由 POINT3 的矩陣）。
- `free_levels(mmax: int=12)` — (2.12.5)：ε_m = m²E。
- `rotor_matrix(V=WEAK_V, mmax: int=60, n: int=ROTOR_N)` — (2.12.12)：⟨m'|H|m⟩ = δ m² + G(m' − m)，G(±n) = −V/2。回傳 (ms, H)。
- `rotor_block(kappa: int, V=WEAK_V, N: int=8, n: int=ROTOR_N)` — (2.12.16)：只含 m = κ + nN'（|N'| ≤ N）的子矩陣（C_n 的 irrep κ）。回傳 (ms, H)。
- `rotor_bands(V=WEAK_V, nb: int=3, n: int=ROTOR_N, N: int=8)` — {κ: 最低 nb 個本徵值}，κ 走第一 Brillouin 區（−5..6）。
- `pair_exact(m: int=5, G=-1, n: int=ROTOR_N)` — (2.12.17)：|m⟩、|m−n⟩ 的 2×2：(m² + (m−n)²)/2 ∓ √(((m−n)² − m²)²/4 + G²)。
- `pair_pert(m: int=5, G=-1, n: int=ROTOR_N)` — (2.12.17) 的近似：m² − G²/((m−n)² − m²)、(m−n)² + G²/(…)。
- `tb_matrix(n: int, H=sp.Symbol('H'), S=sp.Symbol('S'))` — (2.12.23) 第二式：對角 H、最近鄰 −S（環）。
- `tb_levels(n: int, H=sp.Symbol('H'), S=sp.Symbol('S'))` — (2.12.25)：ε_k = H − 2S cos(2πk/n)，k = 0..n−1。
- `tight_fit(V=STRONG_V, n: int=ROTOR_N)` — 最低帶的 n 個能階（每個 irrep κ 一個）→ H = 平均、S = (ε(π) − ε(0))/4、與 H − 2S cos k 的最大偏差。
- `mathieu_edges(Vp: float, nb: int=3, N: int=30)` — 縮放的 Mathieu 方程 ψ'' + (ε + V' cos 2t)ψ = 0 的帶邊：A（週期 π，基底 e^{2ijt}）與 B（反週期，e^{i(2j+1)t}），
- `mathieu_edges_scipy(Vp: float, nb: int=3)` — 同上，用 scipy.special.mathieu_a/b（q = V'/2，帶 j = [a_j, b_{j+1}]）。
- `rotor_edges(V: float, nb: int=3, n: int=ROTOR_N)` — C_n 轉子（E = 1）的帶邊 = (n/2)²·Mathieu 帶邊（V' = V/(n/2)²）：n = 12 時乘 36（Fig. 2.12.5 的縱軸）。
- `square_well_G(m: int, U=1, n: int=ROTOR_N)` — n 方井位能 V(φ) = U(nφ)（|θ| < π/2 為 0、其餘 U）的 Fourier 係數 G(m) = (1/2π)∫e^{−imφ}V dφ（精確）。
- `square_well_G_fft(m: int, U=1.0, n: int=ROTOR_N, M: int=4096 * 3)` — 
- `seg_matrix(E: float, U: float, L: float)` — 常數位能區段的 (ψ, ψ') 傳遞矩陣（ψ'' = (U − E)ψ）。
- `kp_half_trace(E: float, U: float, w: float=KP_W, b: float=KP_B)` — ½ tr(M_位障 M_井) = cos(qL)（Kronig–Penney 方程；|·| ≤ 1 是能帶）。
- `kp_half_trace_ode(E: float, U: float, w: float=KP_W, b: float=KP_B)` — 同上，用 solve_ivp 積分一個週期（單擺式的 monodromy）。
- `kp_edges(U: float, Emax: float=60.0, n: int=6000)` — ½ tr = ±1 的所有根（E < Emax），由小到大。
- `kp_gap_near(U: float, E0: float, half: float=1.5)` — E0 附近最窄的「帶外」區間寬度（偶然簡併處 → 0）：在 [E0 − half, E0 + half] 內 |½tr| > 1 的最長連續段長度。
- `two_level(pE, H=0, S=1)` — (2.12.34)：ε = H ∓ √(S² + (pE)²)。
- `two_level_matrix(pE, H=0, S=1)` — (2.12.33)：[[H − pE, −S], [−S, H + pE]]（基底 up、dn）。
- `tunnel_populations(t, S=1)` — 從 |up⟩ 出發：P_up = cos²(St/ħ)、P_dn = sin²(St/ħ)（拍頻 2S/ħ）。
- `tunnel_numeric(ts, H=0, S=1)` — |ψ(t)⟩ = expm(−i𝓗t)|up⟩ 的 P_dn。
- `nh3_numbers()` — 
- `lz_formula(beta: float, S: float=1.0)` — Landau–Zener（書外）：pE = βt，對角差 2βt ⇒ 跳到另一支的機率 e^{−πS²/β}。
- `lz_numeric(beta: float, S: float=1.0, T: float | None=None)` — 數值積分 iψ' = 𝓗(t)ψ，pE = βt，t ∈ [−T, T]；從下面那支的本徵態出發，最後在上面那支的機率。
- `wavenumbers(e, V=WELL['V'])` — (2.12.37)（更正為 /ħ）：l = √(2mE)/ħ = √(e + V)（井內）、k = √(2m(E − V))/ħ = √e（井外）。
- `s_matrix(e, V=WELL['V'], a=WELL['a'], book: bool=False)` — (2.12.38)：S = [[i(l² − k²) sin 2la, 2lk], [2lk, i(l² − k²) sin 2la]]/D，
- `s_matrix_direct(e, V=WELL['V'], a=WELL['a'])` — 直接解 x = ±a 的四個邊界條件（ψ₁ = I₁e^{ikx} + O₁e^{−ikx}、ψ_w = Ae^{ilx} + Be^{−ilx}、ψ₂ = I₂e^{−ikx} + O₂e^{ikx}）。
- `phase_shifts(e, V=WELL['V'], a=WELL['a'])` — 偶：井內 cos lx、井外 cos(kx + δ₊) ⇒ tan(ka + δ₊) = (l/k) tan la；奇：sin lx、−cos(kx + δ₋) ⇒ tan(ka + δ₋) = −(l/k) cot la。
- `transmission(e, V=WELL['V'], a=WELL['a'])` — |S₂₁|² = 1/(1 + V² sin²(2la)/(4e(e + V)))。
- `resonances(V=WELL['V'], a=WELL['a'], emax: float=70)` — 2la = nπ ⇒ e = (nπ/2a)² − V（T = 1，Ramsauer–Townsend）。
- `crossing_step(g1, g2, x)` — (2.12.50b)：(R₁, L₁) = C (R₂, L₂)，ψ_j = R_j e^{g_j x} + L_j e^{−g_j x}，g_j = [2m(V_j − E)/ħ²]^{1/2}（更正 ħ²）。
- `crossing_well(e, V=WELL['V'], a=WELL['a'])` — (2.12.51)：(I₁, O₁) = C (O₂, I₂)。
- `crossing_well_steps(e, V=WELL['V'], a=WELL['a'])` — 兩個 (2.12.50b) 的乘積：外 → 井（x = −a）、井 → 外（x = a）；g = ik（外）、il（井）。
- `s_from_c(C)` — (2.12.52a)（更正 S₁₁ = C₂₁/C₁₁）。
- `c_from_s(S)` — (2.12.52b)。
- `bound_geometric(V=25.0, w=2.0)` — 習題 2.12.2（ħ = 1、m = ½）：井內 sin(kx + δ)，sin δ = k/√V，k = (nπ − 2δ)/w。回傳 [(k, δ, E = k²)]。
- `bound_fd(V=25.0, w=2.0, pad: float=10.0, N: int=4000)` — 有限差分（Richardson 外插）：−ψ'' + U ψ = Eψ，井 [0, w] 內 U = 0、外面 U = V；回傳 E < V 的本徵值。
- `c_barrier(theta)` — (2.12.54a)。
- `c_valley(alpha)` — (2.12.55a)。
- `wkb_cluster(Ej, hw, theta, n, j=0)` — (2.12.60)：E = E_j − (−1)^j (ħω/π)[1 + θ²]^{−1/2} cos(2πm/n)，m = 0..n−1。
- `wkb_vs_exact(Vp: float=10.0)` — Mathieu（縮放）位能 U(t) = −V' cos 2t（ħ = 1、2m = 1）的最低帶：
- `problem_233(A=1.0, B=0.3, D=1.0)` — H = A/2(p₁² + q₁²) + D/2(p₂² + q₂²) + B(q₁q₂ + p₁p₂)（C = 0）：𝒜 = [[A² + B², B(A + D)], [B(A + D), D² + B²]]；
- `response_234(w, a=3.0, b=1.0, d=0.1, e=0.1, f=(1.0, 0.0))` — 習題 2.3.4 的穩態複振幅 x = (𝒂 − ω² + iω𝒅)⁻¹ f（x(t) = Re[x e^{iωt}]）。
- `response_234_modal(w, a=3.0, b=1.0, d=0.1, e=0.1, f=(1.0, 0.0))` — 同上，用 C₂ 模：x± = f±/(a ± b − ω² + iω(d ± e))，f± = (f₁ ± f₂)/√2。
- `problem_261()` — C₂：[[a, b], [c, d]] 與 R 對易 ⇔ a = d、b = c；C₃：與 r 對易 ⇔ 循環矩陣 [[A, B, C], [C, A, B], [B, C, A]]
- `problem_262(m1=10, m2=1, k1=20, k12=1)` — 座標換尺度 y_i = √m_i x_i 後要求對角相等：(k₁ + k₁₂)/m₁ = (k₂ + k₁₂)/m₂ ⇒ k₂ = (k₁ + k₁₂)m₂/m₁ − k₁₂。
- `problem_271(n: int, a=1, b=0)` — 習題 2.7.1(e)：ω_k = √(b + 2a(1 − cos 2πk/n))，k = 0..n−1。n = 4：0、√2、√2、2；n = 5：0、1.1756、1.9021（各二重）；
- `problem_273(abar=R_(6, 5), alow=R_(4, 5))` — C₈ → C₄（ā = 1.2、a̲ = 0.8、b = 0）：k = 0：ω² = 0、4；k = ±π/2：2 ∓ √2.08 = 0.5578、3.4422（二重）；k = π：2a̲ = 1.6、2ā = 2.4。
- `problem_274(m: int=4, a=1.0, b=0.0)` — 兩牆之間 m 個擺 ＝ C_{2(m+1)} 環中「牆是節點」的那些正弦駐波：ω_j² = b + 4a sin²(jπ/(2(m+1)))，j = 1..m。
- `problem_275(a=-3.0, s=1.2, t=0.8)` — ẍ = 𝓜x，𝓜 = [[a, s, t], [t, a, s], [s, t, a]]（循環 ⇒ 有 C₃ 對稱，但沒有反射對稱）：
- `problem_281(j=1.0, k=2.0, alpha_deg=45.0, m=1.0)` — 「Double Trouble」（模型假設：上下各一條繩子繞過定滑輪，兩根 j 串聯成 j/2）：
- `problem_282()` — n = 8：3（C8、C4×C2、C2×C2×C2）；9：2；12：2；16：5；24：3。
- `problem_283()` — 12 質量環的方形駐波 s_p = −S(2π(p + ½)/12)（前 3 個在下、中間 6 個在上、後 3 個在下）：
- `problem_2121(Emax: float=75.0)` — Kronig–Penney（a = b = π/2）的偶然簡併：E = (2n_w)²、U = E − (2n_p)²（1 ≤ n_p < n_w）。
- `problem_2122(V: float=25.0, w: float=2.0)` — 
- `line_nodes(E, U, nw=4, w=KP_W, b=KP_B)` — 直線上 nw 個方井（外面無限長的 U）：從左邊衰減的解的節點數 ＝ 低於 E 的本徵值個數（Sturm 振盪定理）。
- `line_levels(U: float=60.0, nw: int=4, count: int=16)` — 二分法：第 i 個本徵值 = 使 line_nodes(E) 由 i 變 i+1 的 E（到 1e-12）。
- `problem_2123(U: float=60.0, Emax: float=80.0)` — (a) 圓上 4 個方井（C₄ 環，週期 π）：½tr = cos(2πm/4) = 1、0、−1；(b) 直線上 4 個方井（外面是無限長的 U）：節點計數＋二分。
- `ring4_fd(U: float=60.0, N: int=4000)` — 問題 2.12.3(a) 的有限差分（週期邊界，Richardson）：最低 8 個本徵值。
- `line4_fd(U: float=60.0, pad: float=2.5, N: int=6000)` — 
- `summary()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上算法核對。

## `combkit.sym_coup`

sym_coup.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 7 章

常數：`R_`、`IU`、`HALF`、`SPEC`、`_FACT`、`CG_717_BOOK`、`CG_7119_BOOK`、`HF_C716`、`HF_CONST`、`SLATER_P2`、`ALPHA`、`HARTREE_EV`、`HARTREE_CM`、`CM_TO_GHZ`、`T727_COLS`、`T727_ROWS`、`T727_BOOK`、`BOOK_7292`、`_SAMPLE_J`、`T733_BOOK`、`PAULI_S`、`EXPECT_CARD`、`T7312_BOOK`、`V_TABLES_BOOK`、`V74_DP_BOOK`、`_X40`、`_X44`、`CRYSTAL_COMBO`、`_O_REPS`、`FIG742`、`_D2_FIG743`、`CLUSTER_PATTERNS`、`_DIMS`、`HAM_CARD`、`SIG_C`、`U_RL`、`ELLIPSE_CARD`、`EULER_GAMMA`、`CARD_EXPECTED`

- `S_(x)` — int／float／Fraction／str／sympy → sympy Rational（半整數）。
- `ms(j)` — m = j, j−1, …, −j（sympy Rational）。
- `cg(j1, j2, j, m1, m2, m)` — C^{j1 j2 j}_{m1 m2 m}（sympy 精確，快取；Condon–Shortley）。
- `cg_racah(j1, j2, j3, m1, m2, m3)` — 書上的 Racah 式 (7.2.57)（Schwinger 振子推出；照印，負階乘的項略去）。
- `cg_float(j1, j2, j3, m1, m2, m3)` — 標準 Racah 式（Python 大整數＋Fraction，最後才開根號）；大 J（J = 30 的張量）用。
- `cg_lowering(j1, j2)` — (7.1.8)–(7.1.11)：從 |j1 j1⟩|j2 j2⟩ 用 J− = J−⊗1 + 1⊗J− 降下來；較低的 j 先在 M = j 子空間與已得的態正交化
- `cg_table(j1, j2, method: str='sympy')` — 耦合係數表：(Matrix, rows=[(m1, m2)], cols=[(j, m)])；method = "sympy"、"racah"（(7.2.57)）、"lowering"。
- `stretched_cg(j1, j2, m1, m2, book: bool=False)` — (7.1.17)：J = j1 + j2 的耦合係數（更正：分母 (j2+m2)!）；book=True 照印 (j1+m2)!（j1 − j2 非整數時沒有意義）。
- `d_half_product(beta=None)` — (7.1.6a)：C†(𝒟^½(0β0) ⊗ 𝒟^½(0β0))C（sympy；beta=None 用符號 β）→ 應為 𝒟¹(0β0) ⊕ [1]。
- `hyperfine_matrix()` — J^proton·J^electron 在 |m_p⟩|m_e⟩（(½½),(½,−½),(−½,½),(−½,−½)）的 4×4（精確）。
- `hyperfine_levels()` — (7.1.13)–(7.1.14)：a J^p·J^e 的本徵值（單位 a）→ {1/4: 3, −3/4: 1}。
- `hyperfine_states()` — 三重態 |1 1⟩, |1 0⟩, |1 −1⟩ 與單重態 |0 0⟩（(7.1.10ₓ)、(7.1.11)）＝ (7.1.7) 的四行。
- `hyperfine_numbers()` — (7.1.12)、(7.1.15) 與 21 cm 線：a = (8π/3)g_eβ_e g_pβ_p|ψ(0)|²，|ψ_1s(0)|² = 1/πa0³。
- `two_electron_terms(l1: int=1, l2: int=1, equivalent: bool=True)` — 數值：(反對稱化的) 兩電子空間裡 L²、S² 同時對角化 → {"^{2S+1}L": 態數}（(7.1.20)–(7.1.21)）。
- `terms_by_symmetry(l1: int=1, l2: int=1, equivalent: bool=True)` — (7.2.34)：C^{l l L}_{m2 m1 M} = (−1)^{2l−L}C^{l l L}_{m1 m2 M} ⇒ 軌道偶 L 對稱、奇 L 反對稱；自旋 S = 1 對稱、S = 0 反對稱；
- `hund_order(terms)` — Hund 規則：先比 S（大者低）、再比 L（大者低）→ 由低到高。
- `spin_orbit(l: int, book: bool=False)` — ⟨S·l⟩ 的本徵值（數值對角化 S·l 於 (2l+1)·2 維）→ {j: (值, 簡併)}；book=True 照印 (7.1.26)（沒有 ½）。
- `dirac_level(n: int, j, Z: int=1)` — (7.1.27)：ε_{n,j} = −Z²/2n² − (α²Z⁴/2n³)(1/(j+½) − 3/4n)（hartree）。
- `a_so(n: int, l: int, Z: int=1)` — (7.1.25b)：a_so = Z⁴α²/(2n³ l(l+½)(l+1))（hartree）。
- `fine_structure()` — 1s 的相對論位移、2p_{3/2} − 2p_{1/2}：Dirac (7.1.27) ＝ a_so·(⟨S·l⟩ 的差)（更正的 (7.1.26)）。
- `couple(j1, j2)` — 三角條件：j1 ⊗ j2 = |j1−j2| ⊕ … ⊕ (j1+j2)（Fraction）。
- `ls_jj_counts()` — (2p)² 的 Pauli 允許能階：LS（¹S0、³P0,1,2、¹D2）與 jj（(3/2)²：0、2；(3/2,1/2)：1、2；(1/2)²：0）的 J 值；
- `p2_levels(x: float, F2: float=1.0, zeta: float=4.0)` — 本系列的 (2p)² LS→jj 相關（Fig. 7.1.4 的數值版）：H = (1−x)H_C + xζ(l1·s1 + l2·s2)，H_C 依 (L, S) 取 SLATER_P2×F2。
- `cg_curve(j1, j2, J, M)` — Fig. 7.1.7：C^{j1 j2 J}_{m, M−m, M} 對 m（只取 |M − m| ≤ j2）→ (ms, 值)。
- `classical_range(j1, j2, J, M)` — 向量錐（Fig. 7.1.6）：|J| = √(J(J+1)) 與 z 軸夾 Θ（cos Θ = M/|J|），j1 繞 J 轉、與 J 夾 θ13（餘弦定律）
- `sign_changes(vals, tol: float=1e-12)` — 
- `spherical(v)` — (7.2.3)：A¹₁ = −(A_x + iA_y)/√2、A¹₀ = A_z、A¹₋₁ = (A_x − iA_y)/√2。
- `coupled_product(A, B, k: int, q: int)` — [A¹ × B¹]^k_q = Σ C^{1 1 k}_{m1 m2 q} A¹_{m1} B¹_{m2}（(7.2.4a)、(7.2.5a)、(7.2.6a)）。
- `table_727()` — (7.2.7)：列 A_iB_j（xx, xy, …, zz），行 [ ]^k_q；元素 = [A¹×B¹]^k_q 中 A_iB_j 的係數。
- `check_724()` — (7.2.4b)、(7.2.5b)、(7.2.6b)：[A×B]⁰₀ = −A·B/√3、[A×B]¹₀ = i(A×B)_z/√2、[A×B]²₀ = (−A_xB_x − A_yB_y + 2A_zB_z)/√6。
- `multipole_X(k: int, q: int)` — X^k_q = r^k √(4π/(2k+1)) Y^k_q 寫成 x, y, z 的多項式（Condon–Shortley 號）：
- `rr_tensor()` — (7.2.9)：[r¹r¹]²_q（A = B = r）與 X²_q；比值都是 √(2/3)。BOOK_7292 照印的 [rr]²₁ = −(x − iy)z。
- `ck_formula(k: int, q: int)` — (7.2.11)：C^{k−1 1 k}_{q−1 1 q}、C^{k−1 1 k}_{q 0 q}、C^{k−1 1 k}_{q+1 −1 q}。
- `scalar_cg(j, m)` — (7.2.12)：C^{j j 0}_{m, −m, 0} = (−1)^{j−m}/√(2j+1)。
- `threej(j1, j2, j3, m1, m2, m3)` — (7.2.30a)：(j1 j2 j3; m1 m2 m3) = (−1)^{j1−j2−m3}C^{j1 j2 j3}_{m1 m2 −m3}/√(2j3+1)。
- `threej_symmetries()` — (7.2.31)–(7.2.35) 在 _SAMPLE_J 的所有 m 上逐一核對；也與 sympy 的 wigner_3j 比較。
- `orthogonality_check(pairs=((1, HALF), (1, 1), (R_(3, 2), 1)))` — (7.2.21)、(7.2.22)：耦合係數表是正交矩陣（C^T C = 1 = C C^T）。
- `factorization_check(cases=((1, 1, 2, (1, 0, 1), (0, 1, 1)), (1, HALF, R_(3, 2), (0, HALF, HALF), (1, -HALF, HALF)), (2, 1, 2, (1, -1, 0), (2, 0, 2)), (1, 1, 1, (1, 0, 1), (1, -1, 0))))` — (7.2.25)：∫d(αβγ) 𝒟^{j3*}_{m3 m3'} 𝒟^{j1}_{m1 m1'} 𝒟^{j2}_{m2 m2'} = C^{j1 j2 j3}_{m1 m2 m3}C^{j1 j2 j3}_{m1' m2' m3'}/(2j3+1)
- `su2_irrep_formula(j, u)` — (7.2.38)：𝒟^j_{mn}(u)（u ∈ SU(2) 的 2×2）。
- `boson_coupled(j, m, n, N, half: bool=True)` — (7.2.54)：|j m n N⟩（Fock 佔據數 (n¹₁ n¹₂, n²₁ n²₂) → 係數）；half=False 照印（方括號沒有 ½ 次方）。
- `boson_cg(j1, j2, j3, m1, m2, m3, half: bool=True)` — C^{j1 j2 j3}_{m1 m2 m3} ＝ ⟨(7.2.36)|(7.2.54)⟩（n = j1 − j2、N = 2(j1 + j2)）。
- `boson_norm(j, m, n, N, half: bool=True)` — 
- `det_scalar_state(d: int)` — (7.2.50)：[(2d)!(2d+1)!]^{−½}(a¹₁a²₂ − a¹₂a²₁)^{2d}|00,00⟩ → {佔據數: 係數}；結果 = (2d+1)^{−½}Σ(−1)^r|r 2d−r, 2d−r r⟩。
- `sixj(a, b, c, d, e, f)` — 
- `sixj_from_3j(j1, j2, j3, j4, j5, j6)` — {j1 j2 j3; j4 j5 j6} = Σ (−1)^{Σ(j−m)} (j1 j2 j3; −m1 −m2 −m3)(j1 j5 j6; m1 −m5 m6)(j4 j2 j6; m4 m2 −m6)(j4 j5 j3; −m4 m5 m3)。
- `sixj_from_cg(j1, j2, j3, J, j12, j23, M=None)` — (7.3.49)：Σ C^{j1 j2 j12}C^{j12 j3 J}C^{j1 j23 J}C^{j2 j3 j23}（固定 M）
- `tensor_T(j1, j2, k, q)` — (7.3.1a)／(7.3.2)：T(j1 j2)^k_q = Σ C^{j1 j2 k}_{m1 m2 q}|j1 m1⟩⟨j2 −m2|(−1)^{j2−m2}（列 m1、行 m'）。
- `unit_v(j1, j2, k, q)` — (7.3.17a)：v(j1 j2)^k_q = (−1)^{2j1}Σ C^{j1 j2 k}_{m1 −m2 q}(−1)^{j2+m2}|j1 m1⟩⟨j2 m2|。
- `spin_half_T()` — (7.3.3)：自旋 ½ 的 T^k_q（k = 0, 1）。
- `spin_half_cart()` — (7.3.4)：T_x = −(T¹₋₁ − T¹₁)/√2、T_y = −i(T¹₋₁ + T¹₁)/√2、T_z = −T¹₀ → σ/√2。
- `rotated_Jz_check(beta: float=37.0)` — (7.3.7)–(7.3.8)：R(0β0)J_zR† = sin β J_x + cos β J_z，且 R T¹₀ R† = Σ𝒟¹_{q0}T¹_q（自旋 ½）。→ 最大誤差。
- `expectation_rotation(beta0: float=60.0, beta: float=45.0)` — Fig. 7.3.1、(7.3.9)–(7.3.11)：⟨J_z(轉過)⟩ = sin β⟨J_x⟩ + cos β⟨J_z⟩ ＝ ½(p↑ − p↓)（(7.3.10a) 的機率式）。
- `tensor_j1()` — 
- `j1_relations()` — (7.3.13)：T¹₋₁ = J₋/2、T¹₀ = J_z/√2、T¹₁ = −J₊/2（j = 1）。
- `v_table(j, k)` — Table 7.1–7.2 的濃縮形式：v^k_q（q = −k…k）疊在同一個 (2j+1)² 陣列；每條對角線（q 與 −q 同一個）寫成「分子 / √分母²」。
- `v_tables_check()` — Table 7.1–7.2 的 j ≤ 2 格子：分子與分母都和書上相同（True）；也核對每條對角線正交歸一。
- `elementary_expansion(j, a: int, b: int)` — E_ab（1 起算的列 a、行 b）＝ Σ c^k_q v^k_q，c^k_q = conj(v^k_q[a, b])（v 在跡內積下正交歸一）。→ {(k, q): c}（非零者）。
- `mixed_examples()` — (7.3.17)–(7.3.19) 的例子：√10 V¹₁(pd)（3×5）與 √10 V¹₋₁(dp)（5×3）；Table 7.4 (d)(p) 的 v¹、v²、v³ 分母。
- `mixed_table(j1, j2, k)` — Table 7.4 的長方形 v(j1 j2)^k：(分子, {q: 分母²})（q = m1 − m2；q 與 −q 同分母）。
- `fpoly_check()` — (7.3.30)–(7.3.32) 的多項式（Condon–Shortley 的 X³_q）：xyz = −i(X³₂ − X³₋₂)/√30（書印 −i(X³₋₂ − X³₂)/√30，差整體號）、
- `py_check(theta: float=0.7, phi: float=1.1)` — (7.6.51)：p_x = (−Y¹₁ + Y¹₋₁)r√(4π/3)/√2 = x、p_y = i(Y¹₁ + Y¹₋₁)r√(4π/3)/√2 = y（書印 i(Y¹₁ − Y¹₋₁)，給 −ix）。
- `handedness_check()` — Fig. 7.5.4：|r⟩ = (|x⟩ + i|y⟩)/√2 的 E(Φ) = Re⟨·|r⟩e^{−iΦ} 隨 Φ 逆時針轉（右旋，正文）；圖說把 |l⟩、|r⟩ 說反了。
- `we_v_check(cases=((1, 1, 1), (2, 2, 4), (2, 1, 1), (R_(3, 2), HALF, 1), (3, 3, 6)))` — (7.3.21)–(7.3.22)：⟨j1 m1|v^k_q|j2 m2⟩ = C^{k j2 j1}_{q m2 m1}(−1)^{k+j1−j2}√((2k+1)/(2j1+1))，
- `we_operator_check(j=2, euler=(30.0, 50.0, 70.0))` — 由 J 做成的張量算子：T¹_q = J 的球分量、T²_q = Σ C^{1 1 2}_{q1 q2 q}J_{q1}J_{q2}。核對 (7.3.1b) 的轉動、
- `v4_poly_check()` — (7.3.25)：x⁴ + y⁴ + z⁴ − (3/5)r⁴ = 2(X⁴₄ + X⁴₋₄)/√70 + (2/5)X⁴₀；X⁴_q 也與 multipole_X(4, q) 相同。
- `crystal_matrix(j, k: int)` — ⟨V^{(k)}⟩_j（單位 D⟨j‖X^k‖j⟩；k = 6 時是 E⟨j‖X⁶‖j⟩）＝ Σ c_q v^k_q × (−1)^k√((2j+1)/(2k+1))（(7.3.23)–(7.3.24)）。
- `o_label(P, J)` — 把本徵子空間（投影 P）依 O 的特徵標分解成標記（"T_2"、"A_1+E"…；整數 J）。
- `crystal_levels(j, k: int)` — crystal_matrix 的本徵值（精確，sympy）＋ O 標記（數值投影）→ [(值, 標記, 簡併)]（由低到高）。
- `crystal_direct(j: int, k: int, n_th: int=24, n_ph: int=48)` — 第二種算法：⟨j m'|V^{(k)}(θφ)|j m⟩ = ∫Y*_{jm'} V Y_{jm} dΩ（r = 1，所以單位是 D⟨r^k⟩；Gauss–Legendre × 等分）→ 本徵值（由低到高）。
- `reduced_X(lp, k, l)` — (7.3.40b)：⟨l'‖X^k‖l⟩/⟨r^k⟩ = √((2l+1)/(2l'+1)) C^{k l l'}_{000} = √(2l+1)(−1)^{k−l}(k l l'; 0 0 0)。→ (兩式)。
- `r4_hydrogen(n: int, l: int)` — 氫原子 ⟨r⁴⟩（單位 a0⁴，Bockasten 的式子）。
- `r_power_numeric(n: int, l: int, p: int)` — 第二種算法：用 sympy.physics.hydrogen.R_nl 直接積分 ∫r^{2+p}|R_nl|²dr（精確）。
- `d2_factor(L: int, k: int=4, l: int=2)` — (7.3.50)：兩個電子（l² 組態，總 L）的 ⟨[ll]L‖X^k‖[ll]L⟩/⟨l‖X^k‖l⟩ ＝ 兩項之和
- `d2_factor_numeric(L: int, k: int=4, l: int=2)` — 第二種算法：單電子 X^k_0 的矩陣（約化矩陣元 1：⟨l m'|X|l m⟩ = C^{k l l}_{0 m m'}），X⊗1 + 1⊗X 在耦合態 |[ll]L M=L⟩ 的期望值
- `re_rigid(n, J: int=10, A: float=0.2, B: float=0.4, C: float=0.6)` — (7.4.5)：E = J(J+1)(A n_x² + B n_y² + C n_z²)（n 是 J 方向的單位向量；可為 (…,3) 陣列）。
- `re_rigid_euler(beta, gamma, J: int=10, A: float=0.2, B: float=0.4, C: float=0.6)` — (7.4.4a)＋(7.4.5)：J_x = −|J|sin β cos γ、J_y = |J|sin β sin γ、J_z = |J|cos β。
- `re_rigid_tensor(beta, gamma, J: int=10, A: float=0.2, B: float=0.4, C: float=0.6)` — (7.4.10)：張量形式 J(J+1)[(A+B+C)/3 + (2C−A−B)/6·(3cos²β − 1) + (A−B)/2·sin²β cos 2γ]。
- `rigid_levels(J: int=10, A: float=0.2, B: float=0.4, C: float=0.6)` — H = A J_x² + B J_y² + C J_z²（(7.4.3)，2J+1 維）→ [(E, D2 標記)]（由低到高；標記由 180° 轉動 C2(x)、C2(y)、C2(z) 的期望值，Fig. 7.4.3）。
- `fig742()` — Fig. 7.4.2（J = 10、A, B, C = 0.2, 0.4, 0.6 cm⁻¹）：分隔線 BJ(J+1)、C2(x) 叢（低）與 C2(z) 叢（高）的成對能階與超精細分裂。
- `symmetric_top_check(J: int=10, B: float=0.2, C: float=0.6)` — (7.4.15)：E = BJ(J+1) + (C − B)K²（A = B）＝ 數值對角化。→ 最大誤差。
- `reduced_T2(J: int)` — (7.4.11)–(7.4.13)：⟨J‖T²‖J⟩ = √((2J+3)(2J+2)2J(2J−1))/4 ＝ ⟨J J|T²₀|J J⟩/C^{2 J J}_{0 J J}，T²₀ = (3J_z² − J²)/2。
- `d2_induced()` — (7.4.16)：C2(axis) 的 0₂、1₂ 誘導到 D2（阿貝爾群：χ↑(g) = [G:H]χ(g)（g ∈ H）、否則 0），用 Fig. 7.4.3 的標記分解。
- `rigid_correlation(J: int, Bs, A: float=0.2, C: float=0.6)` — Fig. 7.4.4：每個 B 的 2J+1 個能階（由低到高）→ (len(Bs), 2J+1)。
- `cone_angles(J: int, Ks)` — (7.4.2b)：Θ_{JK} = arccos(K/√(J(J+1)))（度）。
- `octa_re(n)` — (7.4.18)：n_x⁴ + n_y⁴ + n_z⁴ − 3/5（RE 曲面的形狀；t > 0 時 C4 軸是山、C3 軸是谷、C2 是鞍）。
- `octa_re_euler(beta, gamma)` — (7.4.19) 的角函數：(35cos⁴β − 30cos²β + 3 + 5sin⁴β cos 4γ)/20 ＝ octa_re。
- `octa_levels(J: int=30, t: float=1.0, B: float=0.0)` — 
- `octa_clusters(J: int=30, t: float=1.0, n_c4: int=6, n_c3: int=2)` — J 大時的能階叢：由高往下取 n_c4 個 6 態的 C4 叢（K4 = J, J−1, …；標記＝ (K4 mod 4)↑O），由低往上取 n_c3 個 8 態的 C3 叢
- `octa_model_04()` — 0₄↑O 的穿隧模型：六個 C4 方向（SO.SITES）之間的近鄰（夾 90°）穿隧 S = 1 → 本徵值 4（A1）、0（T1）、−2（E）。
- `bottom_sequence(J: int=30, n: int=10)` — (7.4.22)：由最低往上的 n 個標記（J = 30：A1 T1 T2 A2 T2 E T1 T2 E T1 ＝ 0₁↑O 的 24 維）。
- `separatrix_angles()` — 分隔線（過 C2 鞍點的等能量線）與 C4、C3 軸的夾角：沿 γ = 45° 的子午線解 E(β) = E(C2)。數值 ＝ arccos√(2/3) 與 54.7356° − 35.2644°。
- `cluster_cutoffs(J: int=30)` — K4 = √(J(J+1)) cos 35.26°、K3 = √(J(J+1)) cos 19.47°（J = 30：24.90、28.75）。
- `t46_operator(J: int, nu: float, sign4: int=1)` — (7.4.23)–(7.4.25)：T^{4,6}(ν) = √(9/4π)T^{[4]}cos ν + √(13/4π)T^{[6]}sin ν（球諧歸一因子照 (7.4.26)–(7.4.27)），
- `t46_sign_check(J: int=8)` — 與 C3（繞 (1,1,1) 120°）的對易子：T^{[4]} 用 + 號為 0、用 − 號不為 0；T^{[6]} 用書上的 − 號為 0。
- `t46_levels(J: int, nu: float)` — 
- `t46_clusters(J: int=30, nu_units: float=2.0)` — ν = nu_units·π/6 時由最低往上的兩個 12 態叢的標記（C2 叢：0₂↑O、1₂↑O）。
- `re46(n, nu: float)` — (7.4.26)–(7.4.27)：E^{[4]}cos ν + E^{[6]}sin ν（n 為單位向量，可為陣列）。
- `sphere_contours(f, levels, nb: int=91, ng: int=180)` — f(n)（n 是 (…,3) 單位向量）在球面上的等高線 → {level: [n×3 折線, …]}（β ∈ [0, π] 含兩極、γ 週期；閉曲線首尾相同）。
- `resample(path, n: int=60)` — 折線依弧長重取 n 點（場景每條曲線 ≤ 60 點）。
- `rigid_paths(J: int=10, A: float=0.2, B: float=0.4, C: float=0.6, n: int=60)` — J = 10 不對稱陀螺的「量子化路徑」：RE 曲面 (7.4.5) 在每個本徵能量的等高線 → [(E, 標記, [折線…])]。
- `octa_paths(levels, n: int=60)` — 
- `re46_paths(nu: float, levels, n: int=60)` — 
- `spinor_euler(alpha: float, beta: float, gamma: float)` — (7.5.2)：(e^{−iα/2}cos(β/2), e^{iα/2}sin(β/2))e^{−iγ/2}（度）。
- `spin_vector(psi)` — Fig. 7.5.2：(⟨J_x⟩, ⟨J_y⟩, ⟨J_z⟩) = (Re ψ1*ψ2, Im ψ1*ψ2, (|ψ1|² − |ψ2|²)/2)。
- `h_to_omega(H)` — (7.5.5e,f)：ω0 = (H11 + H22)/2（書印 H12）、ω_x = 2Re H21、ω_y = 2Im H21、ω_z = H11 − H22；ω、θ、φ（度）。
- `evolve(H, psi0, t: float)` — 
- `crank_evolve(H, psi0, t: float)` — (7.5.7)：|ψ(t)⟩ = e^{−iω0 t}R[φ θ ωt]|ψ(0)⟩（R 是 RT.su2_darboux）。
- `osc_hamiltonian(A: float, B: float, C: float, D: float)` — (7.5.8)：[[A, B − iC], [B + iC, D]]。
- `jvars(x1, p1, x2, p2)` — (7.5.13)（更正 J_z 的減號）：J_x = x1x2 + p1p2、J_y = x1p2 − x2p1、J_z = (x1² − x2² + p1² − p2²)/2、J_0 = I/2。
- `hc_value(A, B, C, D, x1, p1, x2, p2)` — (7.5.10b) 的 H_c。
- `hc_check(n: int=5, seed: int=7)` — (7.5.14a) ＝ (7.5.10b) ＝ ⟨ψ|H|ψ⟩/2（隨機的 A…D 與相空間點；古典 H_c 是量子期望值的一半）。→ 最大誤差。
- `lissajous(A, B, C, D, psi0, ts)` — 解析解：ψ(t) = e^{−iHt}ψ0 → 陣列 (len(ts), 4)：x1, p1, x2, p2（x = Re ψ、p = Im ψ）。
- `lissajous_numeric(A, B, C, D, psi0, ts)` — (7.5.10a) 的 Hamilton 方程用 solve_ivp 積分（第二種算法）。
- `pauli_bases()` — (7.5.21)–(7.5.22)：C(σ_A,B,C) = σx、σy、σz（圓基底），L(σ) = U C(σ) U†（線基底）→ σz、σx、σy。
- `stokes_of(psi_xy)` — 線基底的態 → Stokes 向量 (S_A, S_B, S_C) = (⟨σz⟩, ⟨σx⟩, ⟨σy⟩)（(7.5.24)）。
- `state_circular(R: float, L: float, phi: float, Phi: float=0.0)` — (7.5.19)：(R e^{−iφ}|r⟩ + L e^{iφ}|l⟩)e^{−iΦ} → 線基底的兩分量（度）。
- `state_linear(X: float, Y: float, vt: float, theta: float=0.0)` — (7.5.17)：(X e^{−iϑ}|x⟩ + Y e^{iϑ}|y⟩)e^{−iθ}（度）。
- `stokes_circular(R: float, L: float, phi: float)` — (7.5.25) 左欄：I(cos 2ψ cos 2φ, cos 2ψ sin 2φ, sin 2ψ)，a = (R+L)/√2 = √I cos ψ、b = (R−L)/√2 = √I sin ψ。
- `stokes_linear(X: float, Y: float, vt: float)` — (7.5.25) 右欄：I(cos 2ν, sin 2ν cos 2ϑ, sin 2ν sin 2ϑ)，X = √I cos ν、Y = √I sin ν。
- `ellipse_field(psi_xy, Phis)` — (7.5.15)：E_x = Re⟨x|Ψ⟩e^{−iΦ}、E_y = Re⟨y|Ψ⟩e^{−iΦ}（Φ 是度；畫偏振橢圓用）→ (len, 2)。
- `ellipse_circular_formula(R: float, L: float, phi: float, Phis)` — (7.5.20a,b)：E_x = a cos Φ cos φ − b sin Φ sin φ、E_y = a cos Φ sin φ + b sin Φ cos φ。
- `ellipse_card()` — 本系列的橢圓：φ = 20°、ψ = 15°、I = 1 → R、L、a、b、Stokes、線基底的 ν、ϑ、X、Y、Euler 角 (7.5.26)、(7.5.27)。
- `faraday_check(C: float=0.3, A: float=1.0, t: float=2.0)` — (7.5.29)：A = D、B = 0、C ≠ 0 ⇒ ω 沿 C 軸，S 繞 C 轉 2Ct：橢圓的 ψ 不變、長軸 φ 轉 Ct。
- `birefringence_check(A: float=1.5, D: float=1.0, t: float=2.0)` — (7.5.30)：B = C = 0、A ≠ D ⇒ ω 沿 A 軸，S 繞 A 轉 (A − D)t：盒角 ν 不變、ϑ 變 (A − D)t/2。
- `beat(A, B, C, D)` — (7.5.31)–(7.5.32)：|ω| = √((A−D)² + 4B² + 4C²)（半拍頻 ＝ 兩本徵頻率差）、ω0 = (A+D)/2；本徵值 ω0 ± |ω|/2。
- `closure(NH: int, NL: int)` — (7.5.33)–(7.5.34)：ω_H/ω_L = N_H/N_L ⇒ ω0/Ω = (N_H + N_L)/(N_H − N_L)（Ω = |ω|/2）、復現時 ω0Υ = π(N_H + N_L)（玫瑰瓣數）、ωΥ = 2π。
- `closure_numeric(NH: int=2, NL: int=1, phi: float=75.0, theta: float=90.0)` — 數值：本徵頻率 ω_H = N_H、ω_L = N_L、曲柄方向 (φ, θ)（Fig. 7.5.10 的 75°、90°）的 Hamilton 量，積分到 Υ = 2π/|ω_L| 回到原點、Υ/2 不回。
- `anharmonic_fixed_points(Omega: float=1.0, ax: float=2.0, I: float=1.0)` — (7.5.39a)–(7.5.40)：H = (I/2)Ω cos β + a_x(I²/4)cos²α sin²β（減去 J0ω0）。α = 0 上 dH/dβ = 0 ⇒ β = 0、π、arccos(Ω/(a_x I))。
- `osc_re(n, omega=(2.0, 0.0, 0.0), ax: float=0.0, I: float=1.0)` — 準自旋 RE 曲面（Fig. 7.5.11–7.5.13）：H = (I/4)(ω_A n_A + ω_B n_B + ω_C n_C) + a_x(I²/4)n_B²（n = (n_A, n_B, n_C)，畫圖時 A→x、B→y、C→z）。
- `local_mode_levels(j=5, Omega: float=1.0, ax: float=0.5)` — 量子版：H = Ω J_A + a_x J_B²（自旋 j = I/2；(7.5.35)–(7.5.38)）的本徵值（由高到低）；a_x I ≫ Ω 時最高的能階成對（局域模雙重態）。
- `h2p_integrals(R: float, Z: float=1.0)` — (7.6.21b)–(7.6.23)（原子單位，更正版）：S = (1 + ZR + Z²R²/3)e^{−ZR}、A = (a|1/r_b|a) = 1/R − (1/R + Z)e^{−2ZR}、
- `h2p_quadrature(R: float, book: bool=False)` — 球體座標 (μ, ν) 的數值積分（(7.6.21a)：d³x = (R³/8)(μ² − ν²)dμ dν dφ）。book=True 用 (7.6.22) 照印的被積函數
- `h2p_V(R: float)` — (7.6.24)：V^± − ε_1s = (1/R)[(1+R)e^{−2R} ± (1 − 2R²/3)e^{−R}]/[1 ± (R²/3 + R + 1)e^{−R}]。
- `h2p_V_generalized(R: float)` — 第二種算法：廣義本徵值問題 (7.6.10b) det|H − ES| = 0（h = ε − A、s = εS − B，ε = −½）＋ 1/R，再減 ε。
- `h2p_matrix_check(R: float=2.0)` — (7.6.12)–(7.6.13)：S⁻¹H 的元素 (h − Ss)/(1−S²)、(s − Sh)/(1−S²)（更正：左因子是 [[1, −S], [−S, 1]]）；照印的 [[1, S], [S, 1]] 不對。
- `h2p_min()` — V⁺(R) − ε_1s 的極小（LCAO）：R = 2.4928 a.u.、−0.0648 a.u.（書：≈ 2.5、−0.065；精確解 2.00、−0.1026）。
- `two_center(f1, f2, R: float, L: int=40, nmu: int=300)` — ∫∫ρ1(1)ρ2(2)/r12（兩個軸對稱電荷，原點在核間中點、核 a 在 +z）：ρ_i = Σ_l ρ_il(r)P_l(cos θ)，
- `h2_integrals(R: float, book: bool=False)` — (7.6.37) 的七個積分（原子單位，Z = 1）：A、B 同 H₂⁺；C = 5/8；D（(7.6.32a)，更正 R³/6；book=True 照印 R³/3）；
- `h2_integrals_numeric(R: float)` — 第二種算法：C、D、E、F 用 two_center（多極展開的數值積分）。
- `hl_energy(R: float, book: bool=False, triplet: bool=False)` — (7.6.30b)／(7.6.31b) 減 2ε_1s 再加 1/R：[−2A + D ∓ 2BS ± E]/(1 ± S²) + 1/R。
- `hl_min(book: bool=False)` — Heitler–London ¹Σg 的極小：R = 1.6425、−0.1160 a.u.（書 1.6、−0.115）；照印 R³/3：1.6220、−0.1275。
- `h2_generalized(R: float)` — 
- `vb_matrix(R: float, book: bool=False)` — (7.6.40)（非 Hermite；列＝行 ¹Σg ion、¹Σg cov、¹Σu ion、³Σu cov）；book=True 用照印的對角分母 (1 − S)²。
- `mo_matrix(R: float, book: bool=False)` — (7.6.47)（列＝行 ¹Σg(σg)²、¹Σg(σu)²、¹Σu σgσu、³Σu σgσu）；book=True 用照印的非對角分母 (1 − S)²（應 1 − S²）。
- `mo_min()` — (σg)² 單一組態（MO）的極小：M[0,0] + 1/R + 1。
- `ci_min()` — 兩個 ¹Σg 組態交互作用（＝ 1s 基底的完整 CI）的極小：最低廣義本徵值 + 1/R + 1。
- `limit_check(R: float=1.4)` — (7.6.48)–(7.6.49)：只留 C、G 時 MO 與 VB 的 ¹Σg 子矩陣本徵值都是 0 與 (1+S²)(C+G)/(1−S²)²（更正；書印 (1+S²)²）。
- `fig763_levels(s: float, R: float=1.4)` — Fig. 7.6.3 的數值版（本系列）：單電子積分 A、B 乘 s、兩電子積分 C D E F G 乘 (1 − s)；4×4 廣義本徵值減 2ε。
- `hybrid_angle(name: str)` — 等價混成的鍵角：s 成分 λ² = 1/(n+1)（sp^n），cos θ = −λ²/(1 − λ²)。sp 180°、sp² 120°、sp³ 109.4712°。
- `sp3_vectors_angle()` — 第二種算法：(7.6.53) 的 C1、C2 的 p 部分 (1,1,1)、(1,−1,−1) 的夾角，且 ⟨C_i|C_j⟩ = δ。
- `xy2_hybrids(two_phi: float, book_sign: bool=False)` — (7.6.54)–(7.6.56)：X1, X2 = λs + a cos φ p_x ± a sin φ p_y、X3 = λ's − a'p_x；a² = 1/(1 − cos 2φ)、λ² = cos 2φ/(cos 2φ − 1)、
- `xy2_config(two_phi: float)` — (7.6.60)–(7.6.61)：x1x2x3²p_z² = s^{2/(1−cos 2φ)} p^{(4−6cos 2φ)/(1−cos 2φ)}（指數由 (7.6.55)–(7.6.56) 加總 ＝ 公式）。
- `huckel_ring(n: int, H: float=0.0, Sx: float=1.0)` — 
- `huckel_chain(n: int, H: float=0.0, Sx: float=1.0)` — 
- `benzene()` — Hückel 苯：ε_k = H − 2S cos(2πk/6) → H − 2S、H − S(2)、H + S(2)、H + 2S；6 個 π 電子 6H − 8S（書印 H − 8S），局域 6H − 6S，共軛 2S。
- `butadiene()` — 丁二烯：C10 環中間四個能階 H − 2S cos(kπ/5)（k = 1…4）＝ 4 原子鏈；4 個電子 4H − 2√5 S，局域 4H − 4S，共軛 (2√5 − 4)S = 0.4721S（書印 0.047）。
- `octa_sigma()` — (7.6.62)–(7.6.63)：Frobenius（SO.frobenius）：0₄↑O、A'↑O_h（C4v）、E↑O_h。
- `ligand_states()` — (7.6.65)：六個 σ 配位基 |1⟩…|6⟩（+z, −z, +x, −x, +y, −y）用對稱化的 |A1⟩、|E,1⟩、|E,2⟩、|T1,x⟩、|T1,y⟩、|T1,z⟩ 展開；
- `card_values()` — plan_ch07.md §6 每張卡上的數字與導讀第 7 章的數字（字串；小數照卡上位數）。
- `summary()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 20 秒）。

## `combkit.sym_dyn`

sym_dyn.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 8 章

常數：`PI`、`LINE_CARD`、`LAGR_CARD`、`FALL`、`FIG816_LEVELS`、`NBASIS`、`COH_CARD`、`ROTOR`、`QE`、`ME`、`HBAR`、`EPS0`、`CLIGHT`、`A0`、`EH`、`TAU_AU`、`TL_CARD`、`FIG853`、`FIG854_DELTAS`、`BLOCH_CARD`、`ADIAB`、`JC_CARD`、`RECAP`、`CARD_EXPECTED`

- `legendre_numeric(f, Xs, xs)` — (8.1.2c)：Y(X) = max_x [xX − f(x)]（凸的 f；xs 是細格點）。
- `line_family(f, xjs)` — Fig. 8.1.2：每點 (x_j, y_j) 生成直線 Y = x_j X − y_j → [(斜率, 截距)]。
- `line_envelope(f, fp, xjs)` — (8.1.2a,b)：切點 X_j = f'(x_j)、Y_j = x_j X_j − f(x_j)（f(x) = x²/2 時 X = x）。
- `circle_envelope(f, fp, R: float, xs)` — 一般接觸轉換 S(x,y:X,Y) = (X−x)² + (Y−y)² = R²（每點生成半徑 R 的圓）的包絡＝平行曲線（取上側）。
- `lagrange_card()` — 
- `fall_v0()` — 
- `fall_state(theta_deg: float, t, v0=None)` — 從原點以仰角 θ（度）、速率 v0 射出：(x, y, ẋ, ẏ)（(8.1.26)）。
- `SP_closed(r0, r1, T: float)` — (8.1.28)：S_P = m(Δx)²/2T + m(Δy)²/2T + (fT/2)(y1 + y0) − f²T³/24m。
- `SP_initial(v0vec, y0: float, T: float)` — (8.1.29)：S_P = mẋ0²T/2 + mẏ0²T/2 + fTy0 + fT²ẏ0 + f²T³/3m。
- `SP_numeric(theta_deg: float, T: float)` — 
- `SH_closed(theta_deg: float, T: float)` — (8.1.24)：s_x = mẋ0(x1 − x0)、s_y = (m²/3f)(ẏ1³ − ẏ0³)。
- `SH_numeric(theta_deg: float, T: float)` — ∫ p·dr = ∫ m v² dt。
- `trajectory(theta_deg: float, n: int=60, tmax: float=8.0)` — 
- `caustic(x)` — 拋物線焦散（安全拋物線）y = ε/|F| − |F|x²/4ε = 32 − x²/128。
- `wavefront(Sval: float, thetas)` — S_H(0:r) = Sval 的點（Fig. 8.1.6(b)）：沿每條仰角 θ 的軌跡解 S_H(t) = Sval。→ (n, 2)。
- `orthogonality(Sval: float, thetas, h: float=0.0001)` — 波前切線（對仰角的中心差分）與動量的夾角餘弦的最大絕對值（(8.1.19)：p·dr = 0 ⇒ 0）。
- `fig816_curves(levels=FIG816_LEVELS, n: int=60, ymin: float=-8.0)` — 畫圖用：{S_H: [段…]}，每段是 (k, 2) 的波前折線（仰角 −90°…270° 取 n 個；y < ymin 的部分切掉，所以一個 level 可能有好幾段）。
- `fall_checks()` — (8.1.28) ＝ (8.1.29) ＝ ∫L dt；(8.1.24) ＝ ∫p·dr；∂S_P/∂T = −H、dS_P/dT = L、∂S_H/∂ε = T；書上 (8.1.27)、(8.1.30b)ₓ 的錯。
- `fall_card()` — 
- `semiclassical_check()` — ψ = ψ0 e^{iS/ħ} 代入 iħψ_t = −ħ²ψ_xx/2m + Vψ：(1/ψ)(…) = −S_t − S_x²/2m − V + iħS_xx/2m（(8.1.38)）；ħ → 0 得 HJ (8.1.17)。
- `wkb_ratio(n: int, x: float)` — ħ|dp/dx|/p²（(8.1.39)）：諧振子 p = √(2(n+½) − x²)（ħ = μ = ω = 1）。
- `huygens_card()` — Fig. 8.1.7 的數字版（自由空間 S_H = p|r|、p = 1）：r_10 = (10, 0)、r_20 = (20, 0)；r'_10 在 10 號圓上偏 θ 使 |r'_10 − r_20| = 20。
- `ladder(N: int=NBASIS)` — 截斷的 a：a[n−1, n] = √n（(8.2.4b)）。
- `xp_ops(N: int=NBASIS)` — 
- `disp_D(x0: float, N: int=NBASIS)` — (8.2.8)：D(x0) = e^{−ix0 p}。
- `boost_B(p0: float, N: int=NBASIS)` — (8.2.10)：B(p0) = e^{ip0 x}。
- `Q_op(x0: float, p0: float, N: int=NBASIS)` — (8.2.11)：Q = e^{i(p0 x − x0 p)}。
- `alpha_of(x0: float, p0: float)` — (8.2.15b)：α = (x0 + ip0)/√2。
- `bch_check(x0: float=1.3, p0: float=-0.7, N: int=NBASIS, k: int=25)` — (8.2.12)–(8.2.15)：Q ＝ e^{−ix0p0/2}BD ＝ e^{ix0p0/2}DB ＝ e^{αa† − α*a} ＝ e^{−|α|²/2}e^{αa†}e^{−α*a}（比左上 k×k，避開截斷邊緣）；
- `coherent_cn(alpha: complex, nmax: int)` — (8.3.1b)：c_n = e^{−|α|²/2} αⁿ/√n!。
- `coherent_numeric(x0: float, p0: float, N: int=NBASIS)` — 
- `coherent_stats(alpha: complex, N: int=NBASIS)` — 
- `ho_eigfun(n: int, x)` — (8.2.7)：ψ_n(x)（μ = ω = ħ = 1），用穩定的遞迴 ψ_{k+1} = √(2/(k+1)) x ψ_k − √(k/(k+1)) ψ_{k−1}。
- `packet_moments(t: float, x0: float, p0: float=0.0)` — (8.2.22)：x_t = x0 cos t + p0 sin t、p_t = p0 cos t − x0 sin t。
- `packet_closed(x, t: float, x0: float, p0: float=0.0, sign: int=1)` — (8.2.28)：χ = e^{ix_tp_t/2 − it/2} N0 e^{−(x−x_t)²/2 + i·sign·p_t(x−x_t)}（sign = −1 是書上 (8.3.3) 的號）。
- `packet_sum(x, t: float, x0: float, p0: float=0.0, nmax: int=60)` — (8.3.1a)：Σ c_n ψ_n(x) e^{−i(n+½)t}。
- `packet_split_operator(t: float, x0: float, p0: float=0.0, L: float=24.0, n: int=1024, steps: int=400)` — FFT 分裂算子數值演化（獨立方法）：從 χ(x, 0) 出發 → (x, χ(x, t))。
- `packet_uncertainty(t: float, x0: float, p0: float=0.0)` — 數值 ΔxΔp（在 x 格點上直接積分）＝ ½。
- `energies(alpha: complex)` — (8.2.25)–(8.2.26)：E_c = |α|²、E_q = E_c + ½（書印 E_0 + ħω/2）。
- `SH_ho(t: float, x0: float, p0: float=0.0)` — (8.2.31)：∫p_t dx_t = (x_tp_t − x0p0)/2 + E_c t。
- `SH_ho_numeric(t: float, x0: float, p0: float=0.0)` — 
- `coh_x0()` — 
- `snapshot_sum(x, x0: float, m: int, N: int, sign: int=1)` — (8.3.3)–(8.3.4)：(1/N) Σ_p χ(x, t_p) e^{i(m+½)t_p}，t_p = 2πp/N（sign = −1 用書上 (8.3.3) 的 −ip_t 號）。
- `snapshot_error(x0: float, m: int, N: int, sign: int=1)` — L2 相對誤差 ‖Σ − c_mψ_m‖/‖c_mψ_m‖（格點積分）與 (8.3.8) 的預測 √(Σ_k c²_{m+kN})/c_m、最大誤差／峰值。
- `ac_ho(t, alpha: complex)` — ⟨α|e^{−iHt}|α⟩ = e^{−it/2} exp(|α|²(e^{−it} − 1))。
- `ac_ho_sum(t, alpha: complex, nmax: int=60)` — 
- `ftac(ts, ac, Es)` — (8.3.15)：(1/2T)∫_{−T}^{T} e^{iEt} AC(t) dt，AC(−t) = AC(t)*；ts 從 0 到 T 等間隔（梯形）。
- `ho_ftac_card(T: float=20 * PI)` — 
- `ho_maslov(T: float=20 * PI)` — 半古典 AC：|MU(t)⟩ = 平移的 Gauss（(8.2.28) 的大括號）× e^{iS_p}、S_p = ∫(p ẋ − E_c)dt ⇒ AC_sc = AC·e^{it/2}，
- `rotor_H(J: int | None=None)` — 
- `rotor_packet0()` — (8.3.24)：R(0, β0, 0)|J J⟩（RT.wigner_D 的第一行，β0 = −30°）。
- `rotor_exact()` — K = 8 叢的兩個量子能量、初始波包 ⟨H⟩、J² 約定的古典能量（書 E_classical = 50.000）、精確 AC 的 FTAC 峰（叢的加權中心）。
- `rotor_localized(K: int=8)` — (8.3.32)：K = 8 叢兩個本徵態中 ⟨J_z⟩ 最大的組合（局域在 +z 的 C 軸）在 K' = 10, 8, …, 0 的分量。
- `rotor_sc(T: float=40.0, nt: int=4001)` — 自旋相干態的時間相依變分（TDVP）半古典傳播（(8.3.25)–(8.3.30) 的做法）：
- `rotor_ebk(K: int=8, Jeff: float | None=None)` — EBK（轉動型軌道沒有轉折點，α = 0）：(1/2π)∮ J_z dφ = K，|J| = J + ½；E(n) = |J|²(A n_x² + B n_y² + C n_z²)。
- `gauge_check()` — 一維偶極近似：ψ^α = e^{iqA(t)x/ħ}ψ^β；若 iħψ^β_t = −ħ²ψ^β_xx/2μ + Vψ^β − qE xψ^β（E = −dA/dt，(8.4.15)），
- `hamiltonian_check()` — (8.4.4)–(8.4.5)：L = μv²/2 − qΦ + q v·A ⇒ p = μv + qA、H = v·p − L = (p − qA)²/2μ + qΦ ＝ (8.4.5c) 的展開。
- `boost_shift(E0: float=1000000.0, omega: float=1e+16, q: float=1.6e-19, mu: float=9.1e-31)` — (8.4.22)：E 座標相對 A 座標的古典平移振幅 qE0/μω²（書上的數：1.758e−15 m，書寫約 2e−15；q 印成 1.6e19）。
- `ap_er_ratio(w_nb: float, w: float)` — (8.4.51) 對 (8.4.35a)：天真的 A·p 矩陣元與 E·r 矩陣元的大小比 = ω_NB/ω（只在共振時相等）。
- `S_amp(Delta, tau)` — (8.4.38a)：∫_0^τ e^{iΔt}dt = (2/Δ)e^{iτΔ/2} sin(τΔ/2)（Δ = 0 時 τ）。
- `S2(Delta, tau)` — 
- `spectral_card(tau: float=4.0)` — |S|² 的峰 τ²（(8.4.41)）、零點 ±2π/τ（(8.4.42)）、面積 2πτ（(8.4.40)）、第一旁帶（u = Δτ/2 解 tan u = u）。
- `first_order_check()` — 弱驅動兩能階（ω_NB = 1、ω = 0.9、E·r 的實餘弦）：數值 |c_N(t)|² ＝ (8.4.37) 兩項（c ≈ (q r E0/2ħ)·i·[e^{iφ}S(Δ⁻) + e^{−iφ}S(Δ⁺)]）；
- `golden_rule_check(V: float=0.001, dD: float=0.01, t: float=30.0)` — 準連續帶（Δ_k = k·dΔ）：P(t) = Σ|V|²|S(Δ_k, t)|² ≈ 2π|V|²t/dΔ（ρ = 1/dΔ）→ 回傳比值。
- `ho_strengths(N: int=30)` — (8.4.47b)：f_{n0} = 2ω_{n0}|x_{n0}|²（μ = ħ = ω = 1）＝ δ_{n1}（(8.4.48)）。
- `h_fn(n: int)` — 氫 1s → np 的振子強度 f_n = 2⁸n⁵(n−1)^{2n−4}/3(n+1)^{2n+4}。
- `h_strengths(nmax: int=20000)` — 
- `h_box_strengths(Rmax: float=80.0, n: int=8000)` — 徑向方程 −½u'' + [l(l+1)/2r² − 1/r]u = Eu 的有限差分（盒 Rmax）：1s → 所有 l = 1 態，f = (2/3)ΔE|∫u_1s r u_p dr|²。
- `classical_factor_check()` — §8.4B(b) 的古典振子：|e^{−ia} − e^{−ib}|² = 4 sin²((a−b)/2) ⇒ |x|² 的係數是 4|(q/μ)E0|²（書印 |4(q/μ)E0|²，大 4 倍）。
- `crank(t, r: float, Om: float, eps: float)` — (8.5.6)：ω = (r cos Ωt, r sin Ωt, ε)。
- `h_lab(t: float, r: float, Om: float, eps: float, rwa: bool=True)` — (8.5.5)（rwa = True）或 (8.5.4a)（實的餘弦曲柄）。
- `evolve_lab(ts, r: float, Om: float, eps: float, psi0=(1, 0), rwa: bool=True)` — 
- `spin_precession(ts, r: float, Om: float, eps: float, S0=(0, 0, 0.5))` — dS/dt = ω(t) × S（數值）。
- `rotating_check(r: float=0.2, Om: float=0.9, eps: float=1.0, t: float=7.3)` — (8.5.11)：轉動座標的 H 常數；(8.5.13a) 的曲柄 (r, 0, Δ)；(8.5.15a) ＝ 數值；SC.h_to_omega 給 (8.5.6)。
- `rabi(Delta, r: float, t)` — (8.5.17)：|c2|² = r²/(Δ² + r²) sin²(√(Δ² + r²) t/2)。
- `rabi_numeric(Delta: float, r: float, ts, eps: float=1.0)` — 
- `rabi_card()` — 
- `floquet_population(Om: float, eps: float, r: float, nst: int=400)` — 全 Hamilton 量 (8.5.4a) 的一週期傳播子（中點、每步閉式指數）→ Floquet 模在一週期內 |1⟩ 布居的時間平均（兩模的較大者）。
- `bloch_siegert(eps: float=1.0, r: float=0.2)` — Floquet 共振（兩個 Floquet 模的時間平均布居都是 ½：|1⟩、|2⟩ 完全混合）→ Ω_res − ε（數值）；
- `bs_sign_check(r: float=0.2, Om: float=0.95, eps: float=1.0, t: float=0.4)` — (8.5.20)：轉動座標中全 H 的 H21 = ½(re^{−2iΩt} + r) ⇒ ω_y = 2Im H21 = −r sin 2Ωt（書 (8.5.21) 印 +）。
- `bloch_rhs(t, S, Delta, r, g1, g2, Sz0)` — (8.5.35a)（轉動座標、ω = (r, 0, Δ)）：dS/dt = ω × S − [γ2 S_x, γ2 S_y, γ1(S_z − S_z0)]。
- `bloch_steady(Delta, r, g1, g2, Sz0)` — 解析穩態：S_z = S_z0/(1 + r²γ2/γ1(γ2² + Δ²))、S_y = −rγ2S_z/(γ2² + Δ²)、S_x = −ΔS_y/γ2。
- `bloch_numeric(Delta, r, g1, g2, Sz0, T: float=600.0)` — 
- `thermal_sz(beta_eps: float)` — ρ = e^{−βH}/Z、H = (ε/2)σ_z（|1⟩ 在上）：⟨S_z⟩ = ½Tr(ρσ_z) = −½ tanh(βε/2)（書 (8.5.36) 後的式子把下能階叫 1：+½ tanh）。
- `dressed(Delta: float, r: float)` — (8.5.39)–(8.5.44)：θ = atan2(r, Δ)、R_D 對角化 H^{R'} 成 diag(ω/2, −ω/2)；|1⟩^D = (cos θ/2, sin θ/2) 的本徵值是 +ω/2。
- `ac_stark(Om, eps: float, r: float)` — (8.5.45)：上 ε + δ/2、下 Ω − δ/2，δ = √(Δ² + r²) − Δ、Δ = ε − Ω（＝ [[ε, r/2], [r/2, Ω]] 的本徵值）。
- `ac_stark_card()` — 
- `circle_geometry(Delta: float, r: float)` — Fig. 8.5.4：δ = ω − Δ、δ' = ω + Δ（δδ' = r²）、δ² + r² = 2ωδ（(8.5.46b)）、tan(θ/2) = δ/r、|1⟩^D 振幅 r/√(2δω)、δ/√(2δω)。
- `adiabatic(r: float=0.5, D0: float=20.0, T: float=400.0)` — 轉動座標 H = Δ(t)σ_z/2 + rσ_x/2，Δ 從 +D0 線性掃到 −D0：數值反轉 ＝ 1 − e^{−πr²/2v}（Landau–Zener，v = 2D0/T）。
- `photon_elements(n: int, N: int=40)` — 
- `mode_energy_check()` — (8.6.3)：⟨ε0E²/2 + B²/2μ0⟩ = 2ε0ω²|a|²（E0 = 2|a|ω、B0 = 2|a|k、k² = ω²ε0μ0）；照印的磁場項 4|a|²/2μ0 少 k²（因次不對）。
- `angular_sum()` — ∫dΩ_k Σ_α |e_α·r̂|² = ∫dΩ (1 − cos²θ) = 8π/3（Σ_α|e_α·r|² = |r|² − |k̂·r|²）。
- `angular_sum_mc(n: int=200000, seed: int=3)` — 同上，用隨機 k̂ 與兩個正交偏振向量直接相加（獨立方法）。
- `einstein_A(omega: float, r2: float, q: float=QE)` — (8.6.32b)（SI）：A = 4ω³q²|r|²/(3ħc³·4πε0)。
- `einstein_A_gauss(omega: float, r2_cm2: float)` — Gauss 單位：A = 4ω³e²|r|²/3ħc³（e = 4.80320e−10 esu、ħ 與 c 用 cgs）。
- `h2p_A()` — 氫 2p → 1s：|⟨1s|z|2p0⟩|² = 2¹⁵/3¹⁰ a0²（閉式 ＝ 徑向數值積分 × 角因子 1/3）、ω = (3/8)E_h/ħ。
- `golden_A_numeric(t_au: float=4000.0)` — (8.6.28)：P(t) = Σ_α∫dΩ∫dω ρ(ω)|⟨H_I⟩|²|S|²/ħ²（原子單位，偶極沿 z）→ P/t ＝ A（t 大時）。
- `number_vs_coherent(n: int=5, N: int=60)` — §8.6B(c)–(d)：數態 ⟨n|a|n⟩ = 0（「無能」）、⟨n|x²|n⟩ = n + ½；相干態 ⟨α|a|α⟩ = α（(8.6.35)）。
- `jc_matrix(Nmax: int, Om: float, eps: float, g: float, book_phase: bool=False)` — H'_JCM（(8.7.10) 的實版本；book_phase = True 用 (8.7.5) 的 i(g/2)(a†σ− − aσ+)）；基底序號 2n + s（s = 0 ↓、1 ↑）。
- `jc_block(N: int, Om: float, eps: float, g: float)` — (8.7.10)：{|N−1⟩|↑⟩, |N⟩|↓⟩}：[[(N−1)Ω + ε, g√N/2], [g√N/2, NΩ]]。
- `jc_block_eig(N: int, Om: float, eps: float, g: float)` — 
- `jc_levels(Nmax: int, Om: float, eps: float, g: float)` — 完整截斷矩陣的數值本徵值 ＝ 各 2×2 塊公式（＋基態 0）；book_phase 的 i 不改本徵值。
- `jc_splitting(N: int, g: float, Om: float=1.0)` — 共振（Δ = 0）時第 N 塊的分裂 g√N（數值 2×2 對角化）。
- `mollow(Om: float, Delta: float, r: float)` — (8.7.11)、(8.7.13)：Ω_F = Ω + ω、Ω_3 = Ω − ω（ω = √(Δ² + r²)）；遠離共振 Ω_F → ε = Ω + Δ、Ω_3 → 2Ω − ε；照印 CSRS 的 −δ/2 給 Ω − Δ。
- `mollow_jc(N: int=100, g: float=0.04)` — JC 本身（r = g√N 隨 N 變）：螢光 ＝ 第 N 塊上能階 − 第 N−1 塊下能階、CSRS ＝ 第 N 塊下能階 − 第 N−1 塊上能階（用 (8.7.10) 的 NΩ1 + … 形式）。
- `jc_ladder(Delta: float, g: float=0.1, Nmax: int=3, Om: float=1.0)` — Fig. 8.7.1：[(N, 下穿衣能量, 上穿衣能量, 裸 |N⟩|↓⟩ 能量, 裸 |N−1⟩|↑⟩ 能量)]（N = 1…Nmax；N = 0 是基態 0）。
- `recap_check()` — 回顧表裡能便宜重算的數：第 1 章 [[4,1],[3,2]] 的本徵值、第 5 章 (50,60,70) 的 Darboux 角、第 7 章 21 cm 的 a、第 2 章拍的時間、第 8 章 c5 與 A。
- `card_values()` — plan_ch08 §6 每張卡上的數字與導讀第 8 章的數字（字串；小數照卡上位數）。
- `summary()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 20 秒）。

## `combkit.sym_mat`

sym_mat.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 1 章

常數：`I`、`EPS`、`LAM`、`HBAR`、`PASS_ANGLES`、`E11`、`E22`、`SERIES`、`SX`、`SY`、`SZ`、`BEAT_AB`、`MB`、`NB`、`DEG4`、`JORDAN`、`RIGHT_VEC`、`LEFT_VEC`、`S_EX`、`OSC`、`DELTA_SIGMAS`、`P121`、`P122`、`P124`、`P125_DATA`、`P126_M`、`P126_N`、`P127_H`、`P127_A`、`P127_U`、`SUMMARY_EXPECT`

- `to_np(M)` — sympy 矩陣 → numpy complex 陣列。
- `dag(M)` — 轉置共軛 M†（(1.1.7)）；sympy 或 numpy 都可以。
- `is_unitary(M, tol: float=1e-12)` — (1.1.8c)：M†M = 1 = MM†。
- `U(theta)` — (1.1.3)：U(θ) = [[cos θ/2, −i sin θ/2], [−i sin θ/2, cos θ/2]]，θ 以度計（精確）。U(a)U(b) = U(a+b)。
- `U_num(theta: float)` — 
- `T(beta, alpha)` — T(β←α)_{ij} = ⟨i_β|j_α⟩：從角度 α 的分析器基底換到角度 β 的（書上 T(b'←b) = U(θ)，θ = α_b − α_b'）。
- `T_real(theta)` — 實的特例（正交矩陣，(1.1.10)）：[[cos θ/2, −sin θ/2], [sin θ/2, cos θ/2]]（繞另一個軸轉；通過率相同）。
- `pass_fraction(theta)` — |⟨1'|1⟩|² = cos²(θ/2)（精確）。
- `analyzer_chain(stages, start=(0, 1))` — 一串分析器。start = (角度, 通道)：入射束全在該分析器的態 |通道⟩。stages = [(角度, mode), …]，最後一個是計數用的分析器。
- `axiom4(theta_in, theta_mid, theta_out, k: int=1)` — Axiom 4：⟨k''|1'⟩ = Σ_j ⟨k''|j⟩⟨j|1'⟩。回傳 (|Σ 振幅|², Σ|振幅|²) = (相干, 不相干)；兩者之差＝干涉項。
- `series_table()` — 0° 的 |1⟩ 進 90° 分析器、再進 180° 分析器：四種情形的 (p1, p2)。
- `problem_111(theta=60)` — 習題 1.1.1：先分成 right/left（90°）、再合併、最後 up/dn（0°）計數；stopper 擋 2''。回傳 {入射態: (無 stopper, 有 stopper)}。
- `problem_111_closed(theta)` — 
- `zeno(n: int)` — 習題 1.1.2：n 個分析器依次轉 (k/n)·90°，每個只留通道 1：最後是 right 的比例 cos^{2n}(π/4n)（n → ∞ 趨於 1）。
- `analyzer_chain_num(stages, start=(0.0, 1))` — analyzer_chain 的浮點版（numpy，角度可為任意浮點；給場景的動畫用）。
- `zeno_chain(n: int)` — 
- `basis(k: int, n: int=2)` — 
- `elementary(i: int, j: int, n: int=2)` — (1.1.21) e_ij = |i⟩⟨j|。
- `ket_primed(theta, j: int=1)` — (1.1.13)：|j⟩ 在 θ 轉過的（primed）基底下的分量 = U(θ) 的第 j 行；j = 1：(cos θ/2, −i sin θ/2)。
- `change_basis(H, Tm)` — (1.1.23c) H' = T H T†。
- `outer(psi, phi)` — (1.1.22) |Ψ⟩⟨Φ| 的矩陣 = Ψ ⊗ Φ†。
- `H2(A=3, B=1)` — 
- `populations(t: float, A: float=3, B: float=1)` — |Ψ(0)⟩ = |1⟩：P₁(t) = cos²(Bt/ħ)、P₂ = sin²(Bt/ħ)（閉式）。
- `psi_t(t: float, A: float=3, B: float=1, psi0=(1, 0))` — e^{−iHt/ħ} 用投影算符：e^{−i(A+B)t}P₊ + e^{−i(A−B)t}P₋，P± = (1 ± σx)/2（(1.2.2)：本徵態只轉相位）。
- `populations_numeric(ts, A: float=3, B: float=1)` — 數值積分 (1.1.17) iħ dΨ/dt = HΨ（solve_ivp，實部虛部分開）。回傳 (P₁ 陣列, P₂ 陣列)。
- `secular(M)` — 久期多項式（首一）S(ε) = det(ε1 − M)；MB：ε² − 6ε + 5。
- `eigen_mult(M)` — {本徵值: 重數}（由久期多項式的根）。
- `distinct_eigs(M)` — 
- `hc_factors(M)` — (1.2.9) 未歸一化的 𝔭_j = Π_{l≠j}(M − ε_l 1)（相異本徵值）；MB：𝔭₁ = [[−1, 1], [3, −3]]、𝔭₅ = [[3, 1], [3, 1]]。
- `projectors(M)` — (1.2.12)/(1.2.15) 𝒫_j = Π_{l≠j}(M − ε_l)/Π_{l≠j}(ε_j − ε_l)（對可對角化的矩陣；重根時就是用最小方程）。
- `minimal_poly(M)` — 最小方程（MEq）：每個本徵值的冪次＝使 nullity((M − ε)^k) 等於代數重數的最小 k。
- `poly_at(poly: sp.Poly, M)` — 把多項式的 ε 換成矩陣（常數項乘 1）——Hamilton–Cayley (1.2.6)。
- `mat_func(M, f)` — 譜分解 f(M) = Σ f(ε_j) 𝒫_j（f 吃 sympy 數）。
- `eig_book(M)` — (1.2.13)：從 𝒫_j 取第一個非零對角元 N_j 所在的列與行，bra = 列/√N_j、ket = 行/√N_j。回傳 {ε: (bra 1×n, ket n×1)}。
- `adj_vectors(M, e)` — (1.2.25)：adj(M − ε_j 1) 的行是 ket、列是 bra 本徵向量；一般 adj(ε_j 1 − M) = 𝔭_j。
- `joint_projectors(A, B)` — (1.2.24)：對易的 A、B：𝒫^{α,β} = 𝒫_α(A)𝒫_β(B)，只回傳非零的 {(α, β): 矩陣}。
- `lagrange_basis(xs, x)` — Lagrange 內插的係數 ℓ_j(x) = Π_{k≠j}(x − x_k)/(x_j − x_k)（(1.2.12) 的「數」版本）。
- `ellipse_points(M, n: int=60)` — 單位圓經 M 變成的橢圓（n 點，2×n 實陣列）。
- `singular_values(M)` — 
- `iterate_angles(M, v=(1, 0), k: int=4)` — M^m v 的方向角（度），m = 0..k：MB、v = (1,0)：0、36.87、43.45、44.69、44.94 → 45°（ε = 5 的方向）。
- `angle_of(v)` — 
- `flux_in_out(S, A_in)` — (1.3.2)：Σ|A^I|² 與 Σ|A^o|²（A^o = S A^I）。
- `rho_pure(psi)` — (1.3.6) ρ = |Ψ⟩⟨Ψ|。
- `purity(rho)` — tr ρ²：純態 1、混合態 < 1。
- `force_matrix(k1, k2, k12)` — (1.4.2) 力矩陣 𝓕（F = −𝓕χ）。
- `mass_matrix(m1, m2)` — (1.4.3) 質量矩陣 m。
- `accel_matrix(m1, m2, k1, k2, k12)` — (1.4.4) 加速度矩陣 a = m⁻¹𝓕（χ̈ = −aχ；書上內文寫 a = −m⁻¹𝓕 是筆誤）。
- `omega2(case: str='equal')` — a 的本徵值 ω²（由小到大，精確）。
- `omega2_generalized(case: str='equal')` — 廣義本徵問題 (1.4.8b)–(1.4.13)：𝓕e = ω² m e，scipy.linalg.eigh(F, m) → (ω², E)，E 的行 m-歸一化（EᵀmE = 1）。
- `osc_analytic(t, case: str='beat', x0=(1.0, 0.0), v0=(0.0, 0.0))` — 用簡正模解 χ(t)：q_j = ⟨e_j|m|x⟩（(1.4.19b)），每個模態以 ω_j 振盪。t 可為陣列；回傳 2×len(t)。
- `osc_numeric(ts, case: str='beat', x0=(1.0, 0.0), v0=(0.0, 0.0))` — 直接數值積分 χ̈ = −aχ（solve_ivp）。
- `osc_energy(x, v, case: str='beat')` — H = T + V = ½⟨ẋ|m|ẋ⟩ + ½⟨x|𝓕|x⟩（(1.4.16)、(1.4.19)）。
- `beat_times(case: str='beat')` — (能量全部到 m₂ 的時間 π/(ω₂ − ω₁), 拍的週期 2π/(ω₂ − ω₁))。
- `modes_symbolic()` — 書上的一般例子：m₁、m₂、k₁、k₂、k₁₂ 符號。回傳 (a 矩陣, [ω²₋, ω²₊])。
- `spring_force_matrix(springs, n: int)` — Fig. 1.4.2 的「目視」規則：⟨i|F|j⟩ = Σ_彈簧 k·c_i·c_j（c_i = 彈簧方向在座標軸 i 上的投影，有正負）。
- `d1_matrix(N: int, L: float=2 * math.pi)` — 週期邊界的中央差分 d/dx（N×N）：反對稱 ⇒ 反 Hermitian；本徵向量是平面波，本徵值 i·sin(kh)/h。
- `d2_matrix(N: int, L: float=2 * math.pi)` — 
- `fourier_matrix(N: int)` — ⟨x_j|k⟩ = e^{2πijk/N}/√N（(1.5.6) 的離散版）：么正，Fourier 就是換基底。
- `d1_eigen(N: int, k: int=1, L: float=2 * math.pi)` — 
- `delta_gauss(x, s: float)` — δ 的高斯序列 δ_σ(x) = e^{−x²/2σ²}/(√(2π)σ)。
- `delta_integral(s: float)` — ∫δ_σ(x) cos x dx：(數值積分, 閉式 e^{−σ²/2}) → σ → 0 時趨於 cos 0 = 1（(1.5.3b)）。
- `delta_deriv_integrals(s: float)` — (1.5.10) 的符號：f = sin、x = 0。回傳 (∫dy f(y) ∂_y δ_σ(x − y), ∫dy f(y) ∂_x δ_σ(x − y))
- `adjoint_fd(mu, nu, lam, N: int=64, L: float=2 * math.pi)` — (1.5.7)、(1.5.12) 的有限差分版：M = diag(μ)D₂ + diag(ν)D₁ + diag(λ)；
- `resolvent(M, lam)` — (1.5.13) Q_λ(M) = (λ1 − M)⁻¹（精確）。
- `resolvent_adj(M)` — (1.5.16) Q_λ = adj(λ1 − M)/S(λ)：回傳 (adj(λ1 − M), S(λ))（λ 是符號 LAM）。
- `contour(M, center: complex, radius: float, f=None, n: int=64)` — (1.5.15)/(1.5.18) 數值圍道積分 (1/2πi)∮ f(λ)(λ1 − M)⁻¹ dλ，圓心 center、半徑 radius、n 點梯形法（對解析函數指數收斂）。
- `contour_PN(M, center: float, radius: float=0.5, n: int=64)` — (1.5.18)、(1.5.20)：(P_a, N_a) = ((1/2πi)∮dλ/(λ − M), (1/2πi)∮(λ − λ_a)dλ/(λ − M))。
- `inverse2(M, book: bool=False)` — 2×2 反矩陣；book=True 照書上左下分母印成 M₁₁M₂₂ − M₂₂M₂₁ 的寫法。
- `cramer(M, y)` — (B.10)：x_j = det(M 的第 j 行換成 y)/det M。
- `sqrt_all(M)` — 相異本徵值時的全部平方根 Σ ±√ε_j 𝒫_j（2ⁿ 個）。
- `secular_coeffs(M)` — 習題 1.2.3：S(ε) = εⁿ + a₁εⁿ⁻¹ + … + a_n。回傳三種算法：(多項式係數, (−1)^k e_k(本徵值), (−1)^k Σ k 階主子式)。
- `problem_125()` — 由本徵解倒推 M = Σ λ |v⟩⟨v|/⟨v|v⟩（Hermitian，向量互相正交）。
- `shift_matrix(n: int)` — 習題 1.2.8：截斷到 n 維的「升階移位」V = Σ|k+1⟩⟨k|（a† 的方向部分）。無限維時 V†V = 1 但 VV† = 1 − |0⟩⟨0|。
- `summary()` — 卡片與講義上的數字（字串）。
- `check()` — 

## `combkit.sym_nonabel`

sym_nonabel.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 3 章

常數：`R_`、`I`、`SQ3`、`EPS3`、`C3V_NAMES`、`C3V_TEX`、`MIRROR_DEG`、`C3V`、`C3V_TABLE_BOOK`、`C3V_CLASS_ALG_BOOK`、`C3V_CHARS_BOOK`、`C3V_IRREP_NAMES`、`D3_AXIS_DEG`、`D3`、`D3_TO_C3V`、`INV3`、`D4_NAMES`、`D4_TEX`、`D4`、`D4_TABLE_BOOK`、`D4_CLASS_TEX`、`D4_CHARS_BOOK`、`D4_IRREP_NAMES`、`D4_2D`、`C4V`、`D4_TO_C4V`、`D6`、`D3_CHARS`、`C2_CHARS`、`D6_CHARS_BOOK`、`D6_IRREP_NAMES`、`SIGMA_H`、`D3H`、`D3D`、`D4H`、`D2`、`D2H`、`D3H_CHARS_BOOK`、`D3H_IRREP_NAMES`、`D3D_IRREP_NAMES`、`D2_CHARS_BOOK`、`D2_LABELS_BOOK`、`D2_LABELS`、`D2H_IRREP_NAMES`、`D4H_IRREP_NAMES`、`QUAT`、`D2D`、`D2D_TO_D4`、`BOB_START`、`C3V_CLASS_IDEMP_BOOK`、`T_CIRC`、`GROUP_VECTORS_BOOK`、`PEND3`、`SIX`、`SIX_ANGLES`、`SIX_MASS`、`ROTOR`、`NH3`、`NH3_NU`、`NH3_OMEGA_BOOK`、`C_LIGHT`、`M_H_BOOK`、`NH3_KJ`、`P341`

- `cosd(deg)` — 
- `sind(deg)` — 
- `rot2(deg)` — xy 平面逆時針轉 deg 度：[[cos, −sin], [sin, cos]]（(3.1.4) 的 ⟨R[2φ]⟩、(3.6.19)）。
- `refl2(deg)` — 對「與 x 軸夾 deg 度的直線」（3D：含 z 軸的 φz 面）反射：[[cos 2φ, sin 2φ], [sin 2φ, −cos 2φ]]（Fig. 3.1.4(b)）。
- `rot3(axis, deg)` — 繞 axis（不必歸一）轉 deg 度的 3×3 精確矩陣（Rodrigues）。
- `e2(deg)` — 平面單位向量 (cos, sin)（浮點）。
- `block3(M2)` — 2×2 → diag(M2, 1)（鉛直軸不動）。
- **class `FiniteGroup`** — 由矩陣（精確 sympy）定義的有限群。乘積 ab = 先 b 後 a（書上的約定，p. 152）。
  - `.from_generators(cls, gens: dict, tex: dict | None=None, label: str='', identity: str='1')` — 由生成元（精確矩陣）做封閉：元素名是最短的字（"r s" 寫成 "rs"，連續相同字母寫成冪 "r^2"；ab = 先 b）。
  - `.mul(*gs)` — g1·g2·…·gk（最右邊的最先作用）。
  - `.inv(g: str)` — 
  - `.conj(h: str, g: str)` — h g h⁻¹（(3.2.4)、Definition p. 159）。
  - `.order_of(g: str)` — 
  - `.table_rows(order=None)` — 乘法表（列 a、行 b、格子 ab）。
  - `.is_abelian()` — 
  - `.is_subgroup(H)` — 
  - `.subgroup_of(g: str)` — g 生成的循環子群 [1, g, g², …]（依冪次）。
  - `.cosets(H, side: str='left')` — 左陪集 gH（side="left"）或右陪集 Hg（"right"），依 names 順序找代表元。
  - `.double_cosets(H, K)` — 雙陪集 HgK（習題 3.1.1(e)(f)）。
  - `.is_normal(H)` — 
  - `.commutator_subgroup()` — 所有 a⁻¹b⁻¹ab 生成的子群（習題 3.5.4）。
  - `.classes()` — 共軛類（依 names 中第一次出現的順序；類內也依 names 順序）。
  - `.class_of(g: str)` — 
  - `.normalizer(g: str)` — N_g = 與 g 對易的元素（p. 160）；Lagrange：°c_g = °G / °N_g（(3.2.7)）。
  - `.class_sizes()` — 
  - `.struct()` — 結構常數 (3.2.11)：struct[i, j][k] = γ^k_{ij}，c_i c_j = Σ_k γ^k_{ij} c_k。
  - `.class_mul(x, y)` — 類代數裡兩個元素（係數向量）相乘。
  - `.class_matrix(i: int)` — 類代數的正則表示（習題 3.5.1）：(C_i)_{kj} = γ^k_{ij}。
  - `.class_alg_rows(symbols=None)` — 類代數表 (3.2.1b)、(3.6.1b)：rows[i][j] = c_i c_j（sympy 式，用 symbols 或 c_1..c_n）。
  - `.min_poly_class(i: int, x=sp.Symbol('c'))` — 類算子 c_i 的最小方程 (3.2.12)–(3.2.13)（首一 Poly）：由 c_i 的冪在類代數中第一次線性相依求得。
  - `.class_idempotents_of(i: int)` — c_i 的最小方程的每個根 λ → 冪等元 Π(c_i − λ')/(λ − λ')（(3.2.14)，類係數向量）。
  - `.class_idempotents()` — 全對易冪等元 ℙ^α（類係數向量）：逐一用各類算子的最小方程分裂、兩兩相乘（§3.2B；D₄ 要兩個方程，p. 191）。
  - `.char_table_class_algebra()` — (維度 ℓ 的列表, 特徵標表 Matrix[α, 類])：只用類代數（§3.5A）。
  - `.ga(d: dict)` — 
  - `.ga_mul(x: dict, y: dict)` — 
  - `.ga_add(*xs, coeffs=None)` — 
  - `.ga_scale(x: dict, c)` — 
  - `.ga_dag(x: dict)` — x† = Σ x_g* g⁻¹（么正表示）。
  - `.ga_eq(x: dict, y: dict)` — 
  - `.el(g: str)` — 
  - `.class_op(i: int)` — 類算子 c_i = 類中元素之和（(3.2.1a)）。
  - `.from_class_vec(v)` — 
  - `.central_idempotents()` — ℙ^α 寫成群代數元素（與 class_idempotents 同順序）。
  - `.regular(g: str)` — 正則表示 (3.4.3)：𝓡_{hf}(g) = δ_{g, h f⁻¹}（列 ⟨h|、行 |f⟩，依 names 順序）；g|f⟩ = |gf⟩。
  - `.regular_of(x: dict)` — 
  - `.vec_of(x: dict)` — |x⟩ = x|1⟩ 的分量（(3.4.6a)：⟨g|x⟩ = x 裡 g 的係數）。
  - `.cyclic_projectors(gen: str)` — 循環子群 {1, g, g², …} 的投影算子 P^{k_n} = (1/n) Σ_p ε_n^{−kp} g^p（g P^{k} = ε_n^k P^{k}，ε_n = e^{2πi/n}）。
  - `.split(P: dict, projs: list)` — ℙ^α 乘上子群的投影算子，留下不為零的（(3.3.2)、(3.3.18)、(3.6.4)、(3.6.6)）。
  - `.elementary(Ps: list, ref: str | None=None)` — 由互相正交的不可約冪等元 [P_1, …, P_ℓ] 做基本算子 P_{ij}（(3.3.11)）：P_{1j} = P_1 g P_j / D_{1j}(g)，
  - `.irrep_from_elementary(Pij: dict)` — P_i g P_j = D_{ij}(g) P_{ij}（(3.3.11)）⇒ {g: D(g)}。
  - `.irreps(sub_gen: str | None=None, ref: str | None=None)` — 全部不可約表示（依 class_idempotents 的順序）：一維的由 gℙ = D(g)ℙ，多維的用 sub_gen 生成的循環子群分裂。
  - `.is_rep(D: dict)` — D(g)D(h) = D(gh) 對所有 g、h（(3.3.13)）。
  - `.char_table_traces(irreps: list)` — χ^α_k = Tr D^α(類 k 的代表元)（(3.5.1)）。
  - `.frequencies(chars, table: sp.Matrix | None=None)` — (3.5.11)：f^α = (1/°G) Σ_k °c_k χ^α_k* χ_k（chars：表示在每個類的跡）。
  - `.correlation(sub: list, sub_chars: sp.Matrix, table: sp.Matrix | None=None)` — 相關（subduction）表 (3.6.15)、(3.6.27)：n_{αβ} = (1/°H) Σ_{h∈H} χ^α(h) χ^β(h)*（sub_chars 的行依 sub 的順序）。
  - `.iso_to(other: 'FiniteGroup', mapping: dict)` — mapping（名稱 → 名稱）是同構：保持乘法表。
- `ga_tex(G: FiniteGroup, x: dict, factor=None, named: dict | None=None)` — 群代數元素 → TeX，提出公因數：{1: 1/3, r: −1/6, r²: −1/6} → \tfrac{1}{6}(2\mathbf 1-r-r^2)（依 names 順序）。
- `reorder(T: sp.Matrix, ref_rows)` — 把特徵標表 T 的列換成 ref_rows 的順序：回傳 perm，T[perm[a], :] = ref_rows[a]。
- `outer_table(T1: sp.Matrix, T2: sp.Matrix)` — 外積的特徵標表 (3.6.17)：kron(T1, T2)（列 (α, β)、行 (類 1, 類 2)，T1 的指標在外）。
- `cn_chars(n: int)` — C_n：χ^{k_n}(g^p) = ε_n^{kp}（列 k、行 p；本系列的約定 𝒟^{k_n}(r) = e^{2πik/n}）。
- `direct_product(G: FiniteGroup, extra: str, M, tex_extra: str, label: str)` — G × {1, X}（X 與 G 全部對易，例如 σ_h、I）：名稱 g、X.g（(3.6.28)–(3.6.30)）。
- `config_after(g: str, pts: dict | None=None, G: FiniteGroup | None=None)` — 群元素 g（鏡面固定在紙上）作用在擺錘的位置上：{名字: g·p}（浮點 2D）。
- `fig312()` — Fig. 3.1.2 的兩列：σ₂σ₁（先 σ₁）、σ₁σ₂（先 σ₂）的中間與最後位置；最後一格＝ r²、r 作用的結果。
- `book_order(G: FiniteGroup, ref_rows)` — (ℓ 列表, 特徵標表) 依書上的列順序（ref_rows）。
- `reflection3(N)` — 法向量 N 的鏡面反射 1 − 2NNᵀ（浮點）。
- `quat(axis, deg)` — 繞 axis 轉 deg 度的單位四元數 (cos θ/2, sin θ/2 n̂)。
- `qmul(p, q)` — Hamilton 乘積 pq（先 q 後 p 的轉動）。
- `quat_to_rot(q)` — 單位四元數 → 3×3 轉動矩陣。
- `su2(q)` — 四元數 → SU(2)：U = q₀1 − i(q₁σ_x + q₂σ_y + q₃σ_z)（第 1 章 (1.1.3) 的半角）。
- `axis_angle(Rm)` — 3×3 轉動 → (單位軸, 角度 0–180°)。
- `rotate_vec(v, axis, deg)` — 
- `hamilton_arc(axis, deg, start)` — Hamilton 轉向量（Fig. 3.1.5）：從 N1 = start（⊥ 軸、單位）沿大圓轉 deg/2 到 N2；R[ω] = σ_{N2} σ_{N1}。
- `arc_points(N1, N2, n: int=30, radius: float=1.0)` — N1 → N2 的大圓弧上 n 點（短弧，球半徑 radius）。
- `arc_quat(N1, N2)` — 弧 N1 → N2 所代表的轉動的四元數 (N1·N2, N1×N2)（轉 2∠(N1, N2)，軸 N1×N2）。
- `hamilton_product(axis1, deg1, axis2, deg2)` — R[ω'']= R[ω']R[ω]（先 R[ω] = (axis1, deg1)，Fig. 3.1.6）：兩段半角弧在兩個大圓的交點 P 頭尾相接，
- `conjugate_rotation(T_axis, T_deg, R_axis, R_deg)` — R' = T R T⁻¹（(3.2.4)、Fig. 3.2.1）：同一個轉角，軸換成 T·n̂。回傳 (四元數算的軸, 角, T·n̂)。
- `body_product(G: FiniteGroup, first: str, second: str)` — 本體座標（§3.1C、Fig. 3.1.8）：先做 first，鏡面跟著物體走；再做「跟著走過的」second，即 T·M₂·T⁻¹（T = first）。
- `d3_arcs()` — Fig. 3.1.7(a)：D₃ 的 Hamilton 弧（單位球）：r：赤道上 60° 弧；ρ_j：過極點、90° 弧（⊥ 轉軸的大圓）。
- `c3v_class_spectral()` — (3.2.17)：c_j = Σ_α c_j^{(α)} ℙ^α 的係數 c_j^{(α)}（列 α、行 j）＝ °c_j χ^α_j / ℓ^α。
- `c3v_split_cv()` — (3.3.1)–(3.3.3)：用 C_v'' = {1, σ₃} 的 P^± = (1 ± σ₃)/2 分裂 ℙ³ ⇒ (P³₁, P³₂)。
- `c3v_elementary(kind: str='cv')` — C3v 的 E 型基本算子 {(i, j): 群代數}（0 起算）：kind="cv" 用 {1, σ₃}、參考元素 σ₂（(3.3.7)、(3.4.8)）；
- `c3v_irreps(kind: str='cv')` — [D^{A1}, D^{A2}, D^E]（{g: Matrix}）：E 由基本算子求得（kind 見 c3v_elementary）。
- `c3v_regular_T()` — (3.4.7)–(3.4.8)、(3.4.16)：行 = |P̂¹⟩、|P̂²⟩、|P̂³₁₁⟩、|P̂³₂₁⟩、|P̂³₁₂⟩、|P̂³₂₂⟩（N^α = ℓ^α/°G，(3.4.20)）。
- `pendulum_rep(g: str)` — 三擺（擺錘在 0°、120°、240°，徑向擺動）的 3×3 置換表示：σ₃|1⟩ = |1⟩、σ₂|1⟩ = |r⟩、σ₁|1⟩ = |r²⟩（p. 189）。
- `pendulum_frequencies()` — (3.5.11)ₓ：三擺座標的 f^{A1}, f^{A2}, f^E（跡 3、0、1）。
- `pendulum_trace_eigs(a=PEND3['a'], b=PEND3['b'])` — (3.5.13)–(3.5.14)：a^α = (1/°G) Σ_類 °c_g χ^α_g* Tr⟨a g⟩（只用第一列）：A1 → b、E → 3a + b（A2 頻率 0，回傳 None）。
- `six_dof_first_row(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — (3.4.21b)：⟨1|a|g⟩，g = 1, r, r², σ₁, σ₂, σ₃（精確）。
- `six_dof_matrix(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — 6×6 加速度矩陣：⟨g|a|h⟩ = ⟨1|a|g⁻¹h⟩（a 與所有 g 對易，只要第一列）。
- `six_dof_geometric(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — 同一個矩陣，直接由彈簧幾何算（浮點）：k₁ 沿三邊、k₂ 從每個質點往外牆，與徑向夾 45° − α；座標方向 SIX_ANGLES。
- `six_dof_blocks(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — (3.4.22)：a^(1) = Σ_g ⟨1|a|g⟩、a^(2) = Σ χ^{A2}(g)⟨1|a|g⟩、a^(3)_{jl} = Σ_g D³_{jl}(g)*⟨1|a|g⟩（2×2，兩份相同）。
- `six_dof_blocks_book(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — (3.4.23b) 照抄。
- `six_dof_omega2(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'])` — {"A1": ω², "A2": ω², "E": (ω²₊, ω²₋)}（精確）。
- `six_dof_mode(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'], which: int=0)` — E 模（which = 0：高頻 ω₊、1：低頻 ω₋）的兩個夥伴 (|x⟩, |y⟩, ω)（浮點 6 維，正則座標）：
- `six_dof_displacements(x6)` — 6 個正則座標 → 三個質點的 2D 位移（列：0°、120°、240° 的質點）。
- `six_dof_circular(k1=SIX['k1'], k2=SIX['k2'], alpha=SIX['alpha'], m=SIX['m'], which: int=0, n: int=24)` — (3.4.24)–(3.4.25)、Fig. 3.4.5：x(t) = Re[e^{−iωt}(|x⟩ + i|y⟩)] = cos ωt|x⟩ + sin ωt|y⟩。
- `fig344_vectors(which: str='fixed')` — Fig. 3.4.4（k₂ = 0）的六個駐波向量（分量依 1, r, r², σ₁, σ₂, σ₃）：which="book" 照印、"fixed" 更正（型 2 兩欄）。
- `dn_label_free(n: int, m: int, kind: str)` — 自由轉子的駐波 cos mφ（kind="c"）／sin mφ（"s"）的 D_n 標記（由 m mod n）：
- `dn_free_levels(n: int, count: int)` — 自由轉子由低到高的前 count 個態 [(能量 m², 標記)]（E 二重算兩個）。
- `rotor_dn(n: int, V, mmax: int=40)` — D_n 位能 −V cos nφ 的轉子：(本徵值, 本徵向量(行), ms)（平面波基底，借 sym_abel.rotor_matrix）。
- `dn_label_vector(n: int, c, ms)` — 由數值本徵向量的對稱定標記：κ = m mod n（R：φ → φ + 2π/n）；κ ≡ 0、n/2 時看 φ → −φ（谷的 ρ）的宇稱。
- `dn_levels(n: int, V, count: int=8, mmax: int=40)` — 最低 count 個能階 [(能量, 由本徵向量定的標記)]。
- `dn_level_paths(n: int, Vs, count: int=8, mmax: int=40)` — 能階隨 V 變化（給「弱 → 強」的連線）：array[len(Vs), count]（第 i 條＝由低到高第 i 個），標記 labels[len(Vs)][count]。
- `dn_gap_2x2(n: int, V)` — 弱位能：m = ±n/2 的 2×2 [[(n/2)², G], [G, (n/2)²]]，G = −V/2 ⇒ (B1, B2) = ((n/2)² + G, (n/2)² − G)，能隙 2|G| = V。
- `dn_density(n: int, V, label: str, phis)` — 標記為 label 的（最低的那個）態的 |ψ(φ)|²（給「B1 波腹在谷、B2 在山頂」）。
- `dn_zeeman(n: int, V, b, count: int=4)` — 加磁場：H = (m − b)² + 位能（C_n 對稱，ρ 不再是對稱）；每個 C_n 不可約表示 κ = m mod n 各自對角化。
- `nh3_geometry(D=NH3['D'])` — 鍵長 1：H_p = (cos D cos 120p°, cos D sin 120p°, 0)、N = (0, 0, sin D)。
- `nh3_force(k, j, D=NH3['D'])` — 12×12 力常數矩陣（H–H 彈簧 k、N–H 彈簧 j，無預應力：Σ k n nᵀ）＝ (3.7.13a)。
- `nh3_mass(m=NH3['m'], M=NH3['M'])` — 
- `nh3_rep(g: str, D=NH3['D'])` — C3v 在 12 個座標上的表示（正交 12×12）：g 作用在原子位置與向量上（σ_j 是含 z 軸、鏡線 MIRROR_DEG 的鉛直面）。
- `nh3_projector(alpha: int, i: int, j: int, D=NH3['D'])` — P^α_{ij} = (ℓ^α/°G) Σ_g D^α_{ij}(g)* 𝓡(g)（(3.4.19)；α = 0, 1, 2 ↔ A1, A2, E；i, j 0 起算）。
- `nh3_cvectors(m=NH3['m'], M=NH3['M'], D=NH3['D'])` — (3.7.8)–(3.7.11) 的 |c¹⟩、|c¹'⟩、|c³₁⟩（圖上的歸一化：(3.7.11) 再除 √2）、|c³'₁⟩（12 維浮點）。
- `nh3_constraints(v, m=NH3['m'], M=NH3['M'], D=NH3['D'])` — (3.7.12)：T = Σ m c、A = Σ m r × c（v：12 維位移）。
- `nh3_reduced(k, j, m=NH3['m'], M=NH3['M'], D=NH3['D'])` — {"A1": (F, m, a), "E": (F, m, a)}：在 {|c¹⟩, |c¹'⟩}、{|c³₁⟩, |c³'₁⟩} 基底的 2×2（數值，由 12×12 直接投影）。
- `nh3_reduced_formulas(k, j, m=NH3['m'], M=NH3['M'], D=NH3['D'], book: bool=False)` — (3.7.20)–(3.7.24)、(3.7.27) 的 2×2：book=True 照印，False 更正後。回傳同 nh3_reduced 的格式。
- `nh3_secular(k, j, m=NH3['m'], M=NH3['M'], D=NH3['D'])` — (3.7.25)、(3.7.29)：S₃ = (3k/2 + Hj)/m、P₃ = (3kj/2m²)(H − cos²D/2)、S₁ = (3k + Kj)/m、P₁ = 3kj(1/m² + 3/mM) sin²D；
- `nh3_numeric(k, j, m=NH3['m'], M=NH3['M'], D=NH3['D'])` — 12×12 的 ω²（由小到大；前 6 個是平移與轉動的 0）。
- `nh3_omega_obs(nu=None)` — ω = 2πcν（10¹⁴ rad/s 為單位前的 rad/s）。
- `nh3_fit(m=NH3['m'], M=NH3['M'], D_assumed=NH3['D'], m_kg=M_H_BOOK, book_k: bool=False)` — (3.7.30)–(3.7.31) 與 p. 216：R = (ω₁₊ω₁₋/ω₃₊ω₃₋)²、D = arcsin √[R(1+3m/M)/(4(1+3m/M) − R(1−3m/M))]；
- `nh3_scan(kj, D=NH3['D'], m=NH3['m'], M=NH3['M'])` — Fig. 3.7.4：ω（以 √(j/m) 為單位）對 k/j：(ω₁₊, ω₁₋, ω₃₊, ω₃₋)。
- `nh3_ratio_fit(D=NH3['D'])` — 在 k/j 上找 ω₃₋/ω₁₋ = 實測比（書：「k/j = 0.3 附近」）；用兩個低頻定尺度，預測 ω₁₊、ω₃₊（10¹⁴ rad/s）。
- `nh3_crossing(D=NH3['D'])` — ω₁₊ = ω₃₊ 的偶然簡併（書：k/j = 0.035）。
- `nh3_modes(kj=NH3_KJ, D=NH3['D'], m=NH3['m'], M=NH3['M'])` — 四種真振動（給動畫）：[(名稱, ω（√(j/m) 單位）, 12 維位移（最大分量 1）)]，名稱 "1+"、"1-"、"3+"、"3-"（E 取型 1 夥伴）。
- `problem_311()` — 3.1.1：C_v = {1, σ₃} 的右／左陪集、C₃ 的右／左陪集；{1,σ₁}、{1,σ₂} 的雙陪集；{1,σ₁}、C₃ 的雙陪集。
- `problem_312()` — 3.1.2：{1, σ₃} 不是正規子群；C₃ 是（左陪集 ＝ 右陪集）。
- `problem_322()` — 3.2.2：°N(g)·°C_g = °G（D₂、D₃、D₄ 每個元素）。回傳 {群: [(g, °N, °C)]}。
- `problem_331()` — 3.3.1：D₄ 的類代數 (3.6.1b)、冪等元 (3.6.2)、譜分解係數 c_j^{(α)}、用 D₂ 與 C₄ 分裂 ℙ^E。
- `problem_332()` — 3.3.2、3.4.2：四元數群：類、類代數、特徵標表（＝ D₄ 的），E 用 {1, k, −1, −k} 分裂（複數）；Q ≇ D₄（2 階元素 1 個 vs 5 個）。
- `d4_square_matrix(e, s, d, m=1)` — 習題 3.4.1：正方形四角的質點（右上 = e(45°)），邊彈簧 s、對角彈簧 d、外接彈簧 e（沿對角線向外）。
- `problem_341(e, s, d, m=1)` — 3.4.1：塊對角約化（(3.4.22)：a^α_{jl} = Σ_g D^α_{jl}(g)* ⟨1|a|g⟩，E 用駐波 (3.6.7)）與 8×8 數值對角化。
- `problem_351()` — 3.5.1：C3v 類代數的 3×3 正則表示 R(c_j)；本徵值 c_j^α ⇒ χ^α_j = ℓ^α c_j^α / °c_j。
- `problem_352()` — 3.5.2：Σ_g D^α_{ij}(g)* D^β_{kl}(g) = (°G/ℓ^α) δδδ；Σ_α χ^α_j* χ^α_k = (°G/°c_j) δ_{jk}；Σ_j °c_j χ^α_j* χ^β_j = °G δ。
- `problem_355()` — 3.5.5：忠實表示的最低維數（群：最小的忠實不可約表示；群代數：Σ ℓ^α）、互相對易的群元素最多幾個、群代數的秩 Σ ℓ^α。
- `problem_356()` — 3.5.6：°G 被子群的階、類的大小、不可約表示的維數整除（C3v、D₄、D₆、Q 逐一核對）。
- `summary()` — {名稱: 字串}：場景卡與導讀上的主要數字。
- `summary_numbers()` — check() 用：轉子能階的 4 位小數字串。
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 5 秒）。

## `combkit.sym_oct`

sym_oct.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 4 章

常數：`R_`、`IU`、`EPS3`、`O_NAMES`、`R_AXES`、`X_AXES`、`I_AXES`、`O_TURNS`、`O_TEX`、`O`、`INV3`、`OH`、`O_CLASS_TEX`、`O_CLASS_SIZES`、`O_PRODUCTS_BOOK`、`COSETS_439_BOOK`、`DIAGONALS`、`O_PERM_BOOK`、`O_CLASS_ALG_BOOK`、`CI_POWERS`、`CI5_BOOK`、`CI_ROOTS`、`O_IRREP_NAMES`、`O_IRREP_ALT`、`O_CHARS_BOOK`、`O_CHARS`、`O_DIMS`、`O_CLASS_IDEMP_BOOK`、`O_CLASS_EIG`、`OH_IRREP_NAMES`、`OH_CHARS`、`OH_CLASS_TEX`、`TD_NAMES`、`TD`、`O_TO_TD`、`T_NAMES`、`TGRP`、`TH`、`T_IRREP_NAMES`、`T_CHARS`、`TH_IRREP_NAMES`、`TH_CHARS`、`GOLDEN`、`S_KOSHER`、`E_BASIS_Q`、`TETRA`、`TETRA_LABELS`、`T2_BOOK`、`T1_BOOK_REPS`、`E_BOOK_REPS`、`E_COSET_CLASSES`、`C4Z`、`D4Z`、`D2N`、`D2P`、`D3T`、`C3T`、`C2T`、`C2Z`、`C1`、`C4V`、`C3V`、`C2V`、`C2V_B`、`D4H`、`D3D`、`D2H`、`SUB`、`D4H_NAMES`、`D3D_NAMES`、`CORR_BOOK`、`P410_BOOK`、`T43`、`TMOV4`、`TMOV2`、`TRIG_U_T1`、`TRIG_U_T2`、`TRIG`、`TRIG_LABELS`、`TRIG_T1_BOOK`、`TRIG_T1_BOOK_TYPOS`、`MOV4`、`MOV4_LABELS`、`MOV4_T1_BOOK`、`D436_BOOK`、`TMOV3`、`_P3_T2`、`MOV3`、`MOV3_LABELS`、`BASES`、`SITES`、`SITE_LEADERS`、`SIX_LEVELS`、`SIX_VECTORS_BOOK`、`SIX_INDUCED_I4_BOOK`、`SIX_INDUCED_1_4_I4_BOOK`、`SIX_E_BOOK_446X`、`VQ_BASIS`、`VQ_LABELS`、`VQ_REDUCED_BOOK`、`VQ_SLOPES`、`VQ_HIGHQ`、`VQ_460_BOOK`、`SCAN`、`VD_BASIS`、`VD_T1U_BOOK`、`VD_LABELS`、`VD_REDUCED`、`VD_469`、`VD_469_RESULT`、`VD_SECOND_ORDER`、`V111_BASIS`、`V111_A`、`V111_E`、`T485_BOOK`、`F_SITES`、`A_LEADERS`、`B_LEADERS`、`XF6_NAMES`、`XF6_TEX`、`XF6_COORDS_LIST`、`XF6`、`UF6_NU`、`SF6_NU`、`SF6_ISO_EXP`、`XF6_DEG`、`PASCAL_464`、`SF6_NU_FIG`、`TD_IRREP_NAMES`、`PROBLEM_421_LABELS`、`CUBE_VERTS`、`CAGES`、`CAGE_GROUPS`、`CUBOCT`、`Y_CHARS`、`Y_CLASS_SIZES`、`Y_IRREP_NAMES`、`D5_CHARS`、`AMU`、`ATOM_MASS`、`XY4_NU`

- `diag_perm(g: str)` — g 對四條體對角線的置換 {k: k'}（g·d_k = ±d_k'）。
- `cycle_str(perm: dict)` — {1: 3, 3: 2, 2: 1, 4: 4} → "[132]"（每個循環從最小的數字寫起，長度 1 的省略）。
- `o_turn(g: str)` — (單位軸 np, 角度, Hamilton 弧長 = 角度/2)。
- `o_quat(g: str)` — 
- `element_of(Rm)` — 3×3 轉動（或 O_h 的矩陣）→ 名稱（O 或 OH）。
- `nomogram_product(a: str, b: str)` — ab（先 b 後 a）的 Hamilton 弧相接（Fig. 4.1.4）：SN.hamilton_product 的結果＋ 'name'（乘積的名稱）＋ 'sign'（四元數的號，自旋 ½）。
- `class_alg_rows()` — 類代數表（列、行依 O_CLASS_TEX）：rows[i][j] = c_i c_j 的類係數向量。
- `ci_powers(n: int=5)` — c_i^0 … c_i^n 的類係數向量（(4.1.7) 的正確版）。
- `o_class_idempotents()` — {irrep: ℙ^α 的類係數向量}（由 c_i 的最小方程 (4.1.8)–(4.1.9) 算，依 CI_ROOTS 命名）。
- `o_char_table_class_algebra()` — 只用類代數（§3.5A）算的特徵標表，列依書 (A1, A2, E, T1, T2)。
- `oh_split(g: str)` — O_h 的元素名 → (轉動部分名稱, 是否帶反演)。
- `chi_A2(g: str)` — 
- `signed_perm(code)` — 
- `sub_chars(key: str)` — 子群的特徵標表（列依 SUB[key][2] 的標記、行依 SUB[key][1] 的元素）。
- `correlation(key: str)` — 相關表 𝒟^α↓H（列 O 或 O_h 的 irrep（書上順序）、行 H 的 irrep，SUB[key][2]）。
- `conjugates(H: list)` — H 的共軛子群 tHt⁻¹（相異的，依出現順序）。
- `ga_proj(G, el: list, chars: list, ell: int=1)` — 子群（或群）的投影算子 (ℓ/°H)Σ χ(h)* h（群代數元素）。
- `o_central(alpha: str)` — ℙ^α（O 的群代數元素，(4.1.10)）。
- `d4_idempotents()` — D4 ⊃ D2 的六個不可約冪等元 (4.2.8)：P^E_1 = P^E p^{B1}、P^E_2 = P^E p^{B2}、P^{A1}、P^{A2}、P^{B1}、P^{B2}。
- `d3_idempotents()` — D3 ⊃ C2 的四個不可約冪等元 (4.2.24)：P^E P^A、P^E P^B、P^{A2}P^B、P^{A1}P^A（C2 = {1, i4}）。
- `chain_split(alpha: str, chain: str='tetragonal')` — ℙ^α 乘上子群鏈的不可約冪等元，留下不為零的（(4.2.9)、(4.2.18a)、(4.2.25)）。
- `elementary_from_irrep(D: dict, ell: int)` — P^α_{jk} = (ℓ/°G) Σ 𝒟_{jk}(g)* g（(3.4.19)、(4.4.1)），{(j, k): 群代數元素}。
- `p418(which: str='A1', prefactor=R_(2, 3))` — (4.2.18b,c)：prefactor·(P^{a}(D4) − ½ r1 P^{a} − ½ r1² P^{a})（書上印 3/4，正確 2/3）。
- `irrep(alpha: str, g: str, basis: str='tetragonal')` — O 的不可約表示矩陣 𝒟^α(g)（basis：tetragonal（O⊃D4⊃D2）、trigonal（O⊃D3⊃C2）、c4（O⊃D4⊃C4）、c3（O⊃D3⊃C3））；
- `site_perm(g: str, G=None)` — g 在六個頂點上的置換矩陣（⟨i|g|j⟩ = 1 若 g·site_j = site_i；O 或 O_h 的元素）。
- `site_of(g: str)` — 
- `induced(G, sub: list, Dsub: dict, leaders: list | None=None)` — 誘導表示 (4.3.34)：ℐ_{ij}(g) = D^a(l_i⁻¹ g l_j)（若在 H 內）否則 0；Dsub = {h: d×d 矩陣}；leaders 預設取每個左陪集
- `induced_char(G, sub: list, chi: dict)` — 誘導特徵標（公式）：χ(g) = (1/°H) Σ_{x∈G, x⁻¹gx∈H} χ_a(x⁻¹gx)，每個類一個值（依 G.classes）。
- `sub_irreps(key: str)` — 子群的不可約表示矩陣 {標記: {h: Matrix}}（一維由特徵標；D4 的 E ＝ T1 的 xy 塊；D3 的 E ＝ 三方 T1 的 E 塊；
- `frobenius(key: str)` — Frobenius 互反 (4.3.29)：對 H 的每個 irrep a，用誘導矩陣的跡與誘導特徵標公式各算一次 f^α(a↑G)，
- `six_H(H=Hs, T=Ts, S=Ss)` — (4.3.3)：對角 H、對面（1-2、3-4、5-6）T、其餘 S。
- `six_accel(k=sp.Symbol('k'), t=sp.Symbol('t'))` — (4.3.2)：古典版 −⟨a⟩ = (h, t, s)，h = 2k + t、s = k/2。
- `six_induced(label: str='0_4')` — D^{label}(C4)↑O，陪集代表 SITE_LEADERS（(4.3.9)–(4.3.12)、(4.3.34ₓ)）。
- `six_projected(alpha: str, j: int)` — P^α_{j,0_4}|1⟩ 歸一化（(4.3.19)、(3.4.19)：(ℓ/24) Σ_g 𝒟_{j a}(g)* g|1⟩；a = 0_4；T1 用行進波 (4.2.35)、E 用四方）。
- `six_coset_vector(alpha: str, j: int)` — 陪集分解 (4.3.44)：(ℓ°H/°G)^{1/2} Σ_l 𝒟_{ja}(l)* |l⟩（只對六個陪集代表求和）。
- `six_coset_energy(alpha: str, H=Hs, T=Ts, S=Ss)` — (4.3.46)：⟨ja|H|ja⟩ = Σ_l ⟨1|H|l⟩ 𝒟_{aa}(l)*（只用 (4.3.3) 的第一列）。
- `six_numeric(H: float=0.0, T: float=0.0, S: float=-1.0)` — 
- `d6_induced_examples()` — (4.3.49)–(4.3.51) 與正則表示：{名: 由誘導特徵標公式算的頻率（D6 依 SN.D6_IRREP_NAMES；O 依書）}。
- `vq_matrix(Q=Qs, H=Hs, S=Ss, T=0)` — (4.3.52)：⟨H + V_Q⟩（|1⟩、|2⟩ 加 Q；書上取 T = 0）。
- `vq_reduced(Q=Qs, H=Hs, S=Ss)` — 
- `vq_levels(Q: float, H: float=SCAN['H'], S: float=SCAN['S'])` — [(能量, 標記)]（由小到大）：E_u 二重、A_{2u}、兩個 A_{1g}（2×2）、B_{1g}。
- `vq_scan(Qs_=None, H: float=SCAN['H'], S: float=SCAN['S'])` — Fig. 4.3.5：{標記: array}（"A1g-"、"A1g+" 是 2×2 的低、高根；"Eu"、"A2u"、"B1g"）。
- `vq_avoided(S=Ss)` — 兩個 A_{1g} 的迴避交叉：最小間距在 Q = 2S，間距 4√2|S|（由 (4.3.60)）。
- `vd_matrix(D=Ds, H=Hs, S=Ss, T=0)` — (4.3.61)：|1⟩ 降 D、|2⟩ 升 D。
- `vd_reduced(D=Ds, H=Hs, S=Ss)` — 
- `vd_scan(Ds_=None, H: float=SCAN['H'], S: float=SCAN['S'])` — Fig. 4.3.7：{"D", "A'1", "A'2", "A'3"（3×3 由小到大）, "E", "B'"}。
- `vb_matrix(beta: float, H: float=SCAN['H'], S: float=SCAN['S'], T: float=SCAN['T'])` — 磁場 C4（本系列的模型，書上只給相關表 (4.3.73)）：赤道環 +x→+y→−x→−y（逆時針看 +z）的穿隧加 Peierls 相位：
- `vb_levels(beta: float, H: float=SCAN['H'], S: float=SCAN['S'], T: float=SCAN['T'])` — {C4 標記: 能量}：1_4 = H + 2S sin β、3_4 = H − 2S sin β、2_4 = H − 2S cos β、T1(0_4) = H − T、
- `v111_matrix(d=ds, H=Hs, S=Ss, T=0)` — (4.3.74)：沿 (111) 的電場，|1⟩|3⟩|5⟩（+z +x +y）降 d、|2⟩|4⟩|6⟩ 升 d。
- `v111_levels(d: float, H: float=SCAN['H'], S: float=SCAN['S'])` — 
- `v111_tunnel(Sp=sp.Symbol("S'"), S=Ss, H=Hs)` — (4.3.87)–(4.3.88)：順場—逆場的鄰居穿隧 S'、同側的鄰居 S（D3d）。{標記: 能量}（sympy）。
- `xf6_coords()` — 由群作用產生的 21 個座標 [(名稱, 原子, 向量)]：|lA⟩ = l·（+z 上的 ẑ）、|lB⟩ = l·（+z 上的 x̂）（＝ XF6_COORDS_LIST）。
- `xf6_force(k: float, j: float, b: float=0.0)` — 由彈簧幾何算 21×21 力矩陣（座標依 XF6_NAMES）：F–F 最近鄰彈簧 k（12 根）、X–F 徑向 j（6 根）、
- `xf6_mass(m: float=19, M: float=240)` — 
- `xf6_rep(g: str)` — O_h 元素 g 在 21 個座標上的表示 Γ(g)（Γ e_c = Σ_c' (g v_c)·v_c' e_c'；有快取，不要改回傳的陣列）。
- `xf6_rows(k, j, b, book: bool=False)` — (4.4.22) 的三列（⟨A|、⟨B|、⟨C1|；只有書上印的部分，其餘 None）。book=True 照印：⟨B|F|C1⟩ = −b/2、⟨C1|F|C1⟩ = 2(j+b)。
- `xf6_force_from_rows(k, j, b, book: bool=False)` — 只用 (4.4.22) 的三列＋對稱 F = Γ(g) F Γ(g)ᵀ ＋ Hermitian 重建整個 21×21（(4.4.2) 的想法）。
- `xf6_translation(axis: int=0)` — 沿 x_{axis+1} 的剛體平移（21 維）：(4.4.20) 的 |t⟩。
- `xf6_numeric(k: float, j: float=1.0, b: float=0.0, m: float=19, M: float=240)` — 21 個 ω²（由小到大；6 個 0）。
- `xf6_secular(k, j=1, b=0, m=19, M=240)` — (4.4.34)：S、P、λ± = ω²（更正後的 (4.4.33)）。sympy 或浮點都可以。
- `xf6_formulas(k, j=1, b=0, m=19, M=240)` — ω²：A1g (j+4k)/m、Eg (j+k)/m、T2g 2(k+b)/m、T2u (k+b)/m、T1u± ((4.4.23)–(4.4.26)、(4.4.34))。
- `xf6_t1u_matrix(k, j=1, b=0, m=19, M=240, book: bool=False)` — (4.4.33) 的 2×2 加速度矩陣（基底 |c⟩、|c'⟩）；book=True 照印右上 2(j+b)(M+6m)/3m。
- `xf6_c_vectors(m=19, M=240)` — (4.4.18)–(4.4.20) 的 x 夥伴（21 維）：|c⟩（軸向 F 1、赤道 F −½、中心 0）、|c'⟩（F 全部 M、中心 −6m）、|t⟩。
- `xf6_t1u_reduced(k, j=1.0, b=0.0, m=19, M=240)` — 用 21 維的 |c⟩、|c'⟩ 直接算 ⟨F⟩、⟨m⟩ 的 2×2（(4.4.29)–(4.4.32)）與 ⟨m⁻¹F⟩。
- `xf6_t1u_vectors(k=0.18, j=1.0, b=0.2, m=19, M=32)` — (4.4.37)–(4.4.41)、Fig. 4.4.7：{'+'/'-': dict(lam, eps, eps_p, center, axial, equatorial)}（m 歸一化 (4.4.40)）。
- `xf6_isotope(k=0.18, j=1.0, b=0.2, m=19, M1=32, M2=34)` — (4.4.36a)：ω±(M1)/ω±(M2)。
- `xf6_scan(kjs=None, m=19, M=240, b=0.0)` — Fig. 4.4.5（UF6：m = 19、M = 240、b = 0）／4.4.6（SF6：M = 32、b = 0.2）：ω 以 √(j/m) 為單位，k/j ∈ [0, 1]。
- `uf6_fit(m=19, M=240)` — 用 UF6 的 ω(T1u+)/ω(A1g) = 624/667 定 k/j（b = 0）；再用 A1g 定 j 的尺度，預測其他頻率（cm⁻¹）。
- `sf6_predict(kj=0.18, b=0.2, m=19, M=32)` — k/j = 0.18、b = 0.2 的 SF6：ω（單位 √(j/m)）與「用 A1g = 775 cm⁻¹ 定尺度」的預測。
- `xf6_isotypic(alpha: str)` — O_h irrep α 的同型投影 P^α = (ℓ/48) Σ χ^α(g) Γ(g)（21×21）。
- `xf6_classify()` — (4.4.3)–(4.4.9)：{軌道: O_h irrep 的頻率}（用 Γ 的跡）。
- `xf6_projected(alpha: str, j: int, k: int, orbit: str='A')` — P^α_{jk}|A⟩ 或 |B⟩（(4.4.10)–(4.4.17)，四方基底的 O_h irreps），歸一化的 21 維向量。
- `xf6_modes(k=0.18, j=1.0, b=0.2, m=19, M=32)` — 動畫用的簡正模（每種取第一個夥伴，座標依 XF6_NAMES）：[(標記, ω, 21 維位移)]；
- `ladder(nmax: int)` — 截斷的 a、a†（(nmax+1)×(nmax+1)，基底 |0⟩…|nmax⟩）：a|n⟩ = √n|n−1⟩（(4.4.62)）。
- `ladder_commutator(nmax: int)` — [a, a†] 在截斷空間：diag(1, …, 1, −nmax)（最後一格是截斷造成的）。
- `osc_levels(nmax: int, omega: float=1.0)` — H = ħω(a†a + ½) 的本徵值 (4.4.52)、(4.4.63)。
- `degeneracy(ell: int, n: int)` — (4.4.65)：d(n) = (ℓ+n−1)!/(n!(ℓ−1)!)。
- `degeneracy_count(ell: int, n: int)` — 數 ℓ 個非負整數和為 n 的組數（(4.4.64) 的另一種算法）。
- `harmonic_states(ell: int, n: int)` — ℓ 重簡併模的第 n 諧波的佔據數組（例：ℓ = 3、n = 2 → 6 組）。
- `state_count(Emax: float, modes: dict | None=None)` — 能量 ≤ Emax（cm⁻¹，零點能不計）的振動態數（含簡併，(4.4.66)）。
- `turning_point(n: int, ell: int, omega: float, hbar: float=1.0)` — (4.4.71)：q(TP) = [(2n + ℓ)ħ/ω]^{1/2}。
- `min_generators(G)` — 最少幾個生成元（習題 4.1.1；暴力搜尋 1、2、3 個）。
- `problem_411()` — 
- `problem_412()` — 
- `subgroups_of_order(G, order: int)` — G 中所有 order 階的子群（由兩個元素生成的；O 的 D2、D3、D4 都是）。
- `problem_413()` — D2（四階、三個 180° 轉）、D3、D4 子群與是否正規。
- `problem_421()` — 習題 4.2.1 的相關表（列：大群的 irrep；行：子群的 irrep；標記見 PROBLEM_421_LABELS）。
- `cage_symmetry(extra: list, field=None)` — 立方體八頂點（都有電荷）再加 extra 上的額外電荷（頂點的 ±1 座標），field = 磁場（軸向量）：
- `split_t1u(el: list)` — p 軌域（O_h 的 T1u）在子群 el 下分成幾個能階：回傳不可約成分的維度（由子群自己的類代數特徵標表）。
- `problem_423()` — {小題: (對稱群的階, 群名, T1u 分裂後的維度)}。
- `problem_431()` — 主 C3 誘導 0_3↑O（8 維）與 A1(D4)↑O（3 維）：{名: (在 1, r1, R1², R1, i1 的跡, 頻率)}。
- `problem_432()` — E(C3v)↑O_h（16 維）與 E_g(D4h)↑O_h（6 維）：{名: (在 1, r1, IR1², IR1, Ii1 的跡, 頻率)}。
- `problem_433()` — 立方體八井（局部 C3v）：{O_h irrep: (重數, 能量)}，最近鄰 −S（稜長 2）、次近鄰 −T（面對角 2√2）。
- `problem_434()` — 立方八面體 12 井（局部 C2v）：最近鄰 −S（距離 √2）、次近鄰 −T（距離 2）。
- `problem_435()` — Y ↓ D5、C5、D3、C3 的相關表（列 A, T1, T3, G, H）。
- `problem_436(F: float, axis=(0, 0, 1), kind: str='dipole', S: float=1.0, T: float=0.0)` — 立方八井加外場的能階（H = 0、穿隧 −S、−T）：dipole 位能 −F (n·r)/|r|、quadrupole −F (n·r)²/|r|²。
- `xy4_omega2(j: float, k: float, mX: float, mY: float)` — 四面體 XY4 的中心力模型（X–Y 鍵 j、Y–Y 稜 k）：{ν1 (A1), ν2 (E), ν3 (T2 高), ν4 (T2 低): ω²}（數值 15×15）。
- `xy4_fit(mol: str='CF4')` — 習題 4.4.1：兩個彈簧常數（j、k，單位 N/m）用最小平方擬合四個基頻（log ω）；回傳 j、k、預測與誤差。
- `problem_442()` — 立方 X8（每個頂點局部 C3v、3 個位移）與立方烷 C8H8 的模分類：{名: O_h irrep 的頻率（扣掉平移 T1u 與轉動 T1g）}。
- `xf6_coset_omega2(k, j=1.0, b=0.0, m=19)` — 陪集和 (4.3.46)、(4.4.23)–(4.4.29)：只用 ⟨A|F|·⟩、⟨B|F|·⟩ 兩列與四方 irreps 的一個分量求 ω²，以及 Q_AA、Q_AB、Q_BB。
- `summary()` — 
- `TRIG_T1_BOOK_R12()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 10 秒）。

## `combkit.sym_prod`

sym_prod.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 6 章

常數：`R_`、`EPS3N`、`_S_KOSHER`、`_S_TRIG`、`GROUP_OF`、`SITES`、`PAIRS`、`KRON_617X_BOOK`、`SIX_STATES`、`TWO_CARD`、`CG_619C_BOOK`、`CG_630_BOOK`、`CG_623_BOOK`、`CG_626_BOOK`、`CG_635A_BOOK`、`CG_635B_BOOK`、`CG_636_BOOK`、`CG_617_N_BOOK`、`PAIR_V`、`D2_TABLE_BOOK`、`TENSOR_COLS`、`TENSOR_KEYS`、`XS`、`POLY_639`、`MANDEL`、`MODULI_CARD`、`CUBIC_KITTEL`、`FIG639_A_BOOK`、`SIX_INV_ROWS_BOOK`、`SIX_CARD`、`WE_617_BOOK`、`FIG651`、`FIG651_RHO_BOOK`、`BOOK_CONST`、`RAYLEIGH_CARD`、`XF6_NU_LABELS`、`JT_CARD`、`D6H_IRREPS`、`LORENTZ_A`、`CARD_EXPECTED`、`KRON_617X_FULL_BOOK`

- `exact(x, maxden: int=96)` — 浮點（或複數）→ sympy：每個實、虛部寫成 ±√(p/q)（q ≤ maxden）；表格上的耦合係數都是這種數。
- `exact_mat(M)` — 
- **class `IrrepSet`** — 有限群 G 的一組不可約表示：elems（元素名）、D[α][g]（ℓ×ℓ complex）、names（irrep 的順序）。
  - `.dim(a)` — 
  - `.stack(a)` — 
  - `.chi(a)` — 
  - `.freq(chi)` — 
  - `.decompose(chi)` — 
  - `.product_chi(a, b)` — 
- `is_rep(G, D: dict)` — D[g]（數值）是否為 SN.FiniteGroup G 的表示（D(a)D(b) = D(ab)）。
- `oset(basis: str='tetragonal')` — O 的五個 irrep（A_1 A_2 E T_1 T_2）：basis = "tetragonal"（O⊃D4⊃D2，kosher T2 (4.2.16)＝ sym_oct 的 TETRA）、
- `ohset(basis: str='tetragonal')` — O_h 的十個 irrep（OH_IRREP_NAMES）；basis 同 oset。
- `d4set()` — D4 = {1, R3², R3, R3³, R1², R2², i3, i4}（O ⊃ D4 ⊃ D2 的中間一環）：A1、B1（x²−y²）、A2（z）、B2（xy）、E（x, y）（書上 (6.1.21) 的順序）。
- `d2set()` — D2 = {1, R1², R2², R3²}：A1、A2（R3² 為 +1，像 z）、B1（R1² 為 +1，像 x）、B2（R2² 為 +1，像 y）（(4.2.10)）。
- `d3set()` — D3 = {1, r1, r1², i2, i4, i5}（O ⊃ D3 ⊃ C2，C2 = {1, i4}）：A1、A2、E（三方基底 (v1, v2) 的塊，𝒟^E(r1) = rot(120°)）。
- `c3set()` — C3 = {1, r1, r1²}：0_3、1_3、2_3，𝒟^{k_3}(r1) = ε^k、ε = e^{2πi/3}（第 2–3 章的約定）。
- `tset()` — T = {1, r_k, r_k², R_k²}（12 個）：A、¹E、²E（ε、ε* 一對複共軛）、T（三維）。
- `qset()` — 四元數群 Q = {±1, ±i, ±j, ±k} 的二維 irrep（第三型的例子）。
- `c3vset()` — C3v（第 3 章的平面幾何 SN.C3V：σ₃ 沿 x 軸、σ₂ 在 60°、σ₁ 在 120°）：A1、A2（行列式）、E（2×2 本身＝ (3.3.7)）。
- `site_perm(g: str)` — O_h 元素 g 在六個頂點上的置換矩陣（⟨i|g|j⟩ = 1 若 g·site_j = site_i；＝ SO.site_perm 的數值版）。
- `outer_op(a: str, b: str)` — (6.1.4)：(a, b)|k, l⟩ = a|k⟩ b|l⟩ → 36×36 = kron(P(a), P(b))（列、行都依 PAIRS）。
- `kron_example()` — (6.1.7ₓ)：𝒟^{T1g}(r2) ⊗ 𝒟^{Eu}(i5)（四方基底；i5 是正常轉動，𝒟^{Eu}(i5) = 𝒟^E(i5)）。
- `outer_frequencies()` — 36 個兩粒子態在外積群 O_h^A × O_h^B（2304 個元素）下的分解：用 (6.1.9) 的 χ^{αβ}(a,b) = χ^α(a)χ^β(b)
- `outer_trace_check(n: int=60, seed: int=6)` — (6.1.9) 的直接驗證：隨機 n 對 (a, b)，Tr kron(𝒟^α(a), 𝒟^β(b)) = χ^α(a)χ^β(b)，且 (a,b)(a',b') = (aa', bb')（(6.1.6)）。
- `six_hamiltonian(H: float=0.0, S: float=1.0, T: float=0.0)` — (6.4.9)／Fig. 6.1.2 的六井 Hamilton 量：第一列 [H, −T, −S, −S, −S, −S]（對谷 −T、近鄰 −S）＝ SO.six_H(H, −T, −S)。
- `six_levels(H: float=0.0, S: float=1.0, T: float=0.0)` — 六井的三個能階（Fig. 6.1.2：T = 0 時 A1g = H − 4S、T1u = H、Eg = H + 2S）。
- `two_particle_levels(SA: float=1.0, SB: float=0.25, H: float=0.0)` — (6.1.10)–(6.1.11)：不相互作用時 ε^{αβ} = ε^α_A + ε^β_B，簡併 ℓ^αℓ^β → [(E, α, β, 簡併)]（由高到低，Fig. 6.1.2(a)）。
- `two_particle_numeric(SA: float=1.0, SB: float=0.25, H: float=0.0)` — 36×36 的 H_A ⊗ 1 + 1 ⊗ H_B 直接對角化（第二種算法）。
- `diag_frequencies(a: str, b: str, S: IrrepSet | None=None)` — (6.1.14)：f^γ(𝒟^{αβ}) = (1/°G)Σ χ^{γ*}χ^αχ^β（對角子群 O_h^{AB}）。
- `char_product_rows()` — (6.1.14) 下面的特徵標乘積表（O_h 十個類 1 r R² R i I Ir IR² IR Ii）：A1g、Eg、T1g、T2g 與 χ^{T1u}χ^{T1u}。
- `cg(S: IrrepSet, a: str, b: str, gammas=None)` — 耦合係數（(6.1.16)–(6.1.18b) 的投影法）→ (C, cols)：C 是 (ℓ^aℓ^b)×(Σℓ^γ) 的么正矩陣，列 = 乘積基底 |a i⟩|b j⟩
- `cg_norms(S: IrrepSet, a: str, b: str, gamma: str)` — (6.1.16c) 的歸一化常數 N^{γab}_{nkl}（依 (n, k, l) 字典序），= (ℓ^γ/°G)Σ 𝒟^{γ*}_{nn}𝒟^a_{kk}𝒟^b_{ll}。
- `cg_reduces(S: IrrepSet, a: str, b: str, C, cols)` — (6.1.19a)／(6.2.13)：對所有 g，C†(𝒟^a ⊗ 𝒟^b)(g)C = ⊕𝒟^γ(g)（依 cols 的塊）。
- `factorization_lemma(S: IrrepSet, a: str, b: str, gamma: str)` — (6.2.15)：(ℓ^γ/°G)Σ_g 𝒟^a_{ik}𝒟^b_{jl}𝒟^{γ*}_{mn} = Σ_ω C^{abγ_ω}_{ijm} C^{abγ_ω *}_{kln}（對所有 i j k l m n）。
- `cg_book(name: str)` — 依書上的鏈算出的表 → (C, cols)。name：
- `compare_table(C, cols, book)` — 把算出的 C 與書上的表逐格比：{"equal": 全同, "flips": 需要整體變號的 γ 份（書上的相位約定不同）, "max": 最大差}。
- `pair_V(v0: float=1.0, v1: float=0.5, v2: float=0.25, K: float=0.0, W: float=0.0)` — 兩粒子相互作用（36×36，列行依 PAIRS）：直接項（對角，依兩井距離 0、√2、2 取 v0、v1、v2）＋ K·近鄰交換 ＋ W·近鄰成對穿隧。
- `t1u_split(v=None)` — V 對 (T1u, T1u) 能階的一階分裂：用 (6.1.19c) 的耦合係數做出 |γ m⟩ = Σ C_{ij,m}|T1u i⟩|T1u j⟩，⟨γ m|V|γ m⟩ → {γ: 位移}。
- `t1u_split_numeric(v=None, lam: float=0.0)` — 第二種算法：V 限制到 9 維 (T1u)⊗(T1u) 子空間後直接對角化（不用耦合係數）；lam > 0 時改為整個 36×36 的
- `pair_V_invariant()` — V 與所有 (a, a) 對易（48 個），但不與 (1, R1) 對易（(6.1.12)：⟨1,1|V|1,1⟩ ≠ ⟨1,3|V|1,3⟩）。
- `d2_table()` — (6.1.20)：D2 一維 irreps 的乘法表 {(α, β): γ}（由特徵標相乘）。
- `solve_6129(book_d: bool=False)` — (6.1.27b)–(6.1.29)：把 T1⊗T1 耦合矩陣的結構（未知數 a…k）代入 (𝒟^{T1}⊗𝒟^{T1})(R1)·C = C·(A1⊕E⊕T1⊕T2)(R1)（kosher）
- `fs_indicator(S: IrrepSet, a: str)` — Frobenius–Schur 指標 (1/°G)Σχ(g²)（用 𝒟(g)² 的跡）：+1 第一型（實）、0 第二型（複共軛成對）、−1 第三型（必複但與共軛等價）。
- `irrep_type(S: IrrepSet, a: str)` — 第二種算法（共軛檢查）：χ 不全為實 → 第二型；否則解不變雙線性型 𝒟(g)ᵀJ𝒟(g) = J，J 對稱 → 第一型、反對稱 → 第三型。
- `sym_squares(S: IrrepSet, a: str)` — 對稱化與反對稱化平方 (6.2.7)：χ_±(g) = (χ(g)² ± χ(g²))/2 → (對稱的分解, 反對稱的分解)。
- `cg_symmetry(S: IrrepSet, a: str)` — 第二種算法：由 C^{aaγ} 的行看 C_{ij,m} = ±C_{ji,m}（(6.2.7a,b)）→ {γ: +1 對稱 / −1 反對稱}。
- `scalar_cg_check(S: IrrepSet)` — (6.2.2)、(6.2.5)–(6.2.6)：C^{α α* 0}_{ij·} = δ_ij/√ℓ（第一型 α* = α；第二型取共軛的那個 irrep）。
- `unitarity_check(S: IrrepSet, a: str, b: str)` — (6.2.12a,b)：C†C = 1（正交歸一）、CC† = 1（完備）；(6.2.14) 的反關係 𝒟^a⊗𝒟^b = C(⊕𝒟^γ)C†。
- `unit_tensors()` — (6.3.6)：𝐓^γ_m = Σ C^{T1uT1uγ}_{ij,m} x̂_i x̂_j（(6.1.19c) 的書上相位，非 kosher T2）→ {(γ, m): 3×3}。
- `tensor_components(Mt)` — (6.3.13a)：𝒯^γ_m = Σ C_{ij,m}𝒯_{ij}（任何 3×3 張量 → 九個對稱化分量）。
- `tensor_from_components(comp: dict)` — (6.3.13b)：𝒯_{ij} = Σ C_{ij,m}𝒯^γ_m。
- `tensor_transform_check()` — (6.3.5)：𝒱(g)𝐓^γ_m𝒱(g)ᵀ = Σ_{m'}𝒟^γ_{m'm}(g)𝐓^γ_{m'}（48 個 g；非 kosher 的 O_h irreps）。
- `poly_from_tensors()` — 把 (6.3.6) 的 x̂_i 換成座標 x_i：p^γ_m = Σ C_{ij,m}x_ix_j（(6.3.7)–(6.3.9)；T1g 的反對稱積為零）。
- `poly_transform_check()` — (6.3.8)：座標跟基向量一樣變：p_m(𝒱ᵀx) = Σ_{m'}𝒟_{m'm}(g)p_{m'}(x)（在隨機點上比，48 個 g）。
- `o3_labels()` — 1 ⊗ 1 = 0 ⊕ 1 ⊕ 2（O3）↓ O_h：0⁺ → A1g、1⁺ → T1g、2⁺ → Eg ⊕ T2g（RT.o3_to_o）；並用 χ¹(ω)² = χ⁰ + χ¹ + χ² 驗證。
- `stress_examples()` — Fig. 6.3.2(b)／(6.3.33) 的例子 → 每個例子的九個對稱化分量（只列非零的）：
- `face_forces(Tm)` — (6.3.12)：𝐅(j) = 𝒯·x̂_j（立方體第 j 面（+x_j 側）每單位面積受的力；相對的面是 −𝐅(j)）。
- `strain_of_displacement(s_func, h: float=1e-06)` — (6.3.14a)：𝒮_{ij} = ∂s_i/∂x_j（數值微分；均勻應變時與點無關）。
- `sym2_rep(Rm)` — 對稱二階張量（應力、應變）上的 6×6 表示：ρ(g)_{ab} = ⟨e_a, 𝒱 e_b 𝒱ᵀ⟩（Mandel 基底）。
- `SN_rot3(axis, deg)` — 
- `elastic_counts()` — 彈性常數的個數（暴力法）：沒有對稱 81 → 36 → 21；O_h：81 個分量裡 4 個不變量、36 → 3、21 → 3；
- `elastic_counts_chars()` — 第二種算法（張量定理）：Sym²(向量) 的頻率 f^γ → 36 個常數變 Σ(f^γ)²、21 個（互易）變 Σf^γ(f^γ+1)/2。
- `elastic_tensor(kA: float, kE: float, kT2: float, kT1: float=0.0)` — (6.3.26)：k_{hi,jk} = Σ_γΣ_m C_{hi,m}C_{jk,m}k^γ（(6.1.19c) 的係數）→ 3×3×3×3。
- `elastic_from_cartesian(K4)` — (6.3.28)（更正的 (6.3.28)ₓ）：k^γ = Σ C_{hi,m}C_{jk,m}k_{hi,jk}（取 m = 第一個）→ {A1, E, T2}；並回傳 Cartesian 的三個數。
- `moduli(kA, kE, kT2)` — (6.3.30)、(6.3.32)、(6.3.36)、(6.3.43)（精確值可給 Fraction）：B = k^{A1}/3、μ = k^{T2}/2、Y = 3k^Ek^{A1}/(2k^{A1}+k^E)、
- `moduli_numeric(kA: float, kE: float, kT2: float)` — 第二種算法：用 Cartesian 張量直接「做實驗」——等向應變 ε·1 → B = 𝒯_11/(ΔV/V)；剪應變 𝒮_12 = 𝒮_21 = ε → μ = 𝒯_12/(2ε)；
- `spring_model(j, k, l_)` — (6.3.41a–c)：SrTiO₃ 的剛性離子彈簧模型（單位 1/4a）：k^{A1} = j + 4k + 4l、k^E = j + k + l、k^{T2} = 2k + 2l。
- `spring_model_bonds(j: float, k: float, l_: float, a: float=1.0)` — 第二種算法（暴力法）：鈣鈦礦晶格（Ti 在原點、O 在 a·x̂_i、Sr 在 (a,a,a)，晶格常數 2a）每根中心力彈簧的能量 ½κ(d̂·𝒮·d)²
- `anisotropy(c11, c12, c44)` — (6.3.43)：A = (k_{11,11} − k_{11,22} − 2k_{12,12})/k_{12,12}（k_{11,11} = C11、k_{11,22} = C12、k_{12,12} = C44）。
- `c3v_field(g: str, kind: str='E')` — C3v 操作對場的 3×3 作用（SN.C3V 的平面幾何＋z 軸）：極向量 E → 𝒱(g)；軸向量 B（電流環）→ det(𝒱)·𝒱(g)。
- `c3v_field_irreps()` — E → A1 ⊕ E（(6.3.46)）、B → A2 ⊕ E（(6.3.48)）：用 3×3 的跡分解。
- `p3_12(kind: str='B', vec=(1, 0, 0))` — (6.3.47)：P^3_{12} = (2/6)Σ𝒟^E_{12}(g)* g 作用在場上；P^3_{12}B̂₁ = −B̂₂、P^3_{22}B̂₁ = B̂₁（更正的 √3/2）。
- `p3_12_book(s=sp.Rational(3, 2))` — 照印 (6.3.47)（把每個 √3/2 印成 3/2）逐項算 → B̂₂ 的係數（s = 3/2 → −3；s = √3/2 → −1）與 B̂₁ 的係數。
- `fig6311()` — Fig. 6.3.11 的四個式子（更正後）：rE₁ = −½E₁ + (√3/2)E₂、σ₃E₁ = E₁、σ₃E₂ = −E₂；rB₁ = −½B₁ + (√3/2)B₂、σ₃B₁ = −B₁、σ₃B₂ = B₂ …
- `striction_counts()` — (6.3.50)–(6.3.52) 的線性電致／磁致伸縮常數個數（暴力法：對稱 ij 的三階張量 ε_{ij,k} 在群下不變；軸向量多乘 det）：
- `striction_counts_chars()` — 第二種算法：Hom(場, 應變) 的維數 Σ_γ f^γ(Sym²V) f^γ(場)（C3v：Sym²V = 2A1 + 2E）。
- `six_invariants()` — (6.4.8)：V^{(A1g)_1} = |A1g⟩⟨A1g|、V^{(A1g)_2} = Σ|T1u j⟩⟨T1u j|/√3、V^{(A1g)_3} = Σ|Eg j⟩⟨Eg j|/√2（6×6）。
- `six_from_levels(eA, eT, eE)` — (6.4.10)：H = ε^{A1g}V₁ + √3ε^{T1u}V₂ + √2ε^{Eg}V₃ 的第一列 ⇒ H = (εA + 3εT + 2εE)/6、T = −(εA − 3εT + 2εE)/6、S = −(εA − εE)/6。
- `param_count(freqs: dict)` — 參數定理 (6.4.11a)：n = Σ f^α(f^α + 1)/2（Hermitian）。
- `commutant_count(mats, hermitian: bool=True)` — 第二種算法：與所有 Γ(g) 對易的（實對稱）矩陣個數 ＝ 解 Γ(g)XΓ(g)ᵀ = X 的線性方程（暴力法）。
- `xf6_param_counts()` — SF6 型 21 個座標：全部 Σf(f+1)/2 = 11（T1u 三次）；去掉平移、轉動後的 15 維振動空間 = 7（書 p. 496）；兩法。
- `six_dipole_elements()` — 六井的電偶極矩陣元 ⟨γ k|x_i|β j⟩（位置算子對角，井在單位距離、q = 1）：{(γ, k, i, β, j): 值}，只列非零。
- `we_check()` — Wigner–Eckart (6.4.15)–(6.4.17)：六井的數值偶極矩陣元 ＝ 耦合係數 × 約化矩陣元（四方基底，投影法的 C）。
- `we_intensities()` — Fig. 6.4.1 的強度比：極化 1（或 2）時 (T1u → Eg 1):(T1u → Eg 2) = 1:3、極化 3 時 1:0；(6.4.19) 的比 (R_A/R_E)²（本模型 = 1/2）。
- `we_complex_check()` — 〔注〕複數耦合係數時共軛在哪裡：T 群（12 個轉動）裡六井的 T1u → T、Eg → ¹E ⊕ ²E（複的一維）；偶極 T ⊗ T → ¹E 的 C 是複數。
- `lorentz_steady(ws: float, w0: float, gamma: float, A: float=1.0)` — (6.5.5e,f)：穩態 x_S = R cos(ω_s t − ρ)（刺激 A cos ω_s t）→ (R, ρ 度，0–180)。
- `lorentz_transient(ws, w0, gamma, A=1.0, x0=0.0, v0=0.0)` — (6.5.7)–(6.5.8)：暫態 x_T = a e^{−Γt}cos(ω_Γ t − α)，ω_Γ = √(ω0² − Γ²)。用 a cos α、a sin α 兩式（共振時書上的 tan 式是 0/0）；
- `lorentz_closed(t, ws, w0, gamma, A=1.0, x0=0.0, v0=0.0)` — (6.5.6) 的實部：x(t) = R cos(ω_s t − ρ) + a e^{−Γt}cos(ω_Γ t − α)。
- `lorentz_numeric(t, ws, w0, gamma, A=1.0, x0=0.0, v0=0.0)` — 第二種算法：scipy 直接積分 ẍ + 2Γẋ + ω0²x = A cos ω_s t（(6.5.1) 更正的 ẍ）。
- `fig651_panels()` — Fig. 6.5.1 的六幅：[(f_s, R, ρ 度, 最大 |x|（0–10 s）)]。
- `lorentz_trace(fs: float, tmax: float=10.0, n: int=801)` — 畫 Fig. 6.5.1 一幅用的取樣：(ts, 刺激 A cos ω_s t, 響應 x(t))（解析式；n 點）。
- `response(w, w0: float=1.0, gamma: float=0.1, A: float=2.0)` — (6.5.5c)：R e^{iρ} = A/(ω0² − ω² − 2iΓω)（Fig. 6.5.2：ω0 = 1、Γ = 0.1、A = 2 ⇒ S = A/2ω0 = 1）。
- `smith(delta: float, S: float=1.0, gamma: float=0.1)` — (6.5.9a)：近共振 R e^{iρ} ≈ S/(Δ − iΓ)。
- `smith_check(S: float=1.0, gamma: float=0.1)` — Smith 圓 (6.5.10)–(6.5.11)：Δ 變動時 S/(Δ−iΓ) 落在直徑 S/Γ（沿虛軸）的圓上；Γ 變動時落在直徑 S/Δ（沿實軸）的圓上；
- `fig652()` — Fig. 6.5.2（ω0 = 1、Γ = 0.1、A = 2）：Im R 的最大（≈ A/2Γω0 = 10）與 Re R 的極大、極小（在 ω0 ∓ Γ 附近，≈ ±5）、DC 響應 A/ω0² = 2。
- `epsilon(w, w0: float=1.0, wp: float=1.0, gamma: float=0.1, book: bool=False)` — (6.5.21)：ε(ω) = 1 + ω_p²/(ω0² − ω² − 2iΓω)；book=True 照印（虛部分子少了 2）。
- `fig654()` — Fig. 6.5.4（ω0 = ω_p = 1）：ε(0) = 2、截止 ω_c = √(ω0² + ω_p²) = √2；Γ = 0 時 Re ε < 0 恰在 (ω0, ω_c)；幾個點的值。
- `book_examples(eps0: float | None=None)` — §6.5B 的數值例（書上的常數）：α_e（氫，ω0 = 2×10¹⁶）、α_v（30 amu、ω0 = 9.4×10¹³＝500 cm⁻¹）與光頻時的 α_v、
- `rayleigh(blue=4.5e-07, red=7e-07, w0=2e+16)` — Rayleigh 散射 (6.5.32)–(6.5.33)：藍／紅的功率比 = (λ_r/λ_b)⁴（ω ≪ ω0）；含共振因子 |ω0² − ω²|⁻² 的完整 Lorentz 比；
- `radiation_integral()` — (6.5.28)：∫₀^π sin³θ dθ = 4/3（數值）⇒ P = q²ω⁴x²/(6πc³ε0)·cos²。
- `phasor_points(eps: float=0.2)` — Fig. 6.5.8 的響應羅盤（e^{−iωt} 約定：複振幅乘 e^{+iδ} 是落後 δ）：振子落後 ρ、下游響應場再落後 π/2：
- `xf6_selection()` — XF6 的選擇規則（特徵標法）：振動 = 全部 21 − 平移 T1u − 轉動 T1g = A1g + Eg + T2g + T2u + 2T1u；
- `xf6_selection_numeric(qF: float=-1.0, qX: float=6.0, d_bond: float=0.0001)` — 第二種算法（數值）：每個 O_h 同型子空間（SO.xf6_isotypic）裡，剛性離子偶極 𝐩 = Σq_a𝐮_a 的導數、
- `sym_power_chi(S: IrrepSet, a: str, n: int)` — Sym^n(𝒟^a) 的特徵標：Newton 遞迴 n·h_n = Σ_{k=1}^n p_k h_{n−k}，p_k(g) = Tr 𝒟(g)^k。
- `sym_power_chi_monomial(S: IrrepSet, a: str, n: int, seed: int=1)` — 第二種算法：在 n 次單項式基底上直接求 𝒟(g) 的作用矩陣（隨機點上最小平方），取跡。
- `overtone(alpha: str, n: int, S: IrrepSet | None=None)` — α 模的 n 次泛音（Sym^n）的 O_h 分解；u 的 irrep 時宇稱 = u^n（如 T2u：n = 2 → A1g + Eg + T2g、n = 3 → A1u + T1u + 2T2u）。
- `pascal(ell: int, n: int)` — (4.4.64)：ℓ 維 irrep 的 n 次泛音簡併度 C(ℓ + n − 1, n)（ℓ = 2：2, 3, 4…；ℓ = 3：3, 6, 10, 15…）。
- `combination(alpha: str, beta: str)` — 組合頻 α ⊗ β（α ≠ β，沒有對易對稱的刪減）：Fig. 6.6.1(b) 的 Eg ⊗ T1u = T1u ⊕ T2u。
- `henon_heiles_check()` — Hénon–Heiles V = x³ − 3xy²（＝ ρ³cos3φ）在 C3v（σ₃ 沿 x）下不變；3x²y − y³（ρ³ sin 3φ）屬 A2（鏡面變號）。
- `jt_matrix(x1: float, x2: float, j: float=1.0, k: float=1.0, r: float=0.0)` — (6.7.6a)：E⊗e 的 2×2 振電位能矩陣（電子基底 |E1⟩、|E2⟩）。
- `jt_sheets(rho, phi_deg, j: float=1.0, k: float=1.0, r: float=0.0)` — (6.7.8)：ε± = kρ² ± (j²ρ² + 2jrρ³cos3φ + r²ρ⁴)^{1/2}（r = 0 時 kρ² ± |j|ρ）。
- `jt_min(j: float=1.0, k: float=1.0)` — 墨西哥帽（r = 0）的谷底：ρ0 = |j|/2k、深度 ε−(ρ0) = −j²/4k。
- `renner_extrema(j: float=1.0, k: float=1.0, r: float=0.2)` — Renner 項（r ≠ 0）：jr > 0 時下層的三個極小在 φ = 0°, 120°, 240°（ρ = j/2(k−r)、ε = −j²/4(k−r)），鞍點在 60°, 180°, 300°
- `jt_eigvec_track(j: float=1.0, phis=None)` — (6.7.9)：r = 0 時沿谷底繞一圈追蹤下層電子態（數值本徵向量，連續選號）→ (φ 陣列, 向量陣列, 繞完的號 −1)；
- `jt_vibronic_levels(j: float=1.0, k: float=1.0, r: float=0.0, nmax: int=26, count: int=12)` — 振電能階（數值對角化）：H = ½(p₁² + p₂²) + kρ² + j(x₁σ_z − x₂σ_x) + r[(x₁² − x₂²)σ_z + 2x₁x₂σ_x]（ħ = 質量 = 1，
- `jt_radial_levels(l_: float=0.5, j: float=1.0, k: float=1.0, count: int=4, L: float=8.0, N: int=2400, book: bool=False)` — 第二種算法：r = 0 時固定半整數 l 的耦合徑向方程（更正的 (6.7.15)：位能在 ½[⋯] 外；交叉項 ∓il/2ρ²，號不影響能譜）有限差分 →
- `jt_effective(rho, l_: float, sheet: int=-1, j: float=1.0, k: float=1.0, book_fig: bool=False)` — Born–Oppenheimer 的有效位能：(l² + ¼)/2ρ² + kρ² ± jρ（更正後，ħ = 質量 = 1）；book_fig=True 是 Fig. 6.7.5 畫的
- `jt_effective_min(l_: float, sheet: int=-1, j: float=1.0, k: float=1.0, book_fig: bool=False)` — 
- `problem_61()` — 習題 6.1：子群鏈標記的 E⊗E 耦合係數表（投影法、書上的相位約定）：
- `problem_62()` — 習題 6.2（x 沿 C2 軸指向頂點、z 垂直平面的 D3）：(a)(b) 一次：(x, y) ∈ E、z ∈ A2；二次：x² + y²、z² ∈ A1，
- `problem_63()` — 習題 6.3（第 3 章 X3 模型，C3v，四個能階 A1、A2、E、E，六個態；ω 由 SN.six_dof_omega2）：從基態 A1 出發，
- `problem_64()` — 習題 6.4（X4 正方形模型，C4v ≅ D4，質點在 x、y 軸上，σ_v 過質點）：8 個平面座標 = A1 + A2 + B1 + B2 + 2E；
- `problem_65()` — 習題 6.5（Fig. 6.4.1 改成三方 O_h ↓ D3d）：T1u → A2u（沿 v3 = (111)/√3）⊕ Eu（v1, v2），A1g → A1g、Eg → Eg（不分裂）。
- `d6h_group()` — D6h 的 24 個元素：(k, s, i) ＝ I^i·R_z(60k)·C2x^s（s = 1 的是繞方位角 30k° 的 C2：k 偶 C2'（過原子）、k 奇 C2''）→ (名, 3×3)。
- `d6h_chi(ir: str, key: tuple)` — D6h 的特徵標（慣例：B1 對 C2'（過原子）為 +1）：key = (k, s, i)。
- `problem_66()` — 習題 6.6（D6h）：(a) 極向量 z ∈ A2u、(x, y) ∈ E1u；(b) 軸向量 R_z ∈ A2g、(R_x, R_y) ∈ E1g；
- `card_values()` — plan_ch06.md §6 每張卡上的數字與導讀第 6 章的數字（字串；小數照卡上位數）。
- `summary()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 15 秒）。

## `combkit.sym_rot`

sym_rot.py —— Harter《Principles of Symmetry, Dynamics, and Spectroscopy》（Wiley 1993）第 5 章

常數：`R_`、`IU`、`INV3`、`SIGMA_XY`、`I_Z`、`N_PARAMS`、`O2_LABELS`、`O2_MOVING_SIGMA_MINUS_BOOK`、`DINFH_CLASSES`、`DINFH_IRREPS`、`PAULI`、`EX538`、`EX539`、`FIG535_BOOK`、`_AX`、`GX`、`GY`、`GZ`、`GZ_PROJ_BOOK`、`SMALL_EPS`、`J2_JX_BOOK`、`D1_SYM`、`D32_SYM`、`D2_SYM`、`D2_BOOK`、`D2_TYPOS`、`ROTOR_CARD`、`O_CLASS_ANGLES`、`O_IRREPS`、`O3_O_BOOK`、`SPEC_LABELS`、`COSET_EULER`、`D2_COL_BOOK`、`T2_SUM_BOOK3`、`D_ORBITALS_BOOK`、`D_ORBITAL_POLY`、`X4`、`V4`、`D6_CLASS_ANGLES`、`D6_CLASS_SIZES`、`D6_IRREPS`、`O3_D6_BOOK`、`WHEEL_D6_EVEN`、`WHEEL_D6_ODD`、`WHEEL_O`、`WHEEL_O_START`、`DIMS`、`D6S_AXES`、`D6S`、`D6S_COLS`、`D6S_ROWS`、`_T574`、`D6S_TABLE_BOOK`、`D6S_CHARS`、`D6S_REPS`、`_FULL_CACHE`、`D6S_FREQ_BOOK`、`D6S_TRACE60_BOOK`、`D2S`、`D2S_TABLE_BOOK`、`D2S_DE_BOOK`、`Q_UNITS`、`OS`、`OS_TABLE_522_BOOK`、`OS_CHARS`、`OS_REPS`、`OS_FREQ_BOOK`、`WHEEL_D6_HALF`、`WHEEL_O_HALF`、`WHEEL_O_HALF_START`、`R4_BASIS`、`R4_REG_LZ_BOOK`、`R4_REG_KZ_BOOK`、`U3_G`、`U3_K`、`COULOMB_SHELLS`、`COULOMB_CLOSURES`、`OSC_SHELLS`、`OSC_CLOSURES`、`MAGIC_ATOMIC_BOOK`、`MAGIC_NUCLEAR_BOOK`、`CIRC_T`、`Y_CLASS_ANGLES`、`CARD_EXPECTED`

- `n_params(n: int)` — (5.1.12)：O_n 的參數個數 n(n−1)/2 = n² − n(n+1)/2。
- `o2_rot(phi_deg)` — (5.1.10a)：det = +1 的轉動。
- `o2_refl(phi_deg)` — (5.1.10b)：det = −1；＝ o2_rot(φ)·diag(1, −1) ＝ 對斜角 φ/2 的直線反射。
- `is_orthogonal(M)` — 
- `o2_standing(m: int)` — (5.2.3)：T·diag(e^{imφ}, e^{−imφ})·T† ＝ [[cos mφ, −sin mφ], [sin mφ, cos mφ]]（φ = PHI 符號）。
- `dinfh_chars(gamma_deg)` — Fig. 5.2.2 的 D∞h 特徵標（行 = DINFH_CLASSES）。Σ^± 的 ± 是對鏡面 IR(−γπγ) 的號；橫向 π 轉 R(−γπγ) 的號 ＝ (±)×宇稱
- `dinfh_trace_check(gamma_deg: float=40.0)` — 第二種算法：向量 (x, y, z) ＝ Σ⁺_u ⊕ Π_u、二次式 (x²+y², z², xz, yz, x²−y², xy) ＝ 2Σ⁺_g ⊕ Π_g ⊕ Δ_g 的跡。
- `rz(a)` — 
- `ry(b)` — 
- `rx(a)` — 
- `euler_matrix(a, b, c)` — (5.3.3)、(5.5.76)：R(αβγ) = R(α00)R(0β0)R(00γ) = R_z(α)R_y(β)R_z(γ)（度；最右邊先做）。
- `euler_matrix_sym()` — (5.5.76) 的 3×3（符號 α, β, γ）。
- `euler_from_matrix(Rm)` — 3×3 → (α, β, γ)（度），0 ≤ β ≤ 180（(5.3.2)）；β = 0 或 180 時只有 α ± γ 有意義，取 γ = 0。
- `euler_alt(a, b, c, sign: int=1)` — (5.3.1)：同一個位置的另一組 Euler 角 (α ± 180, −β, γ ± 180)（sign = ±1 選 ±）。
- `axis_vec(phi, theta)` — (5.5.5c)：Darboux 軸 (cos φ sin θ, sin φ sin θ, cos θ)。
- `darboux_matrix(phi, theta, omega)` — R[φθω]（軸角，Rodrigues ＝ 四元數 SN.quat）。
- `darboux_matrix_cranks(phi, theta, omega)` — §5.5A 的五個曲柄：R(φ00)R(0θ0)R(ω00)R(0−θ0)R(−φ00)。
- `darboux_from_matrix(Rm)` — 3×3 → [φ, θ, ω]（0 ≤ ω ≤ 180）。
- `su2_euler(a, b, c)` — (5.4.30)：𝒟^{1/2}(αβγ)（度）。
- `su2_darboux(phi, theta, omega)` — (5.5.1)／(5.5.4)：𝒟^{1/2}[φθω] = cos(ω/2)1 − i sin(ω/2) ω̂·σ。
- `quat_euler(a, b, c)` — Euler 角的四元數 q(z, α)·q(y, β)·q(z, γ)（SN.qmul）；SN.su2(q) ＝ su2_euler。
- `euler_to_darboux(a, b, c)` — (5.5.2) 的閉式：(α,β,γ) → [φ, θ, ω]（0 ≤ ω < 360）。
- `darboux_to_euler(phi, theta, omega)` — (5.5.2) 的反算（取 cos(β/2) ≥ 0、β 可為負：書上 Fig. 5.3.5 的分支）：γ − α = 180 − 2φ。
- `darboux_sequence(phi, theta, omegas)` — Fig. 5.3.5：R[φθω]|000⟩ 隨 ω 連續變化的 Euler 角（SU(2) 分支：β 在 360° < ω < 720° 為負，720° 才回到原狀）。
- `lab_frame_product(steps)` — 依序繞「實驗室」的軸轉：steps = [(軸, 度), …]，先做第一個 ⇒ R = R_n ⋯ R_1。
- `body_frame_product(steps)` — 依序繞「本體」（跟著轉的）軸轉 ⇒ R = R_1 ⋯ R_n（順序反過來，(5.3.5)）。
- `lab_op(Rm)` — 實驗室算符作用在位置態 |S⟩（S = 3×3）：|S⟩ → |R S⟩。
- `body_op(Rm)` — 本體算符 R̄：|S⟩ → |S R⁻¹⟩（＝ 把宇宙反轉，(5.3.7)）。
- `hamilton_example()` — H507 的例子：先 R_x(90°) 再 R_y(90°) ＝ 繞 (1, 1, −1)/√3 轉 120°（(5.5.10) ＝ 四元數 ＝ 半角弧，SN.hamilton_product）。
- `su2_product(t1, t2)` — (5.5.6)–(5.5.10)：R[t2]R[t1]（t = (φ, θ, ω)）→ (cos(ω''/2), ω̂'' sin(ω''/2))（先 t1）。
- `comm(A, B)` — 
- `gen_exp(axis, w_deg)` — e^{ω n̂·𝒢}（(5.4.6)；scipy expm）＝ 繞 n̂ 逆時針轉 ω。
- `gz_projectors()` — (5.4.7)：𝒢_z 的譜分解 {i: 𝒫₊, −i: 𝒫₋, 0: 𝒫₀}（SM.projectors）；e^{ω𝒢_z} = e^{iω}𝒫₊ + e^{−iω}𝒫₋ + 𝒫₀。
- `small_commutator(eps: float=SMALL_EPS)` — (5.4.8)–(5.4.10)、Fig. 5.4.2：R = e^{−ε𝒢y}e^{−ε𝒢x}e^{ε𝒢y}e^{ε𝒢x}（最右先做）。
- `jmat(j, exact: bool=False)` — (5.4.21)–(5.4.25)：基底 m = j, j−1, …, −j（列＝行同序），ħ = 1。→ {Jx, Jy, Jz, Jp, Jm, J2}（numpy 或 sympy）。
- `wigner_d(j, beta_deg)` — (5.4.45)（α = γ = 0）：d^j_{mn}(β)，列 m、行 n（都從 j 到 −j）。Wigner 的公式（整數運算：j + m 等都是整數）。
- `wigner_D(j, a, b, c)` — (5.4.45)：𝒟^j_{mn}(αβγ) = e^{−imα} d^j_{mn}(β) e^{−inγ}（度）。
- `D_expm(j, a, b, c)` — (5.4.28)：e^{−iαJz}e^{−iβJy}e^{−iγJz}（scipy expm）。
- `D_axis(j, axis, w_deg)` — 繞任意軸轉 ω 的 𝒟^j：e^{−iω n̂·J}。
- `D_schwinger(j, beta_deg)` — (5.4.36)–(5.4.44)：Schwinger 的兩個振子。把 (c a↑ + s a↓)^{j+n}(−s a↑ + c a↓)^{j−n} 展開成 a↑、a↓ 的多項式
- `d_sym(j)` — d^j(β) 的符號矩陣（Wigner 公式，cos(β/2)、sin(β/2)）。
- `d2_book_orthogonality(beta: float=0.7)` — 照印的 (5.4.35)：β = 0.7 rad 時 𝒟𝒟ᵀ 的 (1,1) 元（應為 1）。
- `d_symmetry_checks(j, b: float=37.0)` — 習題 5.4.6 的相位對稱（全部數值核對）。
- `sg_probs(j, n, beta_deg)` — |𝒟^j_{mn}(0β0)|²，m = j…−j（Fig. 5.5.2–5.5.5）。
- `sg_amps(j, n, beta_deg)` — 
- `mean_jz(j, n, beta_deg)` — (5.5.18)：Σ m |𝒟^j_{mn}|² ＝ n cos β。
- `cone_angle(j)` — (5.5.17)：(2 arccos(j/√(j(j+1))), 更正的近似 2 arccos(j/(j+½)), 照印的 arccos(j/(j+½)))（度）。
- `cone_theta(j, m)` — Fig. 5.4.4 的不確定角：cos θ^j_m = m/√(j(j+1))（度）。
- `rotor_ops(J)` — |J m n⟩ 空間（(2J+1)² 維，kron(m, n)）的實驗室 J_x,J_y,J_z 與本體（Jahn 反轉）J_x̄,J_ȳ,J_z̄ ＝ 1 ⊗ J_i*（(5.5.45)）。
- `rotor_levels(kind: str, J, **k)` — 數值對角化（(2J+1)² 維）→ [(能量, 簡併度)]。kind：
- `rotor_formula(kind: str, J, **k)` — (5.5.49)、(5.5.52)、(5.5.56)、(5.5.58) 的公式（含簡併度）。
- `asym_matrix_formula(J, A, B, C)` — (5.5.54)：H_asym 在 |J n⟩（n = J…−J）的矩陣（三對角跳 2）。
- `ylm(l, m, theta, phi)` — (5.5.71)：Y^l_m(θφ)（弧度；含 Condon–Shortley 號，與 (5.6.14) 相同）。
- `ylm_from_D(l, m, theta, phi)` — (5.5.65)：Y^l_m(θφ) = √((2l+1)/4π)·𝒟^{l*}_{m0}(φ, θ, ·)（Euler 角 α = φ、β = θ）。
- `legendre(l, x)` — 
- `addition_theorem(l, b_deg, a_deg, th_deg, ph_deg)` — (5.5.75)：(P_l(cos Θ), (4π/(2l+1)) Σ Y*_m(βα) Y_m(θφ))。
- `rotated_ylm(l, m, euler, th, ph)` — (5.5.70a)：Σ_{m'} 𝒟^l_{m'm}(αβγ) Y_{m'}(θφ) ＝ Y_m 在轉回去的點 R⁻¹r̂ 的值（弧度輸入 θ, φ）。
- `chi(j, w_deg)` — (5.6.4)：χ^j(ω) = sin((j+½)ω)/sin(ω/2)（ω = 0 → 2j+1）。
- `chi_exact(j, w_deg)` — Σ_m e^{−imω}（sympy 精確）。
- `o3_to_o(l)` — (5.6.1)、(5.6.5b)：𝒟^ℓ↓O 的頻率 (f^{A1}, f^{A2}, f^E, f^{T1}, f^{T2})（特徵標公式，精確）。
- `o3_to_o_matrices(l)` — 第二種算法：24 個 O 元素的 𝒟^ℓ(g) ＝ e^{−iω n̂·J}（軸角 SO.O_TURNS），P 算子的跡。
- `o3_label(l, parity=None)` — 𝒟^ℓ↓O_h 的 TeX（宇稱 g/u 預設跟 ℓ 的奇偶，如單電子軌域）。
- `invariant_count(k: int)` — (5.6.20)–(5.6.21)：2k 次立方不變多項式的個數，三種算法：
- `d_orbital_states()` — (5.6.9)–(5.6.13)：只對六個陪集代表求和（係數 E：1/3，T₂：1/2），再歸一化。→ {"E1","E2","T2_1","T2_2","T2_3": 5 維向量（m = 2…−2）}
- `d_orbital_projector_full(label: str)` — 第二種算法：對 24 個元素求和 P = (ℓ/24)Σ 𝒟^α_{tt}(g)* 𝒟²(g)，作用在 |2 2⟩。
- `v4_identities()` — (5.6.28) 方括號的兩種寫法、(5.6.29) 的球座標寫法（sympy 化簡為 0）。
- `oct_potential(r, Q: float=1.0, R: float=1.0)` — (5.6.23)：Σ Q/|R r̂_j − r|（六個電荷在 ±x, ±y, ±z）。
- `oct_multipole(r, Q: float=1.0, R: float=1.0)` — (5.6.28)：6Q/R + (35Q/4R⁵)(x⁴ + y⁴ + z⁴ − 3r⁴/5)。
- `multipole_coeffs(Q: float=1.0, R: float=1.0)` — (5.6.25b)：x^ℓ_m = Σ_j (q_j / r_j^{2ℓ+1}) X^{ℓ*}_m(r_j)（ℓ ≤ 4，用 Y^ℓ_m 數值算）。
- `sf6_levels(J, t: float=1.0, B: float=0.0)` — B J² + t(J_x̄⁴ + J_ȳ⁴ + J_z̄⁴ − (3/5)J⁴) 在 |J n⟩（2J+1 維）的能階 → [(能量, O 標記, 簡併度)]，由低到高
- `sf6_clusters(J: int, gap_ratio: float=0.01)` — 大 J 的能階叢（相鄰能階差 < gap_ratio × 總寬度就算同一叢）→ [[標記…], …]（由高到低）。J = 20：
- `o3_to_d6(l, parity: int=1)` — (5.6.37)：𝒟^ℓ↓D6 的頻率（欄 A1 A2 E2 B1 B2 E1）。
- `o3_to_d6_mrule(l)` — 第二種算法（(5.6.37) 的註）：m = 0 → A1（ℓ 偶）或 A2（ℓ 奇）；±m 對依 |m| mod 6：±1,±5 → E1、±2,±4 → E2、3 → B1+B2、0 → A1+A2。
- `d180(N)` — (5.6.39)–(5.6.40)：𝒟^N(0 180° 0) = δ_{n,−m}(−1)^{N+m}。
- `wheel_read(slots, start: int, step: int, dim: int)` — 從 slots[start] 起，每次走 step 格，讀到維度和 = dim 為止（不剛好就 AssertionError）。
- `wheel_d6(J: int)` — 
- `wheel_o(J: int)` — 
- **class `SpinAlgebra`** — 有限轉動群的自旋 ½「射線代數」：元素用 SU(2) 矩陣代表（名稱 → 2×2），乘積 ab = ±c（先 b 後 a）。
  - `.find(M)` — 
  - `.mul(a, b)` — 
  - `.inv(a)` — 
  - `.signed_tex(sg, g)` — 
  - `.table(rows, cols)` — 
  - `.ga_mul(x: dict, y: dict)` — 
  - `.conj_sum(h)` — 
  - `.class_ops()` — 全對易算子：每個 h 的 Σ g h g⁻¹ 除以 °N_h（＝ 每個不同的共軛元素出現一次），去掉 0 與重複。
  - `.trace_char(j)` — χ^j(g) ＝ 𝒟^j(g) 的跡（g 的 SU(2) 代表 → 轉角 ω ∈ [0, 720)）。
- `d6s_table_check()` — (5.7.4) 的 121 個格子逐一和 SU(2) 乘積比對 → 相符的個數。
- `d6s_class_algebra()` — (5.7.8)–(5.7.13)：c2 = h² − h⁴、c3 = h − h⁵；c2² = 2 − c2、c3² = 2 + c2、c2c3 = c3；c3³ = 3c3；冪等元與特徵標。
- `ray_freq(alg: SpinAlgebra, chars_full: dict, j)` — f^α(𝒟^j↓G) = (1/°G)Σ_g χ^α(g)* χ^j(g)（g 的 SU(2) 代表；乘積與號無關）。
- `cn_spin_chars(n: int)` — Fig. 5.7.5：C_n 的自旋表示 χ(h) = −1 的 n 次方根 e^{iπ(2k+1)/n}（h = 繞 z 轉 360°/n 的 SU(2) 代表，hⁿ = −1）。
- `dn_spin_algebra(n: int)` — 
- `dn_ray_irreps(n: int)` — D_n 射線不可約表示的維度：由類代數的維數（非零的全對易算子個數）與 Σℓ² = 2n 定出；書上：n 偶 n/2 個二維，
- `os_table_check()` — 
- `os_class_algebra()` — (5.7.21)、(5.7.23)：c_r = Σr_k − Σr_k²、c_R = ΣR_k − ΣR_k³；c_r² = 8 + 2c_r、c_rc_R = 4c_R、c_R² = 6 + 3c_r、c_R³ = 18c_R。
- `half_correlation(group: str, j)` — 𝒟^j↓D6（E1, E2, E3）或 ↓O（E1, E2, G）的頻率（半整數 j，SU(2) 跡）。
- `half_d6_mrule(j)` — 半整數 j 的 D6：±m 對依 m mod 6：±½ → E1、±3/2 → E3、±5/2 → E2（Fig. 5.7.3 的輪；圖上 E2、E3 的名字對調了）。
- `wheel_d6_half(j)` — 
- `wheel_o_half(j)` — 
- `r4_generators()` — (5.8.15)：4×4 的 L_x, L_y, L_z（前三維）與 K_x, K_y, K_z（第四維）。
- `r4_regular(name: str)` — (5.8.17b)、(5.8.18)：Lie 代數的正則表示 [X, e_b] = Σ_a 𝓡_ab e_a（基底 L_x, K_x, L_y, K_y, L_z, K_z）。
- `r4_ladders()` — (5.8.20a) 更正：M± = L± + K±、N± = L± − K±（L± = L_x ± iL_y）。
- `coulomb_shell(n: int)` — n 殼層（n² 維）：M = J^{(j)} ⊗ 1、N = 1 ⊗ J^{(j)}，j = (n−1)/2；L = M + N、K = M − N（(5.8.21)）。
- `coulomb_energy(n: int, k: float=-1.0, m: float=1.0)` — (5.8.33)–(5.8.34)：−(k²m/2)/(2M² + 2N² + 1) ＝ −(k²m/2)/n²（本系列取 k = −1、m = 1：−1/(2n²)）。
- `coulomb_weights(jm, jn)` — Fig. 5.8.3：𝒟^{j_m j_n} 的權 (L_z, K_z) = (m + n, m − n)。
- `kepler_check(r0=(1.0, 0.0, 0.0), v0=(0.0, 1.2, 0.3), k: float=-1.0, m: float=1.0)` — 古典 Kepler 軌道（數值積分）：ε = r/r − L×p/km 守恆、ε² = 1 + (2H/k²m)L²（(5.8.3)–(5.8.5)）。
- `runge_lenz_poisson()` — 古典 Poisson 括號（sympy）：{ε₁, ε₂} = −(2H/k²m)L₃ ⇒ 量子 [ε₁, ε₂] = iħ{ , } = −i(2H/k²m)L₃（(5.8.10) 第二行）；
- `quantum_xr_p(point=(0.3, -0.4, 0.5))` — (5.8.8) 更正：[x_i/r, p_j]ψ = i(δ_ij − x_i x_j/r²)ψ/r（p = −i∂；用試驗函數 ψ 在一點數值核對 i = 1、j = 1, 2）。
- `quantum_eps_denominator(point=(0.3, -0.4, 0.5))` — (5.8.6a)（對稱化）＝ (5.8.6b) 的分母 km 版；照印的 2km 不相等。→ (km 版誤差, 2km 版誤差)。
- `E3(i: int, j: int)` — (5.8.39)：基本矩陣 𝓔_ij（1 起算）。
- `su3_cartan(N: int=6)` — (5.8.41)–(5.8.42)：H1 = (E11 − E33)/√(2N)、H2 = (E11 − 2E22 + E33)/√(6N)、E_{±1} = E12/√6, E21/√6 …（N = 6）。
- `su3_weights(n1: int, n2: int)` — 𝒟^{n1 n2} 的權重數：列舉形狀 (n1+n2, n2) 的半標準 Young 表（數字 1–3），權 (a1 − a3, (a1 − 2a2 + a3)/√3)
- `su3_dim(n1: int, n2: int)` — 
- `osc_shell(N: int)` — 三維振子第 N 殼層：(n1, n2, n3) 的個數 ＝ (N+1)(N+2)/2、ℓ = N, N−2, …（Σ(2ℓ+1) 相同）、能量 N + 3/2（ħω = 1）。
- `osc_fock_levels(nmax: int=3)` — (5.8.35a)：截斷的 Fock 空間裡 Σa†a + 3/2 的本徵值與重數（ω = 1）。
- `appendix_e_check(j=Fraction(5, 2))` — (E.2)–(E.6)、(E.17)：[J3, J±] = ±J±、[J+, J−] = 2J3、J² = J+J− + J3² − J3 = J−J+ + J3² + J3、⟨m+1|J+|m⟩ = √(j(j+1) − m(m+1))。
- `problem_511(n: int=4, seed: int=5)` — 習題 5.1.1：正交矩陣的列也正交歸一；φ₁·Ψ₁ = (R²)₁₁、φ₁·Ψ₂ = (R²)₁₂；么正矩陣：⟨φ₁|φ₂⟩ = 0、⟨φ₁|Ψ_k⟩ = (R̄R)₁ₖ。
- `problem_512(phi: float=120.0)` — 習題 5.1.2：跡 = 1 + 2cos φ（任何軸）、非正常 −(1 + 2cos φ)；非正交基底 a = e1+e2、b = e1+e3、c = e2+e3 的
- `problem_513()` — 習題 5.1.3：晶格基底裡的跡是整數 ⇒ 1 + 2cos ω ∈ ℤ ⇒ ω = 0, 60, 90, 120, 180°（非正常 −(1 + 2cos ω) 同理）。
- `dn_chars(n: int)` — D_n 的特徵標表（類：E、C_n^p（p = 1…⌊n/2⌋）、C2'、[C2'']）→ (類 [(種類, 角度或撇, 大小)], {標記: 特徵標})。
- `dinfh_to_dnh(n: int)` — 習題 5.2.1：D∞h ↓ D_nh（n 偶：D_n × C_i，標記加 g/u；n 奇：D_n × σ_h，標記加 ′/″）。特徵標算法：D_n 部分用轉角與
- `dinfh_to_dnh_mrule(n: int)` — 第二種算法：Λ = |m|：m ≡ 0 (mod n) → A1+A2（m = 0：Σ⁺ → A1、Σ⁻ → A2，依 C2' 的號 ±·宇稱）；m ≡ n/2 → B1+B2；其他 → E_{min(m, n−m)}；
- `problem_533()` — 習題 5.3.3：(b) 給定矩陣的 Euler 角；(c) 軸角（兩種算法：矩陣的軸角、(5.5.2) 的閉式）。
- `problem_542(a=20, b=50, c=70)` — 習題 5.4.2：𝒟¹(αβγ) ＝ T†R(αβγ)T（T 的行是圓向量 CIRC_T）；跡相等 1 + 2cos ω。
- `D_from_u(j, u)` — 習題 5.4.5：由 𝒟^{1/2}(u) = [[u11, u12], [u21, u22]] 得
- `problem_553()` — 習題 5.5.3(b)：R1²、r1、i1、R1 的 Euler 角 (α, β, γ) 與軸角 [φ, θ, ω]（度；β = 180° 時只有 α − γ 有意義，取 γ = 0）。
- `problem_554()` — 習題 5.5.4：A = 1/2I₁、B = 1/2I₂、C = 1/2I₃；J = 1、2 的能階（(5.5.54) 的矩陣；rotor_levels 數值核對）。
- `problem_555()` — 習題 5.5.5：B(J_x̄² + J_ȳ²) + C J_z̄² 用 (5.5.83b) 換成 B[(J_α − cos β J_γ)²/sin²β + J_β²] + C J_γ²（sympy：差為 0）。
- `problem_561(lmax: int=12)` — 習題 5.6.1：O(3) ↓ D6 到 ℓ = 12（欄 A1 A2 E2 B1 B2 E1）。
- `o3_to_y(l)` — 習題 5.6.2(b)：𝒟^ℓ ↓ I（欄 A, T1, T3, G, H；SO.Y_CHARS）；I_h 的宇稱照 ℓ±。
- `o3_to_y_mats(l)` — 第二種算法：60 個轉動的 𝒟^ℓ 跡（數值，軸角公式 χ^ℓ(ω)）配 Y 的類（轉角 → 類；72° 與 144° 分開）。
- `o3_to_td(l, parity: int)` — 習題 5.6.2(a)：𝒟^{ℓ±} ↓ T_d（欄 A1 A2 E T1 T2）：+ 與 ↓O 相同；− 時 A1↔A2、T1↔T2（T_d 的 IR、Ii 類變號）。
- `o3_to_td_chars(l, parity: int)` — 第二種算法：直接用 T_d 的元素（SO.TD：I·R_k、I·i_k）與 χ^{ℓ±}(IR) = ±χ^ℓ(R)，配 O 的特徵標（T_d ≅ O）。
- `card_values()` — plan_ch05.md §6 每張卡上的數字（字串；小數照卡上位數）。
- `summary()` — 
- `check()` — 每個上卡片／講義的數字用兩種以上獨立算法核對（約 15 秒）。

## `combkit.symosc`

symosc.py —— 對稱群與耦合振盪子：D3、破缺到 C3／C2、四元群 Q8 的數學（numpy；不 import manim）

常數：`OMEGA`、`D3_NAMES`、`_D3_PERM`、`D3_CLASSES`、`D3_CHARS`、`D3_SIZES`、`P_CYC`、`Q8_NAMES`、`_U`、`Q8_REP2`、`Q8_CLASSES`、`Q8_SIZES`、`Q8_CHARS`

- `perm_matrix(p)` — p[j] = 擺 j 被送到的位置：D[p[j], j] = 1。
- `d3_matrix(name: str)` — 
- `K_d3(k0: float, k: float)` — 
- `d3_chars_of(D)` — 
- `decompose(chars, table: dict, sizes)` — 
- `class_of_d3(name: str)` — 
- `projector_d3(irrep: str)` — 
- `class_sum_values(table: dict, sizes)` — 
- `d3_freqs(k0: float, k: float)` — 
- `beat(t: float, k0: float, k: float)` — 
- `v_m(m: int)` — 
- `H_flux(a: float, b: float, phi: float)` — 
- `flux_levels(a: float, b: float, phi: float)` — 
- `gyro_freqs(k0: float, k: float, g: float)` — m = 1：ω² + 2g sinθ_m ω − κ_m = 0 的正根（x ∝ v_m e^{iωt}）。
- `gyro_freqs_numeric(k0: float, k: float, g: float)` — 一階化 [x, v]' = [[0, I], [−K, −G]] [x, v] 的本徵值 iω → 排序後的正 ω。
- `precession(t: float, delta: float)` — cos(θ_j + ω₊t) + cos(θ_j − ω₋t) = 2cos(ω̄t)·cos(θ_j − δt/2)：m = +1 是 v₁e^{iω₊t}、m = −1 是 v₋₁e^{iω₋t}。
- `K_c2(k0: float, k: float, delta: float)` — 
- `c2_levels(k0: float, k: float, delta: float)` — 
- `q8_mul(a: str, b: str)` — 
- `q8_inv(a: str)` — 
- `L(g: str)` — (L(g)x)_h = x_{g⁻¹h}：把擺 h 送到 g·h。
- `R(g: str)` — (R(g)x)_h = x_{h g}。
- `K_q8(k0: float, ki: float, kj: float, kk: float)` — 
- `q8_class(n: str)` — 
- `q8_levels(k0: float, ki: float, kj: float, kk: float)` — 
- `q8_irrep_vector(irrep: str)` — 一維表示：x_h = χ_α(h)。
- `q8_E_basis()` — [(名稱, 向量, m = L(i) 的本徵值, n = R(i) 的本徵值)]；v3 = L(j)v1、v4 = L(j)v2。
- `frobenius_schur(chars, classes: dict, sizes, square)` — (1/|G|) Σ_g χ(g²)。
- `check()` — 
- `main()` — `python utils/symosc.py`、`python -m symosc` 與安裝後的 `symosc-check` 都跑這裡。

## `combkit.tikz_trees`

tikz_trees.py —— 平面樹（有根、子樹有順序）的 TikZ 畫法（純 Python，不需要 manim）

- `plane_tree_layout(pt, dx: float=0.42, dy: float=0.5, center: bool=True)` — 回傳 (preorder, pos)：preorder[v] = (父節點, 深度)，pos[v] = (x, y)；頂點按前序（第一次到訪）編號。
- `plane_tree_prims(pt, dx: float=0.42, dy: float=0.5, labels=True, color: str='cA', label_style: str='font=\\tiny,gray!70!black')` — 平面樹的圖元。labels=True 標前序編號 0,1,2,…；也可以給一個字串或串列（例如 "ijk"）當頂點名稱；False 不標。

## `combkit.walk1d`

walk1d.py —— 《Elements of the Random Walk》第 1 章的數學（純 Python，不 import manim）

常數：`MAIN_COINS`、`MAIN_STEPS`、`MAIN_POS`、`BOOK_LN10`

- `coins_to_steps(coins: str)` — 
- `positions(steps)` — 
- `walk_count(N: int, d: int)` — 
- `walk_prob(N: int, d: int)` — 
- `distribution(N: int)` — 
- `biased_prob(N: int, n: int, p)` — 
- `all_walks(N: int)` — 
- `sequences_with_heads(N: int, n: int)` — 
- `ensemble(M: int, N: int, seed: int)` — 
- `finals(walks)` — 
- `mean(xs)` — 
- `mean_square(xs)` — 
- `histogram(xs)` — 
- `moment_exact(N: int, k: int)` — 
- `moment_formula(N: int, k: int)` — 
- `moment_via_cosh(N: int, k: int)` — (1.24)：(−2)^k · k! · [y^k] cosh(y/2)^N，冪級數一次一次相乘（不經過分布加總）。
- `mean_heads(N: int)` — 
- `pair_grid_counts(N: int)` — 
- `fourth_moment_pairings(N: int)` — Δ_iΔ_jΔ_kΔ_l 的平均不為 0：四個指標全相同（N 組）或兩兩成對（3·N(N−1) 組）。
- `ln_factorial(n: int)` — 
- `stirling(n: float)` — 
- `ln_stirling(n: float)` — 
- `stirling_ratio(n: int)` — 
- `gamma_integrand(w: float, N: int)` — 
- `gamma_gauss(w: float, N: int)` — 
- `gauss_density(u: float, N: int)` — 
- `gauss_prob(N: int, d: int)` — 
- `ln_binom_approx(N: int, m: float)` — 
- `pascal_rows(Nmax: int)` — 
- `diffusion_step(P: dict)` — 一步：每個人各以 ½ 往左右走。P 的值用 Fraction（或 int）就是精確的。
- `second_difference(f, y)` — 
- `heat_kernel(y: float, N: float, l: float=1.0)` — 
- `chi(q: float)` — 
- `chi_power_coeffs(N: int)` — 把 (e^{iq}+e^{−iq})^N 一次一次乘開，記錄 e^{iqd} 的係數。
- `fourier_count(N: int, d: int, M: int | None=None)` — 
- `saddle_q0(N: int, d: float)` — 
- `saddle_both(N: int, d: int)` — 
- `gz0_coeffs(n: int)` — 
- `d2_gen_coeffs(n: int)` — 
- `d4_gen_coeffs(n: int)` — 
- `cauchy_coeff(f, N: int, r: float, M: int=256)` — 
- `walk_saddle_z(N: int, q: float)` — 
- `sd_exact(N: int, a: float, zc: float=1.0)` — 
- `sd_approx(N: int, a: float, zc: float=1.0)` — 
- `sd_book(N: int, a: float, zc: float=1.0)` — 
- `sd_saddle_exact(N: int, a: float, zc: float=1.0)` — 
- `sd_saddle_book(N: int, a: float, zc: float=1.0)` — 
- `sd_p(z, N: int, a: float, zc: float=1.0)` — 
- `sd_p2(z, N: int, a: float, zc: float=1.0)` — 
- `sd_full(N: int, a: float, zc: float=1.0)` — 
- `check()` — 

## `combkit.walk_bc`

walk_bc.py —— 《Elements of the Random Walk》第 4 章 Boundary conditions, steady state, and the electrostatic analogy 的數學

常數：`PI`、`PREFACTOR`、`CELL_EXAMPLE`、`BP_EXAMPLE`

- `walks(N: int, x0: int)` — 從 x0 出發的全部 2^N 條一維走法（位置 list，長 N+1）。
- `touches_wall(path)` — 這條路線有沒有踩到牆（位置 ≤ 0）。
- `survivors_brute(N: int, x0: int)` — {終點 x: 從 x0 走 N 步、從沒踩到 0 的走法數}（暴力列舉）。
- `survivors_image(N: int, x0: int)` — (4.1) 的一維格子版：C(N; x−x0) − C(N; x+x0)（鏡像點 −x0），x ≥ 1。
- `survivors_dp(N: int, x0: int)` — 逐步遞迴（每步把踩到 0 的刪掉；第三種算法）。
- `reflect_first(path: list[int])` — 把第一次踩到 0 之前的那一段對 0 鏡射（Fig. 4.1）：碰過牆的路線 ↔ 從 −x0 出發的路線。
- `reflection_pairs(N: int, x0: int)` — [(碰過牆、終點 > 0 的路線, 它的鏡像路線)]：鏡像路線從 −x0 出發、終點相同。
- `survival_prob(N: int, x0: int)` — N 步都沒踩到 0 的機率（精確）= P(1−x0 ≤ S_N ≤ x0)，S_N 是無限制的位移。
- `survival_asym(N: float, x0: float)` — 一維 (4.5) 的類比：x0·√(2/(πN))（= 2x0 × 高斯密度 (2πN)^{−1/2}）。
- `survival_3d(N: int, x0: int)` — sc 格子、吸收平面 x = 0、從 (x0, 0, 0) 出發：N 步都沒踩到平面的機率。
- `survival_3d_dp(N: int, x0: int)` — 同上，用「懶惰的一維漫步」逐步遞迴（每步 2/3 不動、1/6 向左、1/6 向右；第二種算法）。
- `survival_3d_asym(N: float, x0: float)` — 書上 (4.5) 之後：√(6x0²/(πN))。
- `wall_reduction_1d(N: int)` — 一維 (4.6)–(4.8)：Σ_{x0≥1} (1 − 存活機率)＝因為牆而少掉的走法數 ÷ 2^N（精確）。
- `wall_reduction_brute(N: int)` — 同上的第二種算法：等於 N 步漫步的最大值的期望值 E[max_k S_k]（暴力列舉）。
- `wall_reduction_asym(N: float, d: int=1)` — (4.8)：一維 √(2N/π)；三維（sc，x 方向的變異數 N/3）√(2N/(3π))。
- `box_dp(N: int, L: int, x0: int)` — {x: 從 x0 走 N 步到 x、從沒碰到 0 或 L 的走法數}（轉移矩陣＝逐步遞迴）。
- `box_images(N: int, L: int, x0: int, x: int)` — Ex. 4.3：無窮多個鏡像點 2kL + x0（正）與 2kL − x0（負）。
- `box_eigs(L: int)` — 內點 1..L−1 的轉移矩陣（鄰居的和）的特徵值 2cos(mπ/L)，m = 1..L−1。
- `box_mode(L: int, m: int)` — 正規化的正弦模態 √(2/L)·sin(mπx/L)，x = 1..L−1（(4.9)(4.10) 的格子版）。
- `box_eig_residual(L: int)` — max |T v_m − λ_m v_m|（核對特徵向量；Ex. 4.4 的格子版）。
- `box_gf(z: float, L: int)` — 兩牆之間的生成函數矩陣 G(z; x1, x2) = [(I − zT)^{−1}]，x1, x2 = 1..L−1（高斯消去；索引 0 ↔ x = 1）。
- `box_amplitude(z: float, L: int, m: int, x2: int)` — (4.12)(4.13) 的 g(z; q)：把 G(z; x1, x2) 對 x1 投影到 √(2/L) sin(qx1)，q = mπ/L（數值）。
- `box_amplitude_formula(z: float, L: int, m: int, x2: int)` — 更正後的 (4.13)：√(2/L) sin(q x2)/(1 − 2z cos q)（書上印成 sin q x1）。
- `box_image_points(L: int, x0: int, K: int=1)` — Ex. 4.3 的鏡像點（k = −K..K）：(正的 2kL + x0 們, 負的 2kL − x0 們)。
- `box_modes(N: int, L: int, x0: int, x: int)` — (4.15) 的格子精確版：Σ_m (2/L) sin(mπx0/L) sin(mπx/L) (2cos(mπ/L))^N。
- `box_total(N: int, L: int)` — 從 1..L−1 任一點出發、終點不限、N 步都沒碰牆的走法總數（逐步遞迴，精確整數）。
- `box_total_modes(N: int, L: int)` — 同上 ÷ 2^N，用格子模態：Σ_{m 奇} (2/L) cot²(mπ/2L) cos^N(mπ/L)（第二種算法）。
- `box_total_cont(N: float, L: float, mmax: int=2000)` — (4.19) ÷ 2^N：(8L/π²) Σ_m e^{−N((2m+1)π/L)²/2}/(2m+1)²（連續極限，Δ = 1）。
- `box_long(N: float, L: float, book: bool=False)` — (4.24) ÷ 2^N：(8L/π²)e^{−π²N/(2L²)}；book=True 用書上印的 e^{−π²N/L²}。
- `box_survival(N: int, L: int)` — 存活比例：box_total ÷ (2^N·(L−1))（起點在 L−1 個內點上平均）。
- `box_survival_long(N: float, L: int, kind: str='lattice')` — 存活比例的長時間直線（對數圖用）：'lattice' 格子最低模態 (2/L)cot²(π/2L)cos^N(π/L)/(L−1)；
- `decay_rate(L: float, kind: str='lattice')` — 長時間每步的衰減率（存活比例 ∝ e^{−rate·N}）：'lattice' −ln cos(π/L)、'cont' π²/(2L²)、'book' π²/L²。
- `odd_inv_squares(M: int)` — Σ_{m<M} 1/(2m+1)² ＋ Euler–Maclaurin 尾巴（→ π²/8，(4.22)）。
- `zeta2(M: int)` — Σ_{n≥1} 1/n²（部分和＋Euler–Maclaurin 尾巴；(4.21)）。
- `j0(x: float)` — (4.32) sin x / x。
- `j0_zero(n: int)` — j0 的第 n 個零點（二分法；應等於 nπ，(4.33)）。
- `sphere_weight(n: int, R: float)` — (4.34) 的係數 8R³/(πn²)：第 n 個 l = 0 模態對「走法總數」的權重。
- `sphere_weight_numeric(n: int, R: float)` — 同上，數值積分：ψ = c·j0(kr)，c² = 1/∫|j0|²d³r；權重 (∫ψ d³r)²（第二種算法；更正後的 (4.31)）。
- `sphere_weight_book(n: int, R: float)` — 書上 (4.31) 印的前因子 (8π/R)·(2/R)·(∫r²j0 dr)²（和 (4.34) 對不上）。
- `sphere_total_cont(N: float, R: float, book: bool=False, nmax: int=4000)` — (4.35) ÷ 6^N：(8R³/π) Σ e^{−π²n²N/(6R²)}/n²（更正）；book=True 用書上的 e^{−π²N/R²}（沒有 n²、沒有 1/6）。
- `sphere_rate(R: float, book: bool=False)` — 球的長時間衰減率（每步）：π²/(6R²)（更正後的 (4.37)）；書上 π²/R²。
- `ball_sites(R: float)` — sc 格點裡 x²+y²+z² < R² 的點（踩到球外就被吸收）。
- `ball_totals(R: float, Nmax: int)` — 格子球：起點在球內任一點、N 步都沒出球的走法數 ÷ 6^N，N = 0..Nmax。
- `ball_rate(R: float, Nmax: int=120)` — 格子球的長時間衰減率：−½ ln(T_{N+2}/T_N)（偶數 N，夠大；第二種算法，核對 1/6）。
- `cube_total(N: int, L: int)` — Ex. 4.5：立方體 0 < x, y, z < L（邊長 L、牆在 0 與 L）裡 N 步都沒碰牆的走法總數（逐步遞迴，精確）。
- `cube_total_modes(N: int, L: int)` — 同上 ÷ 6^N，用格子模態：Σ_{三個奇數 m} Π (2/L)cot²(m_iπ/2L) · ((Σ cos(m_iπ/L))/3)^N。
- `cube_total_cont(N: float, L: float, mmax: int=60)` — Ex. 4.5 的連續極限 ÷ 6^N：(8L/π²)³ Σ_{奇 n1,n2,n3} e^{−π²N(n1²+n2²+n3²)/(6L²)}/(n1n2n3)²。
- `cube_long(N: float, L: float)` — 立方體長時間 ÷ 6^N：(512L³/π⁶) e^{−π²N/(2L²)}。
- `same_volume_rate_ratio()` — 同體積的立方體與球：衰減率之比 (π²/2L²)/(π²/6R²) = 3(R/L)² = 3(3/(4π))^{2/3}。
- `laplacian_coeff(a: float=0.001)` — (4.40) 的係數：[(1/6)Σ_i f(r+a ê_i) − f(r)] / (a²∇²f) → 1/6（書上印 1/3）。測試函數 f = e^{−(x²+2y²+3z²)/2}。
- `tent_lattice(L: int, x0: int)` — 一維格子的穩定態：每步在 x0 放 1 個走路人，(4.39) 的一維版 w(x) = ½(w(x−1)+w(x+1)) + δ_{x,x0}，牆在 0、L。w(0..L)。
- `tent_formula(L: int, x0: int, x, s0=1, beta=Fraction(1, 2))` — (4.49)：(s0/β)·x(L−x0)/L（x < x0）、(s0/β)·x0(L−x)/L（x > x0）。格子上 β = a²/(2τ) = ½。
- `absorb_left(L: int, x0: int)` — Ex. 4.6：從 x0 出發、先被左牆 0 吸收的機率 (L − x0)/L（賭徒破產）。
- `absorb_left_solve(L: int, x0: int)` — 同上，解調和方程 h(x) = ½(h(x−1)+h(x+1))，h(0) = 1、h(L) = 0（第二種算法；Fraction）。
- `absorb_left_dp(L: int, x0: int, tol: float=1e-14)` — 同上，逐步遞迴到幾乎全被吸收（第三種算法；float）。
- `wall_fluxes(L: int, x0: int)` — 穩定態流進左、右牆的走路人數（每步）：w(1)/2、w(L−1)/2（Fraction；和 = 1）。
- `heat_kernel(r: float, t: float, D: float=1.0)` — 三維高斯 (4πDt)^{−3/2} e^{−r²/(4Dt)}（(4.54) 的正規化版本）。
- `time_integral(r: float, D: float=1.0)` — ∫_0^∞ heat_kernel dt（數值；t = e^v）→ 1/(4πDr)（(4.52)(4.55)(4.56)）。
- `K_const(a: float)` — (4.56) 的 K = ∫_0^∞ T^{−3/2} e^{−1/(aT)} dT = √(πa)。
- `K_const_numeric(a: float)` — 
- `ine_list(nmax: int, y: float)` — [e^{−y} I_n(y)]，n = 0..nmax（y < 1 用級數；否則 Miller 逆向遞迴，用 e^{y} = I_0 + 2Σ I_k 正規化）。
- `sc_green(n1: int, n2: int, n3: int)` — sc 格子的 Green 函數 G(z_c; r)：從原點出發、到訪 r = (n1, n2, n3) 的平均次數
- `sc_green_asym(r: float)` — 遠處 G ≈ 3/(2πr)：更正後的 (4.41) −∇²w = (6/a²)S 的解 6/(4πr)·(a = 1)。
- `sc_green_book(r: float)` — 照書上 (4.41) 的 3/a² 會得到 3/(4πr)（差 2 倍）。
- `sc_laplace_residual(n1: int, n2: int, n3: int)` — 離散 (4.39)：G(r) − (1/6)Σ G(鄰居) − δ_{r,0}（應為 0；核對 sc_green 的數值）。
- `conc_wall(x: float, y: float, z: float, l: float)` — (4.57)(4.58)：點源在 (l,0,0)、鏡像源在 (−l,0,0)、吸收牆 x = 0；β = s0 = 1：(1/4π)(1/|r−r_s| − 1/|r−r_i|)。
- `conc_dipole(x: float, y: float, z: float, l: float)` — (4.58) 末行（r ≫ l）：(1/4π)·2lx/r³。
- `plane_flux(X: float, l: float)` — 穿過平面 x = X（往 +x）的總流量 ∫∫ j_x dy dz（j = −∇c，精確的 (4.58) 第一行；數值積分，ρ = X·tan u）。
- `plane_flux_exact(X: float, l: float)` — 同上的閉式：每個點源穿過無限大平面的流量是 ±½（立體角 2π）。X > l：½ − ½ = 0；0 < X < l：−½ − ½ = −1（全部流向牆）。
- `dipole_flux_integral(x: float)` — (4.60)：∫_0^∞ [3x²/(x²+R²)^{5/2} − 1/(x²+R²)^{3/2}] R dR（數值；應為 0）。
- `conc_sphere(r: float, r0: float, cinf: float=1.0)` — (4.61) c∞(1 − r0/r)。
- `sphere_flux(r0: float, D: float=1.0, cinf: float=1.0)` — (4.63)/(4.65)(4.66)：吸收球每單位時間吃掉 4πDr0c∞（D 就是書上的 β）。
- `sphere_flux_numeric(r: float, r0: float, D: float=1.0, cinf: float=1.0)` — 在半徑 r 的球面上算 4πr²·D·dc/dr（數值微分；任何 r 都相同）。
- `kelvin_image(R_: float, r0: float)` — 吸收球外的點源（距球心 R）：鏡像源在 ρ = r0²/R、強度 −r0/R（Fig. 4.6）。
- `hit_prob(R_: float, r0: float)` — 從距球心 R 處出發，終究撞到半徑 r0 吸收球的機率 r0/R（逃脫 (R − r0)/R）。
- `sphere_surface_max(R_: float, r0: float)` — 點源＋鏡像源在球面上的濃度最大絕對值（應為 0：球面是吸收邊界）。
- `hit_prob_numeric(R_: float, r0: float)` — 點源（每單位時間 1 個）＋鏡像源：球面上的流入量 ∮ ∂c/∂r dS（數值積分；第二種算法，→ r0/R）。
- `cell_rmax(D: float, cinf: float, q: float)` — 營養流入 4πDr0c∞ ＝ 需求 (4/3)πr0³q → 最大半徑 √(3Dc∞/q)。
- `capacitance_sphere(r0: float)` — (4.66) C = r0（高斯單位）。
- `disk_capacitance(a: float)` — （書外）半徑 a 的導體圓盤 C = 2a/π；嵌在不吸收表面上的受器只有一面，等效 a/π，流量 4Dac∞。
- `disk_capacitance_numeric(a: float)` — 同上，數值：面電荷 σ ∝ 1/√(a²−ρ²)，總電荷 Q 與圓心電位 φ(0) 都用數值積分（ρ = a sin u），C = Q/φ(0)。
- `disk_potential_offcenter(a: float, p: float)` — 同一個電荷分布（σ = 1/√(a²−ρ²)）在圓盤上距圓心 p 處的電位（數值；等於 φ(0) = π²，導體是等電位）。
- `berg_purcell(n: float, a: float, r0: float)` — J/J_max = na/(na + πr0)：n 個半徑 a 的圓盤受器（Berg & Purcell 1977）。
- `berg_purcell_series(n: float, a: float, r0: float)` — 同上，用「電容串聯」：1/C_eff = 1/r0 + 1/(n·a/π)，J/J_max = C_eff/r0（第二種算法）。
- `bp_half(a: float, r0: float)` — 半飽和的受器數 n = πr0/a。
- `bp_area_fraction(n: float, a: float, r0: float)` — 受器占球面的比例 n·πa²/(4πr0²) = na²/(4r0²)。
- `half_counts(N: int, x: int, y: int)` — 書上 §4.3 的約定：格點 n ≥ 0，牆緊鄰在 −1 的左邊；從 x 走 N 步到 y、從沒踩到 −1 的走法數（逐步遞迴）。
- `half_counts_image(N: int, x: int, y: int)` — 同上，鏡像點在 −x − 2（對 −1 鏡射）：C(N; y−x) − C(N; y+x+2)。
- `half_gf(z: float, x: int, y: int)` — Ex. 4.1 的格子版（同一約定）：G'(z; x, y) = (λ^{|x−y|} − λ^{x+y+2})/√(1−4z²)，λ = (1−√(1−4z²))/(2z)。
- `half_gf_residual(z: float, x: int, ymax: int=6)` — (4:S-1)–(4:S-3)：在 y ≥ 0 都滿足 G(x,y) = z(G(x,y−1)+G(x,y+1)) + δ_{x,y}，且 G(x,−1) = 0（回傳最大誤差）。
- `catalan(k: int)` — 
- `half_visits_critical(x: int, y: int)` — z → z_c = ½：從 x 出發、到訪 y 的平均次數 2(min(x, y) + 1)（帳篷在 L → ∞ 的極限）。
- `check()` — 

## `combkit.walk_cond`

walk_cond.py —— 《Elements of the Random Walk》第 10 章 Walks and the O(n) model: mean field theory and spin waves 的數學

常數：`F`、`TORI`、`STEP_DIRS`、`STEP_CHARS`、`SAW_INF`、`GRAY_Q4`、`TORUS_SEED`、`TORUS_START`、`TORUS_EXAMPLE_N`、`TORUS_EXAMPLE_WRAPS`、`MF_U`、`MF_N0`、`MF_R`、`MF_N0S`、`MF_NS`、`MF_NS_NEG`、`WALL_X`、`SW_M`、`SEED_SW`、`GAUSS_M`、`GAUSS_H`、`HOLE_L`、`HOLE_R`、`HOLE_D`、`SC_EWALD`、`CONT_RHO`、`VIS_N`、`VIS_R`、`RHO_Z`、`MEAN_N_Z`、`MF_RHO_N0`、`MF_RHO_R`、`RR_CASES`

- `torus_neighbors(L)` — 
- `torus_counts(L)` — 8 重對稱：第一步向東、第一次轉彎向北的 ×8；一路直走的 ×4。
- `torus_counts_plain(L, Nmax=None, sym=1)` — sym=1：完全不用對稱；sym=4：只固定第一步向東（×4）。Nmax 限制深度。
- `torus_counts_brute(L, N)` — 
- `torus_vs_inf(L)` — 
- `torus_first_diff(L)` — 
- `torus_walks(L, N)` — 
- `torus_steps(path, L)` — 
- `torus_lift(path, L, start=(0, 0))` — 
- `torus_segments(path, L, start=(0, 0))` — 
- `torus_wraps(path, L, start=(0, 0))` — 
- `torus_example(name, L=4, start=TORUS_START, seed=TORUS_SEED)` — 
- `torus_moments(L, z)` — 
- `torus_mean_N(L, z)` — 
- `book101_mean_N(L, z, n_comp=0)` — 照書上 (10.1) 的 z^n：G = Σ_N c_N z^n 與 N 無關，⟨N⟩ = z ∂ln G/∂z 只剩 n = 0。
- `mf_f(lam, r, u=MF_U, N0=MF_N0)` — 
- `mf_I(r, u=MF_U, N0=MF_N0, method='erfcx')` — 
- `mf_lnI(r, u=MF_U, N0=MF_N0)` — 
- `gamma_mf(r, u=MF_U, N0=MF_N0, book=False)` — 
- `ln_gamma_mf(r, u=MF_U, N0=MF_N0)` — 
- `gamma_n(n, r=MF_R, u=MF_U, N0=MF_N0, method='quad')` — 有限 n：N0⟨M1²⟩ = (N0/n)·∫λ^{n/2} f / ∫λ^{n/2−1} f（徑向座標 d^nM ∝ λ^{n/2−1}dλ）。
- `gamma_n1_direct(r=MF_R, u=MF_U, N0=MF_N0)` — 
- `gamma_gauss_check(r=MF_R, N0=MF_N0, ns=(3, 1, 0.01))` — u → 0（高斯）：任何 n 都是 1/r；更正的 (10.15) 給 1/r、書上 1/(2r)。
- `s_integral(lam, eps=0.1)` — 
- `mf_asym_I(r, u=MF_U, N0=MF_N0)` — 
- `gamma_mf_asym(r, u=MF_U, N0=MF_N0)` — 
- `ln_gamma_mf_asym(r, u=MF_U, N0=MF_N0)` — 
- `fn3_check(r=MF_R, u=MF_U, N0=MF_N0)` — 
- `wall_M(x, r=WALL_R, u=WALL_U)` — 
- `wall_bvp(r=WALL_R, u=WALL_U, X=20.0, n=2001)` — 
- `wall_errors(r=WALL_R, u=WALL_U)` — 
- `ex101_grad_check(seed=SEED_SW, n=2, sites=12, r=-0.8, u=1.3, hstep=0.5)` — 一維週期鏈、n 分量：ℋ = Σ hstep[(r/2)|M|² + (u/4)|M|⁴ − h·M + ½|ΔM/hstep|²]；數值梯度 vs (10.17) 的離散版。
- `ex102_coeffs(seed=SEED_SW, n=3, r=-0.7, u=1.3)` — g(ε) = (r/2)|M+εσ|² + (u/4)|M+εσ|⁴ − h·(M+εσ) 的 ε^k 係數：數值（5 點內插）vs 公式；回傳最大誤差。
- `sw_lambdas(r=SW_R, q2=SW_Q2, u=SW_U, M2=None)` — 
- `sw_X(M=SW_M, r=SW_R, q2=SW_Q2, u=SW_U)` — 
- `sw_Xinv(M=SW_M, r=SW_R, q2=SW_Q2, u=SW_U)` — 
- `sw_check(ns=(1, 2, 3, 4, 5, 6), seed=SEED_SW, trials=4)` — 
- `mode_table(ns=(3, 2, 1, 0))` — 
- `gauss_n2_check(npts=1201, lim=12.0)` — 
- `D_box(L=HOLE_L, r=HOLE_R, u=1.0)` — 
- `D_axis(L=HOLE_L, r=HOLE_R)` — 
- `D_direct(L, r, delta, u=1.0)` — 
- `D_book46(L=HOLE_L, r=HOLE_R)` — 照書上 (10.46)：Σψψ'/(q²(2|r|+q²))（少了分子的 2|r|）＝ D/(2|r|)。
- `D_inf(delta, r=HOLE_R, u=1.0)` — 
- `box_shift(L, r, delta)` — 
- `coulomb_fit(L=HOLE_L, r=HOLE_R, dmin=4, dmax=16)` — 
- `D_cont(rho, r=HOLE_R)` — 
- `D_cont_quad(rho, r=HOLE_R)` — 
- `case1_check(L=32, r=1.0, dmax=4)` — 
- `sc_prob(n, R)` — 
- `sc_visits(N, R)` — 
- `sc_visits_dp(N, R)` — 
- `erfc_form(R, N)` — 
- `erfc_form_book(R, N)` — 
- `eq1080_check(k2=0.7, zc=6.0, N=12)` — (10.80)：z^N 的係數（級數摺積）＝ z_c^{−N} ×（書上的左邊）＝ z_c^{−N}·[1 − (1+k²/z_c)^{−(N+1)}]/(z_c k²)。
- `eq1083_check(R=3.0, a=2.0)` — (10.82)→(10.83) 的 Fourier 積分：∫d³k/(2π)³ e^{ik·R}(1 − e^{−ak²})/k² ＝ erfc(R/(2√a))/(4πR)（QAWF 數值 vs 閉式）。
- `book45_argmax(r=HOLE_R, u=1.0)` — (10.45) 的指數：照書上 −(rλ + uλ²) 的極大在 |r|/(2u)；(10.15) 的 f 在 |r|/u（數值極值）。
- `green_sc(R)` — 
- `green_sc_series(R, N=400)` — 第二種算法：精確的 Σ_{n≤N} P_n(R) ＋ 連續極限的尾巴 (3/(2πR)) erf(R√(3/(2(N+1))))。
- `torus_rho(L, z)` — 
- `torus_rho_fd(L, z, h=1e-05)` — ⟨N⟩ = −∂lnΓ/∂r，e^{−r} = z/z_c（中央差分）。
- `torus_var_fd(L, z, h=0.0001)` — 
- `torus_peak(L)` — 
- `torus_curve(L, n=60, zmax=1.2)` — 
- `mf_rho(r, u=MF_U, N0=MF_RHO_N0)` — 
- `mf_chi(r, u=MF_U, N0=MF_RHO_N0)` — 
- `mf_rho_fd(r, u=MF_U, N0=MF_RHO_N0, h=None)` — 
- `mf_chi_fd(r, u=MF_U, N0=MF_RHO_N0, h=None)` — 
- `mf_rho_quad(r, u=MF_U, N0=10 ** 4)` — 
- `book55_rho(r, u=MF_U)` — 書上 (10.55)：Ω = (1/V){ln V − V r²/(4u) + …} → ρ = −∂Ω/∂r = r/(2u)（r < 0 時是負的）。
- `lnG53(r, N0, u=MF_U, sign=+1)` — (10.53) r < 0（照書上的 N0/4 正規化，Γ 多乘一個 N0）：ln[(1/2)(π/u)^{1/2} N0^{3/2} e^{±N0r²/4u}]。
- `rr_exponents(nu, d)` — 
- `table101(alpha, nu)` — 
- `rr_laws(e, d)` — 
- `dc_mean_field()` — 
- `summary()` — 
- `check()` — 

## `combkit.walk_crit`

walk_crit.py —— 《Elements of the Random Walk》第 9 章 Scaling of walks and critical phenomena 的數學

常數：`F`、`NU_KNOWN`、`TC_2D`、`MS_X`、`BETA_CITED`、`PAIR`、`RING4`、`TORUS4`、`TINY_K`、`FIELD_H`、`MC_SEED`、`MC_SWEEPS`、`MC_REPLICAS`、`MC_M0`、`MC_T`、`MC_L`、`MC_TRACES`、`SNAP_T`、`K1D`、`FISHER1D_K`、`EXPONENTS`、`BOOK_3D`、`CITED_3D`、`WIDOM_CASE`、`LANDAU_R`、`LANDAU_H`、`EX94_U`、`ZMF_V`、`ZMF_H`、`ZC_2D`、`SAW_Z`、`S_GAMMA`、`S_GAMMA2NU`、`COEF_N`、`SAW_N`、`WEAK_N`、`WEAK_U`、`GAUSS_T`、`GAUSS_H`、`TROTTER_EPS`、`TROTTER_M`、`LADDER`、`LADDER_K`、`LADDER_PAIRS`

- `flory_p(d)` — (9.2) p = 6/(d+2)（Fraction）＝ 2ν_Flory。
- `flory_gap_3d(nu3=0.588)` — 三維 Flory 3/5 比已知 ν 大多少（相對）。
- `mf_m(T)` — 平均場 m = tanh(m/T) 的正根（T_c = 1；brentq）。
- `mf_m_newton(T, n=500)` — 同上，用牛頓法從 m = 1 往下（g(m) = m − tanh(m/T) 在 m > 0 是凸的，單調收斂；第二種算法）。
- `mf_m_asym(T)` — T → 1⁻ 的領頭項 √(3(1−T))（β = 1/2）。
- `onsager_M(T)` — 二維 Ising（正方格子，J = 1）的自發磁化 (1 − sinh^{−4}(2/T))^{1/8}。
- `onsager_lowT(T)` — 低溫團簇展開：翻一個自旋（能量 8，磁化少 2）＋相鄰兩個（能量 12、兩種方向、少 4）：1 − 2u² − 8u³，u = e^{−4/T}。
- `beta_slope(fn, Tc, eps=0.0001)` — ln M 對 ln(1 − T/T_c) 的局部斜率（取 1 − T/T_c = eps 與 eps/10）。
- `ms_table()` — 
- `ring_edges(Ns)` — 
- `torus_edges(L)` — L×L 週期正方格子（L ≥ 3）的 2L² 條鍵；格點 (x, y) 的編號 x + L·y。
- `on_energy(S, edges)` — (9.4) H = −Σ_{鍵} S_i·S_j（J = 1）。
- `rotation_spins(n=ROT_N, Ns=ROT_NS, seed=902, spread=ROT_SPREAD)` — Ns 個單位 n 維向量排成一圈（大致朝 x 方向、亂度 spread）與「全部一起轉」之後的樣子 → (S, S')，各 (Ns, n)。
- `rotation_check(n=ROT_N, Ns=ROT_NS, seed=902, spread=ROT_SPREAD)` — (E, E', M⃗, M⃗')：能量 (9.4) 不變、|M⃗| 不變、方向轉了。
- `angle_deg(v)` — 二維向量的方向角（度，−180..180）。
- `ising_exact(graph, K, H=0.0)` — 逐一列舉 2^Ns 個 ±1 組態（Ns ≤ 16），權重 e^{K Σ_鍵 s s + H Σ s}。
- `ring4_closed(K)` — 4 個自旋一圈（H = 0）的閉式：(Z, ⟨|m|⟩, P(M=4), P(M=2), P(M=0))。
- `pair_closed(K)` — 2 個自旋（一條鍵）：(Z, ⟨|m|⟩) = (2e^K + 2e^{−K}, e^K/(e^K + e^{−K}))。
- `mc_run(L, T, sweeps=MC_SWEEPS, replicas=MC_REPLICAS, seed=MC_SEED)` — 二維 Ising L×L 週期、棋盤式 Metropolis（全部從 +1 開始）→ (m, e)：各 (sweeps, replicas)，每個 sweep 之後記一次。
- `mc_snapshot(L=8, T=MC_T, sweeps=400, seed=MC_SEED)` — 單一系統（從 +1 開始）跑 sweeps 次之後的自旋 (L, L)（±1；畫箭頭用）。
- `flip_count(ms, m0=MC_M0)` — m 在「> m0」與「< −m0」兩個狀態之間換邊的次數（從 + 開始；0 附近的抖動不算）。
- `flip_tau(L, T=MC_T, sweeps=MC_SWEEPS, replicas=MC_REPLICAS)` — 平均翻轉等待時間（sweep）＝ 總 sweep 數 × 複本數 ÷ 翻轉次數。
- `trace_replica(L, T)` — 畫哪一個複本：T < T_c 時取第一個「翻轉次數 = 50 個複本的平均（四捨五入）」的複本（典型的一條），否則第 0 個。
- `mc_trace(L, T, n_pts=400)` — 畫圖用：trace_replica 那一條的 m(t)，均勻取 n_pts 點 → (sweep 陣列, m 陣列)。
- `trace_flips(L, T)` — 
- `mc_check_L4(T=MC_T)` — L = 4 的 Metropolis 平均 vs 精確列舉 → ((⟨|m|⟩_MC, 精確), (e_MC, 精確))。
- `ising1d_G(r, K)` — 
- `ising1d_xi(K)` — 
- `ising1d_chi(K)` — 
- `ising1d_chi_partial(K, R)` — 
- `ising1d_transfer_G(r, K)` — 轉移矩陣 [[e^K, e^{−K}], [e^{−K}, e^K]] 的兩個本徵值之比的 r 次方（無窮長鏈）。
- `ising1d_ring_closed(N, K)` — N 個自旋的環：G(r) = (t^r + t^{N−r})/(1 + t^N)（r = 0..N−1）、χ_N = Σ_r G(r)。
- `ising1d_xi2_moment(K)` — Σ r² G / Σ G（r 從 −∞ 到 ∞）＝ 2t/(1−t)²（Table 9.1 的比值是長度的平方）。
- `ising1d_xi2_moment_sum(K, R=4000)` — 
- `laws(e)` — 
- `from_nu_eta(nu, eta, d)` — 由 ν、η、d 推出 (α, β, γ, δ)：hyperscaling、Fisher、β = ν(d−2+η)/2。
- `book3d_alpha()` — 
- `widom_exponents(alpha, Delta)` — Ex. 9.1–9.2：β = 2 − α − Δ、γ = Δ − β、δ = Δ/β。
- `widom_G(r, h, alpha, Delta)` — 一個具體的 Widom 形式 (9.23)：G = |r|^{2−α} f(h/|r|^Δ)（r < 0；f 見 _wf，q = (2−α)/Δ）；r = 0 取極限 −(|h|^q + h|h|^{q−1})。
- `widom_slopes(alpha, Delta, method='exact')` — 數值斜率 (β, γ, 1/δ)：M(r, 0) ∝ |r|^β、χ(r, 0) ∝ |r|^{−γ}、M(0, h) ∝ h^{1/δ}（r = −1e-2 與 −1e-4、h = 1e-3 與 1e-5）。
- `fourier_check(N=8, seed=904)` — (9.30)–(9.32) 的離散版（一維 N 點環，S(q) = (1/N)Σ S_i e^{iqr_i}）→ 三個相對誤差。
- `landau_F(m, r, u=1.0, h=0.0)` — 
- `landau_min(r, u=1.0, h=0.0)` — ∂F/∂m = u m³ + r m − h = 0 的實根中 F 最小者（h = 0、r < 0 時取正根）。
- `landau_closed(r, u=1.0)` — h = 0 的閉式：(m*, F_min) = (√(−r/u), −r²/4u)（r < 0），r ≥ 0 為 (0, 0)。
- `ex94_m(r, h)` — 
- `ex94_closed(r)` — (9.35) h = 0、r < 0：(m*, F_min) = (√(−3r), −3r²/4)。
- `ex94_scaling(x, r)` — m/|r|^{1/2}，h = x|r|^{3/2}（r < 0）→ ℳ(x)（與 r 無關）。
- `ex94_slopes()` — 數值：(β, Δ, α, 1/δ, γ)。β 由 m*(r)、Δ 由「固定 m/|r|^β 所需的 h」、α 由 F_min(r)、1/δ 由 m(0, h)、γ 由 χ = 1/F''(m*)。
- `zmf_avg(V, r=-1.0, u=1.0, h=0.0, what='m', method='quad')` — (9.34) n = 1 的有限體積積分 ∫ dm e^{−V F(m)}：what = 'm'（⟨m⟩）、'absm'（⟨|m|⟩）、'f'（−ln Z/V）。
- `saw_gf(z, Nmax=14)` — (9.37)：χ = Σ c_N z^N、分子 Σ c_N⟨R²⟩_N z^N、ξ² = 分子/χ（walk_saw 的精確列舉）。
- `saw_gf_realspace(z, Nmax=8)` — 由終點：G(z; r⃗) = Σ_N (N 步到 r⃗ 的自迴避條數) z^N，再算 Σ r²G / Σ G（全部 4^N 條暴力篩選；第二種算法）。
- `saw_table(Nmax=6)` — 
- `coef_log(N, s, zc)` — (z_c − z)^{−s} 的 z^N 係數的對數：ln[z_c^{−N−s} Γ(N+s)/(Γ(s) N!)]（lgamma）。
- `coef_log_rec(N, s, zc)` — 同上，由遞迴 a_0 = z_c^{−s}、a_k = a_{k−1}(k−1+s)/(k z_c) 的對數和（第二種算法）。
- `coef_exact(N, s, zc)` — 
- `coef_rec(N, s, zc)` — 
- `coef_asym(N, s, zc)` — 
- `coef_ratio(N, s)` — exact ÷ N^{s−1} z_c^{−N−s}/Γ(s)（z_c 消掉）。
- `coef_ratio_book(N, s, delta)` — 照書上 N^{δ−1}（而不是 N^{s−1}）的比值。
- `r2_prefactor(gamma, nu, zc, K1=1.0, book=False)` — (9.41) 的前因子 K₁Γ(γ)/(z_c^{2ν} Γ(γ+2ν))；book=True 照書上 z_c^γ。
- `r2_ratio(N, gamma, nu, zc, rec=False)` — [z^N]K₁K₂(z_c−z)^{−γ−2ν} ÷ [z^N]K₂(z_c−z)^{−γ} ÷ (前因子·N^{2ν})（K₁ = K₂ = 1；→ 1）；rec=True 用遞迴的係數。
- `zc_pow()` — (z_c^{−2ν}, z_c^{−γ}) = (μ^{3/2}, μ^{43/32})（二維）。
- `saw_amp(N)` — 
- `hs_lhs(s, u=1.0)` — 
- `hs_rhs(s, u=1.0)` — (9.46) 右邊：(1/√(πu))∫_{−∞}^{∞} dφ e^{−φ²/u} cos(φs)（虛部是奇函數，積分為 0）。
- `hs_book72(u=1.0, method='quad')` — 照書上 (9.72) 從 0 積到 ∞：(1/√(πu))∫_0^∞ e^{−φ²/u − iφ} dφ（複數）。閉式 ½e^{−u/4}(1 − i·erfi(√u/2))。
- `visit_weight(k, u=1.0)` — 一個格點被走過 k 次：(1/√(πu))∫dφ e^{−φ²/u − ikφ} = e^{−uk²/4}（(9.72) k = 1、(9.73) k = 2）。
- `intersection_factor(u=1.0)` — 
- `pair_sum(k)` — 同一格點走過 k 次時 Σ_{l,l'} δ（有序對、含 l = l'）＝ k²＝ k ＋ 2·C(k,2)。
- `weak_saw_count(N, u, sign=-1)` — Σ_w e^{sign·(u/2)·n_w}（N 步、全部 4^N 條；u = inf 時只剩 n_w = 0 的自迴避）。
- `heat_kernel(x, eps)` — e^{½∇²ε} 的核（一維）：(2πε)^{−1/2} e^{−x²/2ε}。
- `heat_kernel_fourier(x, eps)` — (1/2π)∫dk e^{ikx} e^{−εk²/2}（數值；第二種算法）。
- `heat_kernel_book(x, eps, d=1)` — 照書上 (9.70)：(2π)^{−d/2} exp(−(d/2)x²/ε)（一個方向的因子）。
- `kernel_mass_book(eps, d=1)` — 書上的核對 d 維空間積分 = (ε/d)^{d/2}（數值，d = 1 的一維積分取 d 次方）。
- `gauss_I(T=GAUSS_T, h=GAUSS_H)` — (9:S-11)：(2π)^{N/2}/√det T · exp(½ hᵀT^{−1}h)。
- `gauss_I_quad(T=GAUSS_T, h=GAUSS_H)` — 同上，2×2 直接數值積分（dblquad）。
- `gauss_moment(T, a, b, n=20)` — ρ ∝ exp(−½xᵀTx) 下的 ⟨x₁^a x₂^b⟩（2×2；Gauss–Hermite 張量求積，對多項式精確）。
- `wick4(C, r, s, t, u)` — (9:S-20)：⟨x_r x_s x_t x_u⟩ = C_rs C_tu + C_rt C_su + C_ru C_ts。
- `wick4_book(C, r, s, t, u)` — 照書上 (9:S-19) 第三行括號裡的 ⟨x_r x_s⟩（應為 ⟨x_t x_u⟩）再微分得到的結果。
- `wick_count(m)` — 
- `wick_count_enum(m)` — 
- `trotter_mats(N=6, seed=906)` — X = −iφ（對角，φ ~ N(0,1) 固定種子）、Y = ½∇²（N 點環的離散 Laplacian）。
- `trotter_err(eps, order=1)` — ‖e^{(X+Y)ε} − e^{Xε}e^{Yε}‖（order=1）；order=2 再乘 e^{−(ε²/2)[X,Y]}（Ex. 9.6 的下一項）。
- `trotter_lead(eps)` — 一階誤差的領頭項 (ε²/2)‖[X,Y]‖。
- `trotter_global(M, t=1.0)` — 
- `ex95_check(r=1.0)` — Ex. 9.5：V = (r/2 − ½∇²)I + iφ（6 點環、φ 同 trotter_mats）；V⁻¹ = ∫_0^∞ e^{−tV}dt（quad_vec）
- `ladder_G(n, K, pair)` — 閉式（|S|² = n、權重 Π(1 + K S·S)）：只有「端點度數 1、其他度數 0 或 2」的鍵集合存活；每個閉圈乘 n。
- `ladder_G_ising(K, pair)` — n = 1：64 個 ±1 組態暴力加總，權重 Π(1 + K s_i s_j)。
- `ladder_G_xy(K, pair, m=8)` — n = 2：S = √2(cos θ, sin θ)，6 個角度各取 m 個等距點（被積函數是低次三角多項式，梯形法精確）。
- `ladder_saws(a, b)` — 梯子上從 a 到 b 的自迴避路徑長度（DFS）。
- `ladder_saw_sum(K, pair)` — 
- `summary()` — 
- `check()` — 

## `combkit.walk_eps`

walk_eps.py —— 《Elements of the Random Walk》第 12 章 More on the renormalization group 的數學（全書最後一章）

常數：`F`、`TOY`、`SHELL_CASES`、`LC`、`LC_MAT`、`NU_3D_CITED`、`NU_EPS`、`FLOW_U0`、`U_CASES`、`BOOK_RECAP`、`BOOK_SCENES`、`RECAP_EXPECT`

- `toy_H(s, Sp, r, u)` — 一個「內區」變數 s、一個「殼層」變數 S'：ℋ = r s²/2 + u s⁴ + (r+1) S'²/2 + 6u s² S'²（只留 (12.11) 的 6u 頂點）。
- `toy_A(r)` — (12.7) 的 A：S' 的 Gaussian 正規化 √(2π/(r+1))。
- `toy_Heff(s, r, u)` — 閉式：ℋ'(s) = −ln[∫dS' e^{−ℋ} / A] = r s²/2 + u s⁴ + ½ ln(1 + 12u s²/(r+1))。
- `toy_Heff_numeric(s, r, u)` — 同上，對 S' 做數值積分（第二種算法）。
- `toy_coeffs(r, u)` — ℋ' 的 Taylor 係數 → (r', u', c6)：r' = r + 12u/(1+r)、u' = u − 36u²/(1+r)²、s⁶ 的新項 c6 = (12u/(1+r))³/6。
- `toy_coeffs_numeric(r, u)` — 數值 ℋ'(s) 對 s² 做多項式擬合（第二種算法）→ (r', u')。
- `toy_Z(r, u)` — (12.8)(12.9)：Z = ∫∫ e^{−ℋ} ds dS'（二重積分）與 A·∫ e^{−ℋ'(s)} ds（先積 S'）。
- `shell_fraction(d, dl)` — 殼層 e^{−ℓ−δℓ} < |q| < e^{−ℓ} 佔布里淵區的比例 1 − e^{−dδℓ}（(12.20)–(12.25)）。
- `shell_linear(d, dl)` — 
- `ball_count(d, R)` — |x| ≤ R 的整數格點數（d = 2、3；d = 3 對 z 用解析計數）。
- `shell_fraction_lattice(d, dl, L=200)` — 布里淵區裡「正好 N 個 q」：用半徑 L 的球內格點數數殼層的比例（第二種算法）。
- `shell_dq(ell, dl)` — (12.24)：dq = e^{−ℓ} − e^{−ℓ−δℓ} ≈ e^{−ℓ}δℓ → (精確, 一階)。
- `rescale_factors(d, dl, book=False)` — (12.32)–(12.34) 之後各項的乘數：q² 項、r 項、u 項。正確：q → q'e^{−δℓ}、S → S''e^{δℓ}、N → N'e^{dδℓ}
- `matchings(n)` — n 個四次片段（每個 2 條殼層腳：片段 i 的腳是 2i、2i+1）的全部配對：(2n−1)!! 個。
- `ring_sizes(m, n)` — 一個配對把片段連成的環（每個片段剛好 2 條腳 ⇒ 每個連通塊都是環）的大小，由大到小；
- `type_census(n)` — {環的大小 tuple: 配對數}（腳的層次）；總和 = (2n−1)!!。
- `diagram_census(n)` — {環的大小: 不同的「片段圖」數}（只看哪些片段相連、不管用哪一條腳；第二種數法）。
- `sizes_of(n1, n2=0, n3=0)` — 
- `pair_factor(n1, n2)` — (12:S-4)：n!/(n1! n2! (2!)^{n2})，n = n1 + 2n2。
- `ring_aut(k)` — f_k：k 個頂點排成環的「換頂點但配對不變」的方法數（列舉 k! 個排列）：1、2、6、8、10（k ≥ 3 是 2k）。
- `general_factor(counts)` — (12:S-6)：n!/Π(n_k! f_k^{n_k})，counts = {k: n_k}。
- `ring_weight(k)` — 一組 k 個（標了號的）片段連成一個環的腳配對數 w_k = 2^{k−1}(k−1)!：1、2、8、48。
- `egf_check(x, nmax=40)` — 連結叢集的最簡形：Σ_n (2n−1)!! xⁿ/n! = exp(Σ_k w_k x^k/k!) = (1−2x)^{−1/2} → (左, 右, 閉式)。
- `lc_moments(a, sig2, nmax)` — ⟨Xⁿ⟩ = (−a)ⁿ (2n−1)!! σ^{2n}（全部的圖，含不連通的），n = 0..nmax。
- `lc_cumulants_rings(a, sig2, kmax)` — 連通的圖＝k 個片段連成一個環：κ_k = (−aσ²)^k · 2^{k−1}(k−1)!（k = 1..kmax）。
- `lc_cumulants_moments(a, sig2, kmax)` — 由動差遞迴求累積量（第二種算法）：κ_n = m_n − Σ_{j<n} C(n−1, j−1) κ_j m_{n−j}。
- `lc_average(a, sig2, method='closed')` — ⟨e^X⟩：'closed' (1+2aσ²)^{−1/2}、'quad'（Gauss–Hermite）、'series'（Σ⟨Xⁿ⟩/n!）、'linked'（exp Σ κ_k/k!）。
- `lc_order(a, sig2, n)` — 第 n 階：(全部 ⟨Xⁿ⟩/n!, 連通 κ_n/n!, 不連通 = 差)。
- `lc_multi(A=LC_MAT, method='det')` — 三個殼層變數（σ² = 1）：⟨e^{−sᵀAs}⟩：'det' det(1+2A)^{−1/2}、'quad'（三維 Gauss–Hermite）、'linked'（環的和 Σ(−1)^k 2^{k−1} tr A^k / k）。
- `on_counts(n)` — O(n) 頂點 (S·S)(S·S) 的縮併（Ex. 12.5–12.7，numpy einsum）→ (2n+4, 4(n+8))：
- `AB(n=1)` — (12.66)(12.67) 的係數 A = 4(n+2)、B = 4(n+8)；Ising（n = 1）12、36。
- `beta(r, u, eps, n=1)` — dr/dℓ = 2r + A u/(1+r)、du/dℓ = εu − B u²/(1+r)²（本系列：書上 (12.37) 的 d 吸收進 u）。
- `beta_book(r, u, eps)` — 書上 (12.37)＋更正後的 (12.38)：d = 4 − ε 留在式子裡：2r + 12ud/(1+r)、εu − 36d u²/(1+r)²。
- `fixed_point(eps, n=1)` — Wilson–Fisher 不動點（精確，eps 可用 Fraction）：r* = −Aε/(2B+Aε)、u* = 4Bε/(2B+Aε)²；
- `fixed_point_numeric(eps, n=1)` — 
- `fixed_point_lowest(eps, n=1)` — 一階：r* ≈ −Aε/(2B)、u* ≈ ε/B（Ising −ε/6、ε/36；書上 u* → ε/36d）。
- `jacobian(r, u, eps, n=1)` — 
- `jacobian_fd(r, u, eps, n=1, h=1e-06)` — 
- `eigen(eps, n=1, which='wf')` — 線性化矩陣的特徵值（由大到小）：'gauss' 在 (0, 0) 是 (2, ε)；'wf' 在 Wilson–Fisher 點。
- `lambda_r_closed(eps, n=1)` — WF 點的相關特徵值（閉式）：矩陣 [[2 − Aε/B, ·], [·, −ε]]、非對角乘積 2Aε²/B。Ising ε = 1：(1+√22)/3 = 1.8968。
- `nu_first(eps, n=1)` — (12.59)(12.68)：ν = ½ + (n+2)ε/(4(n+8))（Fraction）。
- `nu_unexpanded(eps, n=1)` — 書上 (12.52)–(12.54) 不展開：ν = 1/(2 − Aε/B)（Ising 1/(2 − ε/3) = 3/(6−ε)，剛好是 Flory 3/(d+2)）。
- `nu_eigen(eps, n=1)` — 完整的一迴圈線性化：ν = 1/λ_r（數值特徵值）。
- `nu_from_flow(eps, n=1, deltas=(1e-05, 1e-07), Dist=0.05)` — (12.53)–(12.58) 的數值版：從 WF 點沿相關方向偏 δ 出發，量到 |r − r*| = Dist 的 ℓ*；ν = −Δℓ*/Δln δ（第三種算法）。
- `flow_line(r0, u0, eps=1.0, n=1, lmax=10.0, npts=40, box=(-0.45, 0.45, -0.001, 0.06))` — (r, u) 平面的一條流線（solve_ivp），跑出 box = (rmin, rmax, umin, umax) 就停 → (k, 2) 陣列（≤ npts 點）。
- `flow_end(r0, u0, eps=1.0, n=1, lmax=40.0)` — 往上（+1）或往下（−1）跑出去，或停在 WF 附近（0）。
- `r_critical(u0, eps=1.0, n=1)` — Ex. 12.4：臨界的 r(0)（二分法：再大就流向 +∞、再小就流向 −）。
- `r_critical_linear(u0, eps=1.0, n=1)` — Gaussian 點的 ε 特徵向量 (A, ε − 2)：r_c ≈ −A u0/(2 − ε)；ε → 0 為 −(A/2)u0（Ising −6u0，書上 r = t − 6ud）。
- `flow_starts(eps=1.0, n=1)` — 
- `u_analytic(ell, u0, eps, B=36)` — (12.41)（d 吸收進 u）：u(ℓ) = εu0 e^{εℓ}/(ε + B u0(e^{εℓ} − 1))（r = 0）。
- `u_numeric(ell, u0, eps, B=36)` — 
- `u_from_1240(u0, u1, eps, B=36)` — (12.40) 左邊的積分 ∫_{u0}^{u1} du/(εu − Bu²)（數值積分；應等於 ℓ）。
- `int_1248(u0, u1, eps, d, book=False)` — (12.47)–(12.49)：12d∫u dℓ = 12d∫ u du/(εu − 36du²) = −⅓ ln(ε − 36ud)|；book=True 照 (12.48) 少了分子的 u'。
- `int_1249(u0, u1, eps, d)` — 
- `t_closed(ell, t0, u0, eps, n=1)` — (12.51)（O(n)：指數 A/B）：t = t0 e^{2ℓ}[ε/(ε + B u0(e^{εℓ} − 1))]^{A/B}；大 ℓ ∝ e^{(2 − Aε/B)ℓ}。
- `t_numeric(ell, t0, u0, eps, n=1)` — (12.45) dt/dℓ = 2t − A u(ℓ) t 與 du/dℓ = εu − Bu² 一起數值積分。
- `flory_nu(d)` — 
- `flory_series(k)` — 3/(6 − ε) = ½ Σ (ε/6)^k 的第 k 個係數：½、1/12、1/72、1/432。
- `g_numeric(q0, g0, eps, B=36, eta=0.0)` — (12.86) q0 dg/dq0 = (−4 + d − 2η)g + B g²（−4 + d = −ε），從 q0 = 1 積到 q0。
- `pc_exponent(d, book=False)` — (12.89) 的冪次：u² r^{d−3}（＝ r·(u r^{−(4−d)/2})²）；book=True 照書上 u² r^{2d−3}。
- `pc_rhs(d, book=False)` — 右邊的組合：r^{+1}·(u r^{−(4−d)/2})² 的 r 冪次 d − 3；書上 r^{−1}·(…)² 是 d − 5。
- `watermelon(d, r)` — Fig. 12.27 的二迴圈（外動量 0）＝ ∫ dᵈx G(x)³（實空間；d < 3 時收斂）。
- `pc_numeric(d, r1=1.0, r2=4.0)` — 數值冪次：ln(I(r2)/I(r1))/ln(r2/r1) → d − 3（d = 2：−1、d = 2.5：−0.5）。
- `recap_numbers()` — 全書回顧表每一列的數字（字串）；check() 用兩種算法核對。
- `scene_range(ch)` — 第 ch 章的場景代號範圍，例如 'W101–W112'、'W1201–W1206'。
- `saw_size(N, nu)` — R ∝ N^ν 的例子：N = 1000 → √N = 31.62、N^{0.5625} = 48.70、N^{0.588} = 58.08。
- `summary()` — 場景卡與講義用的全部數字（字串已照卡上的小數位格式化）。
- `check()` — 

## `combkit.walk_gf`

walk_gf.py —— 《Elements of the Random Walk》第 2 章 Generating functions I 的數學（純 Python，不 import manim）

常數：`SQ_STEPS`、`SQ_NAMES`、`SC_STEPS`、`BCC_STEPS`、`FCC_STEPS`、`STEPS`、`UV_STEP`、`CHI_SMALL_Q`、`DIM`、`CELL_VOLUME`、`BIPARTITE`、`STEP_SCALE`、`CHI3_MEDIAN`、`SC_BASIS`、`BCC_BASIS`、`FCC_BASIS`、`MORTAL_MEAN`、`MORTAL_AT_ORIGIN`、`MORTAL_TRAIL_SEED`、`MORTAL_TRAIL_LEN`、`SQ_RATIOS`、`UNIV_TABLE`、`UNIV_SCALED10`、`CAUCHY_MEDIAN`、`MEDIANS`

- `exp_coeffs(n: int)` — 
- `half_binom_series(n: int)` — (1 − t)^{−1/2} = Σ C(2k,k)/4^k t^k。
- `legendre_gf(nmax: int)` — 把 (1 − t)^{−1/2}、t = 2xz − z² 展開，收集 z^n 的係數（x 的多項式）。
- `poly_mul(p, q)` — 
- `poly_eval(p, x)` — 
- `poly_int11(p)` — ∫_{−1}^{1} p(x) dx（精確）。
- `legendre_bonnet(nmax: int)` — 
- `ex22_integral(z: float)` — ∫_{−1}^{1} dx/sqrt(1 − 2xz + z²)，閉式：[−sqrt(1−2xz+z²)/z]_{−1}^{1} = (1 + z − |1 − z|)/z。
- `ex23_closed(y: float, z: float)` — 
- `ex23_coeffs(n: int)` — ex23_closed = Σ_k 2/(2k+1) (yz)^k → ∫P_n P_m dx = 2δ_{nm}/(2n+1)。
- `bessel_series(n: int, rho: float, terms: int=40)` — 
- `bessel_integral(n: int, rho: float)` — Ex. 2.4：J_n(ρ) = (1/2π)∫_{−π}^{π} e^{i(ρ sinθ − nθ)} dθ（取實部；虛部為 0）。
- `ideal_gas_Z(N: int, x)` — 
- `grand_Q(z: float, x: float)` — 
- `grand_probs(x: float, z: float, Nmax: int)` — 巨正則系綜裡粒子數 N 的機率 Z(N) z^N / Q（理想氣體 → Poisson(zx)）。
- `survival_prob(n: int, z)` — 
- `survival_mean(z)` — 
- `mortal_walkers(M: int, z: float, seed: int)` — 
- `mortal_trail(z: float, seed: int, nmax: int=10 ** 6)` — 
- `ellipk(k: float)` — K(k) = π / (2·AGM(1, sqrt(1−k²)))。
- `mortal_origin_exact(z: float)` — Σ_N z^N(1−z) C(N;0)/4^N = (1−z)·Σ_n C(2n,n)² (z/4)^{2n} = (1−z)(2/π)K(z)。
- `mortal_origin_series(z: float, nmax: int=4000)` — 
- `count1(N: int, d: int)` — 
- `count_grid(N: int, steps=SQ_STEPS)` — (2.7)：C(N; y) = Σ_i C(N−1; y − s_i)（最後一步來自鄰居），C(0; y) = δ_{y,0}。
- `count_grids(Nmax: int, steps=SQ_STEPS)` — 
- `count_sq(N: int, x: int, y: int)` — 
- `returns_sq(N: int)` — 
- `returns_sc(N: int)` — 
- `returns_bcc(N: int)` — 
- `returns_fcc(N: int)` — 
- `uv_walk(ucoins: str, vcoins: str)` — 
- `trail_sq(n: int, seed: int)` — 
- `chi_1d(q: float)` — 
- `chi_sq(qx: float, qy: float)` — 
- `chi_sc(qx: float, qy: float, qz: float)` — 
- `chi_bcc(qx: float, qy: float, qz: float)` — 
- `chi_fcc(qx: float, qy: float, qz: float)` — 
- `chi_from_steps(q, steps, scale: float=1.0)` — Σ_i e^{iq·w_i}，w_i = scale·steps[i]（對稱的鄰居 → 實數）。
- `g_fourier(z: float, chi: float)` — 
- `plane_wave_check(q: float, xs=range(-3, 4))` — max |cos(q(x−1)) + cos(q(x+1)) − 2cos q·cos(qx)|。
- `fourier_count_sq(N: int, x: int, y: int, M: int | None=None)` — (2.17)（二維、a = 1）：均分 M×M 點求積；M > N + |x| + |y| 時精確（三角多項式）。
- `fourier_count_sc(N: int, r=(0, 0, 0), M: int | None=None)` — 
- `chi2_terms()` — χ² = (Σ_i e^{iq·w_i})² 的 16 項：(第一步, 第二步, 是否為常數項＝回到原點)。
- `gauss_sc(N: float, r2: float=0.0)` — (2.20)，a = 1（書上版）。
- `gauss_sq(N: float, r2: float=0.0)` — 
- `gauss_lattice(N: int, kind: str, r2: float=0.0, parity: int=0)` — Z^N Ω (d/(2πN))^{d/2} e^{−d r²/(2N)} × 奇偶因子。r2 用最近鄰距離 1 的單位；
- `asym_returns_sq(N: int)` — 
- `asym_returns_1d(N: int)` — 
- `eq222(N: float, r: float, a: float, b: float)` — (2.22) 右邊（指數 −3r²/(2a²N)；書上印成 2aN）。
- `eq221_first_order(N: float, r: float, a: float, b: float)` — ∫ e^{−iqr − a²Nq²/6}(1 − bNq⁴) dq 的數值積分（核對 (2.22)）。
- `coin_quartic(N: int, d: int)` — 一維擲幣：cos^N q = e^{−Nq²/2 − Nq⁴/12 − …}，即 (2.21) 的 a² = 3、b = 1/12；兩個鞍點 → 乘 2/(2π)。
- `p_gauss_book(r: float, r0: float=1.0)` — 
- `p_gauss(r: float, r0: float=1.0)` — 
- `radial_integral(p, r0: float=1.0, R: float=60.0, n: int=20000)` — 
- `p_cauchy(r: float, r0: float=1.0)` — 
- `chi_cauchy(q: float, r0: float=1.0)` — 
- `chi_cauchy_numeric(q: float, r0: float=1.0, R: float=400.0, n: int=200000)` — (2.27) 第二行：(4/π)∫_0^∞ sin(qr) r r0 / (q (r²+r0²)²) dr。
- `P_cauchy_N(r: float, N: int, r0: float=1.0)` — 
- `P_cauchy_N_numeric(r: float, N: int, r0: float=1.0)` — (2.31) 第一行：(4π/(2π)³)∫_0^∞ e^{−Nqr0} q sin(qr)/r dq。
- `cauchy_radial_cdf(R: float, r0: float=1.0)` — 
- `cauchy_radial_median(r0: float=1.0)` — 
- `steps_gauss_3d(n: int, r0: float, rng: random.Random)` — 
- `steps_cauchy_3d(n: int, r0: float, rng: random.Random)` — 三維各向同性 Cauchy（多變量 t，自由度 1）：r0·Z/|W|，Z ~ N(0, I₃)、W ~ N(0,1)；密度正是 (2.26)。
- `trail_3d(steps)` — 
- `median_distance(kind: str, N: int, M: int, seed: int, r0: float=1.0)` — 
- `prob_grid_sq(N: int)` — 
- `diffuse_sq(P: dict)` — (2.32)，d = 2：P(N; y) = ¼ Σ_nn P(N−1; y − w_i)。
- `variance(P: dict)` — 
- `ex28_sides(N: float, r: float)` — P = (3/(2πN))^{3/2} e^{−3r²/(2N)}（a = 1）：∂P/∂N 與 (1/6)∇²P（徑向拉普拉斯 P'' + 2P'/r）。
- `g_cont(z: float, q: float, d: int=1)` — 
- `g_cont_book(z: float, q: float, d: int=1)` — 
- `G1_exact(z: float, x: int)` — 
- `G1_fourier(z: float, x: int, n: int=4000)` — 
- `G1_continuum(z: float, x: float)` — 
- `gauss_moment(n: int, a: float=1.0)` — 
- `gauss_moment_frac(n: int, a=1)` — (2:S1-6) 去掉 sqrt(π/a) 的有理部分 n!/((n/2)!(4a)^{n/2})。
- `gauss_moment_numeric(n: int, a: float=1.0)` — 
- `gauss_gen(a: float, b: float)` — 
- `square_complete(a: float, b: float, x: float, book: bool=False)` — 
- `dot(u, v)` — 
- `cross(u, v)` — 
- `triple(a1, a2, a3)` — 
- `reciprocal(a1, a2, a3)` — 
- `roots_sum(m: int, n: int)` — 
- `series_coeffs_power(alpha: float, zc: float, nmax: int)` — (zc − z)^{−α} = zc^{−α} Σ_n [α(α+1)…(α+n−1)/n!] (z/zc)^n。
- `coeff_power_exact(N: int, alpha: float, zc: float)` — 
- `coeff_power_asym(N: int, alpha: float, zc: float)` — 
- `coeff_two_poles(N: int, a: float, b: float, zc: float, zd: float)` — 
- `coeff_log_exact(N: int, alpha: float, zc: float)` — −(zc − z)^{−α} ln(zc − z)；ln(zc − z) = ln zc − Σ_{k≥1} (z/zc)^k / k。
- `digamma(x: float)` — 
- `coeff_log_asym(N: int, alpha: float, zc: float, book: bool=False, terms: int=3)` — (2:S3-11)：d/dα [zc^{−N−α} N^{α−1}/Γ(α)] = [−ln zc − ψ(α) + ln N]·zc^{−N−α}N^{α−1}/Γ(α)；
- `check()` — 

## `combkit.walk_rec`

walk_rec.py —— 《Elements of the Random Walk》第 3 章 Generating functions II（回家、維度、走過幾個點）的數學

常數：`EULER_GAMMA`、`Z`、`BIPARTITE`、`STEPS`、`_I0_ASYM`、`B1_SC`、`G_C`、`SIM_SEED`、`SIM_RETURNED`、`SIM_RETURNED_100`、`SIM_COUNTS_SEED`、`SIM_COUNTS_SC`、`SIM_DISTINCT_SEED`、`SIM_DISTINCT`

- `CT(m: int)` — 一維 m 步回到原點的走法數 C(m, m/2)（m 奇數為 0）。
- `walks_1d(N: int, x: int=0)` — 
- `visit_counts_1d(N: int, w: int, x: int=0)` — 
- `visit_counts_gf_1d(N: int, w: int, x: int=0)` — 
- `visit_poly(counts: dict)` — 
- `counts_from_poly(p)` — 
- `virial_terms_1d(N: int, w: int, x: int=0)` — 
- `never_return_1d(N: int)` — 
- `never_return_1d_brute(N: int)` — 
- `never_return_1d_gf(nmax: int)` — 
- `escape_1d(N: int)` — 
- `escape_1d_asym(N: float)` — 
- `escape_1d_book(N: float)` — 
- `escape_biased(p: float)` — 
- `escape_biased_series(p: float, N: int)` — [z^N] (1/(1−z))/G_p(z; 0)，G_p = Σ C(2k,k)(pq)^k z^{2k}（機率的生成函數，z_c = 1）。
- `escape_biased_dp(p: float, N: int)` — 
- `returns_count(kind: str, N: int)` — 
- `return_probs(kind: str, nmax: int)` — 
- `first_return_probs(kind: str, nmax: int)` — f_0..f_nmax：第一次回家在第 n 步的機率（更新方程 u_n = Σ_{k≤n} f_k u_{n−k}，f_0 = 0）。
- `no_return_probs(kind: str, nmax: int)` — 
- `return_count_dist(kind: str, N: int, mmax: int)` — N 步內剛好回家 m 次的機率（m = 0..mmax）：P(≥ m) = Σ_{t≤N} f^{*m}(t)。
- `never_return_dp(kind: str, N: int)` — 
- `no_return_sq_book(N: float)` — 
- `no_return_sq_refined(N: float)` — 
- `log10_steps_for(p: float)` — 
- `K_d(d: int)` — 
- `K_d_gauss(d: int)` — 
- `tail_integral(d: int, eps: float)` — 
- `tail_integral_numeric(d: int, eps: float)` — 
- `i0e(y: float)` — e^{−y} I_0(y)（y ≥ 0）：小 y 用級數、大 y 用漸近展開。
- `G0(kind: str, x: float)` — 
- `G0_series(kind: str, x: float, nmax: int=4000)` — 
- `sc_correction_coeff()` — 
- `watson_gamma(kind: str)` — 
- `watson_sc_bessel()` — 
- `hurwitz_zeta(s: float, a: float, M: int=12)` — ζ(s, a) = Σ_{j≥0} (a+j)^{−s}（Euler–Maclaurin）。
- `watson_bcc_series(K: int=200)` — Σ_{k≤K} c_k³ + 尾巴；c_k = C(2k,k)/4^k ≈ (πk)^{−1/2}(1 − 1/8k + 1/128k² + 5/1024k³)，
- `watson_fcc_series(nmax: int=60)` — Σ_{n≤nmax} u_n + (1/√2)(3/2π)^{3/2} ζ(3/2, nmax+1)（高斯尾巴；fcc 單胞體積 2、每座標變異數 2/3）。
- `escape(kind: str)` — 
- `f_m(kind: str, m: int)` — 
- `new_site_probs_brute(kind: str, n: int)` — 
- `distinct_sites_brute(kind: str, n: int)` — 
- `distinct_sites_1d(n: int)` — 
- `distinct_sites(kind: str, nmax: int)` — 
- `footprint_1d_asym(n: float)` — 
- `footprint_2d_book(n: float)` — 
- `footprint_2d_mw(n: float)` — 
- `simulate_first_returns(kind: str, M: int, N: int, seed: int)` — 
- `simulate_return_counts(kind: str, M: int, N: int, seed: int)` — 
- `simulate_distinct(kind: str, M: int, n: int, seed: int)` — 
- `trail(kind: str, n: int, seed: int)` — 
- `check()` — 

## `combkit.walk_rg`

walk_rg.py —— 《Elements of the Random Walk》第 11 章 Scaling, fractals, and renormalization 的數學

常數：`F`、`LN4_LN3`、`LN2_LN3`、`_ROT60`、`MASS_GEN`、`MASS_K`、`BOX_OFFSET`、`BOX_K`、`BOX_K_CANTOR`、`REPLICAS`、`RW_SEED`、`RW_NS`、`RW_M`、`SAW_NEFF`、`COLLAPSE_R`、`COLLAPSE_HMAX`、`COLLAPSE_POINTS`、`A_BOOK`、`COBWEB_R0`、`COBWEB_N`、`SMALL_R`、`SMALL_H`、`ISING_K0`、`ISING_FIELD`、`RING_NS`、`XI_K`、`FREE_K`、`ALPHA_N`、`FLOW_R0`、`FLOW_H0`、`FLOW_ELL`、`TOY_A`、`TOY_STARTS`、`RG_CASES`、`XI_R`

- `koch_points(gen, p0=(0.0, 0.0), p1=(1.0, 0.0))` — 第 gen 代 Koch 曲線的頂點 (4^gen + 1, 2)，從 p0 到 p1（突起在前進方向的左邊）。
- `koch_counts(gen)` — (段數 4^gen, 總長 (4/3)^gen, 每段長 3^−gen)（Fraction）。
- `cantor_intervals(gen)` — 第 gen 代 Cantor 集的區間 [(a, b)]（Fraction，2^gen 個）。
- `cantor_counts(gen)` — (段數 2^gen, 總長 (2/3)^gen, 每段長 3^−gen)（Fraction）。
- `koch_mass(R, gen=MASS_GEN)` — 以起點 (0, 0) 為圓心、半徑 R 的圓內 Koch 曲線的份量（第 gen 代，以段中點判斷；全長 = 1）。
- `cantor_mass(R, gen=12)` — 以左端點 0 為中心、半徑 R 的區間內 Cantor 集的份量（第 gen 代；全部 = 1）。
- `mass_table(kind='koch')` — [(k, R = 3^−k, n(R))]：Koch n = 4^−k、Cantor n = 2^−k（(11.1)：n(R/3) = n(R)/4）。
- `mass_radius_dim(kind='koch', Rs=None)` — ln n(R) 對 ln R 的最小平方斜率。Rs 預設 3^−k（k = 0..4，正好 ln4/ln3、ln2/ln3）；
- `box_counts(kind='koch', ks=BOX_K, base=3)` — [(ε, 盒子數)]：邊長 ε = base^−k 的方格（Cantor 用線段）裡有曲線的個數（格線位移 BOX_OFFSET）。
- `box_dim(kind='koch', ks=None, base=None)` — 盒子計數的維度：ln N(ε) 對 ln(1/ε) 的斜率（第二種算法；Koch 1.2619、Cantor 0.6309 附近）。
- `replica_dim(nrep, mag)` — (11.6)：p = ln ℛ / ln ℳ（ℛ 個複本，每個放大 ℳ 倍變回原物）。
- `rw_paths(Nmax=RW_NS[-1], M=RW_M, seed=RW_SEED)` — M 條二維正方格子普通漫步（每步 E N W S 各 1/4），形狀 (M, Nmax+1, 2)。
- `rw_rg2(N, M=RW_M, seed=RW_SEED)` — 前 N 步（N+1 個點）的迴轉半徑平方 R_g² 的 M 條平均（固定種子）。
- `rg2_exact(N)` — ⟨R_g²⟩ = N(N+2)/(6(N+1))（每步長 1，任何維度；Σ_{i<j}|i−j| / (N+1)²）。
- `rw_df(Ns=RW_NS, exact=False)` — 普通漫步的質量－半徑：ln N 對 ln R_g 的斜率（→ 2）。exact=True 用 rg2_exact。
- `rw_trail(N=1500, seed=RW_SEED)` — 畫圖用的一條普通漫步（N+1 個格點座標）。
- `saw_df(d)` — 自迴避漫步的碎形維度 1/ν（d = 1、2、4 是 Fraction；三維 ν = 0.588 只引用）。
- `saw_df_eff(N)` — 有限 N 的有效維度 2/(2ν_eff)（walk_saw.two_nu_eff，二維精確列舉）。
- `landau_m(r, h, u=1.0)` — F/V = r m²/2 + u m⁴/4 − h m 的平衡磁化（h ≥ 0；u m³ + r m = h 的實根中 F 最小者）。
- `landau_m_bisect(r, h, u=1.0)` — 同上，用二分法（第二種算法；h > 0 時根唯一的正根，或 r < 0、h = 0 時的正根）。
- `landau_scaled(x, sign=1)` — 標度函數 f(x)：y = m/|r|^{1/2}、x = h/|r|^{3/2}（u = 1）滿足 sign·y + y³ = x（(11.16)）。
- `landau_curve(r, n=40, hmax=COLLAPSE_HMAX)` — 畫圖用 (h, m(r, h))，n 點。
- `collapse_curve(r, n=40, hmax=COLLAPSE_HMAX)` — 同上換成 (h/r^{3/2}, m/r^{1/2})（r > 0）。
- `mf_rescale_check(r=0.3, h=0.2, d=3, bs=(0.5, 2.0, 3.0))` — (11.10)–(11.11)：m(r, h) = b^{−d/4} m(r b^{d/2}, h b^{3d/4})（u 不變）的最大誤差。
- `book1116_limit(r, h)` — p. 263：書上說 f(h/r^{3/2}) → h^{1/3} r；應為 → h^{1/3} r^{−1/2}（即 f(x) ≈ x^{1/3}）。
- `gauss_m(r, h)` — Gaussian 模型（只有 q = 0 的模態受外場）：m = h/r（r > 0）。
- `gauss_rescale_check(d, r=0.3, h=0.2, bs=(0.5, 2.0, 3.0))` — (11.22)：m = b^{(2−d)/2} m(r b², h b^{(d+2)/2}) 的最大誤差（m = h/r）；(11.24) 的形狀 r^{(d−2)/4}F(h/r^{(d+2)/4}) 也核對。
- `gauss_exponents(d)` — Ex. 11.4：由 (11.24) 讀出 β = (d−2)/4、γ = 1、δ = (d+2)/(d−2)、α = 2 − d/2、ν = 1/2（Fraction）。
- `scaling_combos(d)` — (Landau 的 h/r^{3/2} 指數, Gaussian 的 h/r^{(d+2)/4} 指數, Landau 前因子 1/2, Gaussian 前因子 (d−2)/4)：d = 4 時兩兩相同。
- `square_check(r=0.5, a=1.0, S2=0.3, S4=-0.7, h=0.0)` — (11.28)→(11.30)：∫dS3 exp(−r/2 S3² − a/2(S3−S2)² − a/2(S3−S4)² + h S3) 的 (數值積分, 完成平方的閉式)。
- `square_book_form(r=0.5, a=1.0, S2=0.3, S4=-0.7)` — (11.30) 照書上的寫法：√(2π/(r+2a)) exp[−a²(S2−S4)²/(2(r+2a)) − ar(S2²+S4²)/(2(r+2a))]（對的）。
- `decimate_gauss(r, a=A_BOOK, h=0.0)` — 抽掉奇數格點、還沒重新標度：(r', a', h')＝(r(r+4a)/(r+2a), a²/(r+2a), h(r+4a)/(r+2a))（(11.32)(11.33)）。
- `spin_factor(r, a=A_BOOK)` — (11.34) 的自旋因子 √((r+2a)/a)（舊自旋 = 因子 × 新自旋；小 r → √2）。
- `rr_gauss(r, a=A_BOOK)` — (11.36)：r'' = r(r+4a)/a（a'' = a）。
- `hh_gauss(r, h, a=A_BOOK)` — (11.38)：h'' = h(r+4a)/√(a(r+2a))（小 r → 2^{3/2} h）。
- `gauss_rescaled(r, a=A_BOOK, h=0.0)` — decimate_gauss 再乘自旋因子的平方 / 一次方（第二種算法）→ (r'', a'', h'')。
- `ring_matrix(N, r, a=A_BOOK)` — 
- `ring_lnZ(N, r, a=A_BOOK, h=0.0)` — N 個自旋的 Gaussian 環（週期，N ≥ 3）：ln Z = (N/2)ln 2π − ½ ln det M + ½ h²·1ᵀM⁻¹1。
- `ring_decimation_check(N=16, r=0.5, a=A_BOOK, h=0.2, book=False)` — ln Z_N(r, a, h) ＝ (N/2)[½ ln(2π/(r+2a)) + h²/(2(r+2a))] + ln Z_{N/2}(r', a', h') 的誤差（精確，不必小 r）。
- `gauss_m_invariance(r, h, a=A_BOOK)` — m = h/r 與「自旋因子 × h''/r''」（一次抽格＋重新標度）相同：回傳差。
- `book1147_m(r, h, n)` — 照書上 (11.47) m = 2^{−n/2} f(r 2^{2n}, h 2^{3n/2})（f 用精確的 h/r）；正確是 2^{+n/2}。回傳 (書上, 更正, 精確 h/r)。
- `gauss_iterate(r0, n, a=A_BOOK)` — [r^(0), …, r^(n)]（(11.39)）。
- `cobweb_points(r0=COBWEB_R0, n=COBWEB_N, a=A_BOOK)` — Fig. 11.7 的蛛網折線：(r0, 0) → (r0, r1) → (r1, r1) → (r1, r2) → …（2n 個轉角＋起點）。
- `rr_curve(n=40, xmax=None, a=A_BOOK, ymax=10.0)` — Fig. 11.6 的曲線 (r, r'')，r 從 0 到 r''= ymax 的地方（a = 1：r = −2 + √14 = 1.7417）。
- `small_r_table(n=5, r=SMALL_R, h=SMALL_H, a=A_BOOK)` — [(k, r^(k)/(4^k r), h^(k)/(2^{3k/2} h), 自旋總因子/2^{k/2})]（(11.41)–(11.45)；小 r 時都 ≈ 1）。
- `gauss_fixed_points(a=A_BOOK)` — r'' = r 的根：0（不穩定，斜率 4 = 2^{Λ_r}）與 −3a（r + 2a < 0，積分發散，不物理）。
- `chain_lambdas()` — 抽格 b = 2 的特徵值 → Λ = ln(特徵值)/ln b：r 的 4 → 2、h 的 2^{3/2} → 3/2、自旋 √2 → 1/2 = (2−d)/2（d = 1）。
- `ising_A(K)` — (11.54)：A = 2√cosh 2K。
- `ising_Kp(K)` — (11.55)：K' = ½ ln cosh 2K。
- `ising_pair_check(K)` — (11.52)–(11.57)：Σ_{S3} e^{K S3(S2+S4)} 與 A e^{K' S2 S4} 在四種 (S2, S4) 的最大差。
- `book1153(K)` — (11.53) 照書上：左邊 cosh K(S2+S4)（同號時 cosh 2K），右邊 2cosh 2K → (左邊, 書上的右邊, 正確左邊 2cosh K(S2+S4))。
- `ex118(K, H)` — Ex. 11.8：Σ_{S3} e^{S3(K(S2+S4)+H)} = A e^{K1' S2S4 + K2'(S2+S4)}
- `ex118_check(K, H)` — Ex. 11.8 的四種 (S2, S4) 直接加總 vs 公式的最大差。
- `ising_iterate(K0=ISING_K0, n=5)` — (11.62)：K_{n+1} = ½ ln cosh 2K_n 的 [K_1 = K0, K_2, …]（→ 0）。
- `ring_Z_brute(N, K, H=0.0)` — N 個 ±1 自旋的環：逐一列舉 2^N 個組態（N ≤ 16）。
- `ring_Z_transfer(N, K)` — 轉移矩陣：Z = (2cosh K)^N + (2sinh K)^N。
- `ring_Z_decimation(N, K, book=False)` — (11.63)：N = 2^k 的環一直抽格到剩 1 個自旋：Z = Π_{j=1}^{k} A(K_j)^{N/2^j} · 2e^{K_{k+1}}
- `x_of_K(K)` — x = e^{−2K}（Ex. 11.9）。
- `x_map(x)` — Ex. 11.9：x' = 2x/(1+x²)（由 e^{2K'} = cosh 2K）。
- `x_iterate(x0, n)` — 
- `x_fixed_points()` — x' = x 的根：x = 0（T = 0，不穩定，斜率 2）與 x = 1（T = ∞，穩定，斜率 0）（x = −1 不物理）。
- `x_slope(x)` — dx'/dx = 2(1−x²)/(1+x²)²：x = 0 → 2、x = 1 → 0。
- `xi_exact(K)` — 一維 Ising 的關聯長度 ξ = −1/ln tanh K（＝ walk_crit.ising1d_xi）。
- `xi_ratio(K)` — ξ / e^{2K}（→ ½；Ex. 11.9(b)：ξ ∝ 1/x = e^{2J/k_BT}）。
- `tanh_rg_check(K)` — tanh K' = tanh² K（精確）⇒ ξ(K') = ξ(K)/2：回傳 (|tanhK' − tanh²K|, ξ(K)/ξ(K') − 2)。
- `free_energy_ratio(K, N=None)` — −F/(NJ) = ln Z/(NK)：N = None 用無限長 ln(2cosh K)/K；N = 2^k 用抽格的連乘（→ 1：F → −NJ）。
- `free_energy_series(K, jmax=60)` — 無限長鏈：ln Z/N = Σ_j 2^{−j} ln A(K_j)（(11.63) 取 N → ∞；第二種算法）。
- `alpha_q(q, a=A_BOOK, c=1.0)` — 一維 Gaussian 鏈的 Fourier 係數：(a/2)Σ(S_i − S_{i+1})² = Σ_q a(1 − cos qc)|S(q)|²；小 q：a c² q²/2。
- `alpha_q_book(q, a=A_BOOK, c=1.0)` — p. 273 書上 α(q) = −2a(1 − cos qc)（負號與 2 倍都不對）。
- `alpha_table(N=ALPHA_N, a=A_BOOK)` — [(k, q = 2πk/N, a(1−cos q), aq²/2, 書上)]，k = 0..N/2。
- `fourier_check(N=ALPHA_N, r=0.5, a=A_BOOK, seed=1105)` — 隨機實數自旋：Σ_i [r/2 S_i² + a/2(S_i−S_{i+1})²] ＝ Σ_q [r/2 + a(1−cos q)]|S(q)|²（S(q) = N^{−1/2}Σ S_i e^{−iqi}）的相對誤差；
- `shell_fraction(d, b=2)` — 布里淵區每邊縮成 1/b：拿掉的模態比例 1 − b^{−d}（d = 2：3/4）。
- `shell_rg(r, h, d, b=2)` — (11.71)(11.72) 與 b 的推廣：(b² r, b^{(d+2)/2} h)。
- `shell_rg_steps(r, h, d, b=2)` — 一步一步：q → q/b（Jq² → Jq²/b²）、N → b^d N、S → bS（保持 J）→ (r', h')（第二種算法）。
- `flow(r0, h0, d, ell)` — (11.82)(11.83) 的解：r(ℓ) = r0 e^{2ℓ}、h(ℓ) = h0 e^{(d+2)ℓ/2}。
- `flow_numeric(r0, h0, d, ell)` — 同上，solve_ivp（RK45）數值積分（第二種算法）。
- `flow_line(r0, h0, d, ell_max, n=40)` — 畫圖用流線 [(r(ℓ), h(ℓ))]，ℓ ∈ [0, ell_max]，n 點（h ∝ |r|^{(d+2)/4}）。
- `m24(r, h, d, Ffn=None)` — (11.24)：m = r^{(d−2)/4} F(h/r^{(d+2)/4})（預設 F(x) = x/(1+x²)，任意可微的 F 都行）。
- `pde_residual(d, r=0.4, h=0.3, eps=1e-05, Ffn=None)` — (11.87)：[−(d−2)/2 + 2r ∂_r + (d+2)/2 h ∂_h] m 用中央差分（m = (11.24)）；回傳殘差 ÷ |m|。
- `char_check(d, r=0.4, h=0.3, ells=(0.5, 1.0, 2.0), Ffn=None)` — 特徵線：e^{(2−d)ℓ/2} m(r(ℓ), h(ℓ)) 與 ℓ 無關（(11.86)）：最大相對差。
- `toy_eigen(A=TOY_A)` — 2×2 Fraction 矩陣的特徵值與特徵向量（特徵多項式 λ² − tr λ + det，判別式要是完全平方）：
- `toy_flow(c1, c2, ell)` — (11.95)：δκ(ℓ) = c1 e^{2ℓ}(2, 1) + c2 e^{−3ℓ}(1, −2)。
- `toy_flow_expm(c1, c2, ell)` — 同上，用矩陣指數 e^{Aℓ} δκ(0)（第二種算法）。
- `toy_flow_line(c1, c2, ell_max=1.0, n=30)` — 畫圖用流線（n 點）；c1 = 0 的線留在臨界線（無關方向）上流進不動點。
- `rg_exponents(Lr, Lh, d)` — (11.105)–(11.108)＋ν = 1/Λ_r → dict(alpha, beta, gamma, delta, nu)（Fraction）。
- `rg_exponents_widom(Lr, Lh, d)` — 第二種算法：α = 2 − d/Λ_r、能隙 Δ = Λ_h/Λ_r 代入 walk_crit.widom_exponents（第 9 章 Ex. 9.1–9.2）。
- `lambdas_from(nu, beta, d)` — 反過來：Λ_r = 1/ν、Λ_h = d − β/ν（第三種核對）。
- `gauss_lambdas(d)` — (11.109)(11.110)：Λ_r = 2、Λ_h = (d+2)/2。
- `xi_from_flow(r, Lr, D=1.0, Kc=1.0)` — (11.99)–(11.103)：K r e^{Λ_r ℓ*} = D ⇒ ξ = e^{ℓ*} = (K r/D)^{−1/Λ_r} → (ℓ*, ξ)。
- `xi_from_flow_numeric(r, Lr, D=1.0, Kc=1.0)` — 同上，數值積分 dδ/dℓ = Λ_r δ（δ(0) = K r）直到 δ = D（事件）（第二種算法）。
- `summary()` — 
- `check()` — 

## `combkit.walk_saw`

walk_saw.py —— 《Elements of the Random Walk》第 7 章 Path integrals and self-avoidance 的數學

常數：`F`、`STEPS2`、`EX_SPLIT`、`EX71`、`SAW_KNOWN_C`、`SAMPLE_SEEDS`、`SECOND_ORDER`、`EX_CONV`、`DYSON_TOY`、`MF_CASES`、`NU_QUOTED`

- `count1(N, x)` — 一維 N 步從 0 到 x 的走法數。
- `count2(N, x, y)` — 二維正方格子 N 步從原點到 (x, y) 的走法數（u = x+y、v = x−y 兩個獨立的一維）。
- `count_brute(N, end, dim)` — 暴力列舉 (2·dim)^N 條路線，數終點 = end 的條數。
- `split_terms(N, M, end, dim)` — (7.3) 的每一項：(中間點 r, C(M; 0→r), C(N−M; r→end))，只列乘積非零者。
- `split_sum(N, M, end, dim)` — 
- `multi_split(N, Ms, end, dim)` — (7.5)：依序走 Ms[0], Ms[1], … 步，對每個中間點求和（逐段摺積的分布）。
- `gauss1(M, x, d=3, a=1.0, book=False, book711=False)` — (7.6) 一段高斯的一個直角分量。
- `gauss_conv(M1, M2, x, **kw)` — 
- `conv_ratio(M1, M2, x=0.7, book=False)` — 
- `lattice_var(M, d=3)` — 簡單立方格子 M 步後 x 分量的變異數：x 分量每步以 1/d 的機率 ±1。
- `saw_enum(Nmax=14)` — 二維自迴避漫步：c_N 與 Σ R²（N = 0..Nmax）。對稱：第一步向東（×4），第一次轉彎向北（×2）。
- `saw_c(N)` — 
- `saw_R2(N)` — 
- `saw_frac(N)` — 
- `all_walks_xy(N)` — 
- `n_coinc(N)` — 
- `saw_brute(N)` — 
- `rw_R2_brute(N)` — 
- `saw_list(N)` — 
- `sample_rw(N, seed=702)` — 
- `sample_saw(N, seed=702)` — 
- `two_nu_eff(N)` — 
- `mu_eff(N)` — 
- `ie_terms(N)` — 
- `ie_terms_subsets(N)` — 
- `order1_formula(N)` — 
- `ie_partial(N, terms=None)` — 
- `weighted_count(N, f)` — 
- `order2_types(N)` — 
- `matchings(k)` — 
- `chord_count(k)` — 
- `is_irreducible(m)` — 
- `irreducible_count(k)` — 
- `reversal_classes(k)` — 
- `dyson_counts(kmax=4)` — 
- `convolve(a, b)` — 
- `gf_value(coeffs, z)` — 
- `gf_check(z=0.1, n=80)` — 
- `fourier_conv_error(n=16, seed=7)` — 
- `dyson_partial(g0, sigma, n)` — 
- `dyson_closed(g0, sigma)` — 
- `dyson_book35(g0, sigma)` — 
- `dyson_book33(g0, sigma)` — 
- `g0_1d(z, q)` — 
- `sigma1_1d(z)` — 
- `sigma1_1d_quad(z)` — 
- `g_dyson_1d(z, q, f, n=None)` — 
- `power_exp(l, d)` — 
- `power_exp_count(l, d)` — 
- `power_table(dmax=5, lmax=4)` — 
- `loop_slope(case, e1=1e-06, e2=1e-08)` — 
- `G_MF(x, Ns, f)` — 
- `G_MF_quad(x, Ns, f)` — 
- `cmf_ratio(m, Ns, f, method='quad')` — Ns z_c^m C_MF(m) = (1/m!) ∫_0^∞ t^m e^{−t − f t²/Ns} dt（t = Ns λ）。
- `cmf_asym(m, Ns, f)` — 
- `cmf_book59(m, Ns, f)` — 照書上 (7.59)：指數 −Ns(fλ² + Ns λ) + m ln λ，乘 Ns^{m+1}/m!（對數空間數值積分）。
- `ex75(Ns, x=0.5, f=1.0)` — 
- `ex75_big(Ns, x=1.5, f=1.0)` — 
- `gmf_series(x, Ns, f, mmax=400)` — 
- `eq755(z0, z, zc)` — 
- `G1_lattice(z, x)` — 
- `eq747_residual(z, x, book=False)` — 一維、a = 1、z_c = 1/2、δ 換成 δ_{x,0}：更正版 (1−z/z_c)G − zΔG − δ；書上 zΔG − (z/z_c − 1)G − δ。
- `eq748_eigen(q, book=False, x=0.3, h=0.001)` — 
- `flory_F(R, m, d, u=1.0, c=1.0)` — (7.63) 的指數；c = 熵那一項的係數（書上 1；高斯的正確係數 d/2 不改變指數）。
- `flory_Rstar(m, d, u=1.0, c=1.0)` — 
- `flory_Rstar_num(m, d, u=1.0, c=1.0)` — 
- `flory_nu(d)` — 
- `flory_slope(d, m1=1000.0, m2=1000000.0, c=1.0)` — 
- `summary()` — 
- `check()` — 

## `combkit.walk_scale`

walk_scale.py —— 《Elements of the Random Walk》第 8 章 Properties of the random walk: introduction to scaling 的數學

常數：`WIDTH_N`、`WIDTH_M`、`WIDTH_SEED`、`TRAIL_SEED`、`SUBST_Y`、`SUBST_EPS`、`COLLAPSE_EPS`、`COLLAPSE_X`、`XI_CASES`、`SD_N`、`SD_X`、`FISHER_NU`、`ETA_2D`、`NU_2D`、`MU_QUOTED`、`GAMMA_QUOTED`、`RATIO_N`、`SA_N`、`UPPER_CRITICAL_D`

- `pareto_alpha(y, d=3)` — 
- `delta_of(y)` — 
- `moment2_trunc(y, R)` — (8.2)：⟨r²⟩_R = ∫_1^R y r^{−1−y} r² dr。
- `moment2_trunc_quad(y, R)` — 
- `chi_pareto(q, y, d=3)` — χ(q) = ⟨e^{iq⃗·r⃗}⟩：d = 3 平均 sin(qr)/(qr)、d = 1 平均 cos(qr)（QAWF 振盪積分）。
- `B_pareto(y, d=3)` — 
- `B_pareto_quad(y, d=3)` — B = y ∫_0^∞ u^{−1−y}(1 − sin u/u) du（d = 3）或 y ∫_0^∞ u^{−1−y}(1 − cos u) du（d = 1）。
- `chi_small(q, y, d=3)` — 
- `chi_ratio(q, y, d=3)` — 
- `chi_slope(y, d=3, q=0.001)` — 
- `sample_steps(rng, M, y, d=3)` — 
- `median_width(N, y, M=WIDTH_M, seed=WIDTH_SEED)` — 
- `width_pred(N, y)` — 
- `width_slope(y)` — 
- `trail(N, y, seed=TRAIL_SEED)` — 
- `levy_density(r, s, y)` — 
- `levy_tail(r, s, y)` — 
- `levy_slope(r, s, y, h=1.05)` — 
- `tail_exponent(delta, d=3)` — 
- `subst_integral(r, y, eps)` — 
- `subst_closed(r, y, eps)` — 
- `subst_book(r, y, eps)` — 照書上 w = k/r（k = rw、dk = r dw）：r^{y+1}∫ e^{i r² w} w^y e^{−ε r w} dw——r 仍留在指數裡。
- `zc(d)` — 
- `AB(d)` — 
- `chi_axis(q, d)` — 
- `g_lattice(z, q, d=1)` — 
- `g_88b(z, q, d=1)` — 
- `g_88c(z, q, d=1)` — 
- `collapse(eps, x, d=1)` — 
- `collapse_limit(x, d=1)` — 
- `g_general(tau, q, nu, eta)` — 
- `gamma_numeric(nu, eta, q=1e-12)` — 
- `G_lattice_sc(z, n)` — ∫_0^∞ e^{−t} I_n(2zt) I_0(2zt)² dt（1/(1 − zχ) = ∫ e^{−t(1−zχ)} dt，每個方向 ∫ e^{2zt cos q} e^{iqn} dq/2π = I_n(2zt)）。
- `xi_88c(x)` — 
- `xi_lattice(x)` — 
- `F3(u)` — 
- `G_yukawa(x, r)` — 
- `G_ratio(x, n)` — 
- `G_yukawa_ft(x, r)` — (8.8c) 的數值逆 Fourier：(1/(2π²r))∫ q sin(qr) (a/b)/(q² + m²) dq（第二種算法；m = 1/ξ）。
- `ratio_exact(N, x)` — 
- `sd_K(nu, B, zc_, book=False)` — 
- `sd_ratio(N, r, nu, B, zc_, book=False)` — 
- `sd_delta(N, r, nu, B, zc_)` — 
- `sd_delta_num(N, r, nu, B, zc_)` — 
- `sd_phi_num(N, r, nu, B, zc_)` — 
- `sd_K_dim(d)` — 
- `fisher(nu)` — 
- `fisher_num(nu, N=10 ** 6, r=1000.0)` — 
- `gamma_exp(eta, nu)` — 
- `growth(gamma, zc_)` — 
- `log_coeff_zN(N, gamma, zc_)` — ln[z_c^N × (z^N 在 (z_c − z)^{−γ} 的係數)] = ln Γ(N+γ) − ln Γ(γ) − ln N! − γ ln z_c。
- `p_numeric(gamma, zc_, N=1000)` — 
- `ratio_method(N, mu=MU_QUOTED)` — 
- `theta(gamma, nu)` — 
- `ex85_s(d, nu)` — 
- `return_slope(d, N)` — 
- `kin_factor(l1, l2)` — 
- `int_factor(l1, l2, d)` — 
- `int_factor_book86(l1, l2, d)` — 
- `rescale_check(l1, l2, d, n=200, w=0.5, seed=804)` — 離散高斯鏈 r_i（s_i = iΔs）；K = Σ|Δr|²/Δs、V = Σ_{i≠j} Δs² g_w(r_i − r_j)（g_w 是寬 w 的歸一高斯＝δ 的近似）。
- `sa_amplitude(N, d)` — 
- `flory_l1(l2, d)` — 
- `flory_l1_num(l2, d)` — 
- `nu_rescale(d)` — 
- `ex88_factor(l1, l2, q)` — 
- `ex88_amplitude(N, q)` — 
- `ex88_relevant(d, q)` — 
- `ex88_nu(q)` — 
- `summary()` — 
- `check()` — 

## `combkit.walk_shape`

walk_shape.py —— 《Elements of the Random Walk》第 6 章 The shape of a random walk 的數學

常數：`SEED_CLOUD`、`SEED_FIG63`、`SEED_FIG64`、`SEED_FIG65`、`SEED_GAUSS`、`SEED_D`、`N_D`、`M_D`、`ZETA3`、`SMALL_STEPS`

- `sites(steps)` — 
- `center_of_mass(pts)` — 
- `tensor(pts)` — 
- `tensor_pairs(pts)` — (6:S2-2b)：T = (1/(2n²)) Σ_i Σ_j (r_i − r_j)(r_i − r_j)ᵀ，n = 點數。
- `com_frac(pts)` — 
- `tensor_frac(pts)` — 
- `kramers(N)` — 
- `kramers_np(N)` — 
- `kramers_book(N, kind='nodiag')` — 照書上：'nodiag' (6.22) 只寫 α<β 與 α>β（對角線 0）；'noscale' (6.23) x<(N+1−x>)（少 1/(N+1)²）。
- `kramers_count(N)` — 補充 3：(r_k − r_j)² 展開後含 η_m η_n（m ≤ n）的點對 j < k 個數（點編號 0..N）；回傳 N×N 整數矩陣。
- `tensor_kramers(steps)` — 
- `tensor_kramers_frac(steps)` — 
- `principal(T)` — 
- `rotate(pts, theta)` — 
- `extent(T, phi)` — 
- `small_example()` — 
- `lattice_steps(M, N, d, seed)` — 
- `batch_tensor(steps)` — 
- `batch_eigs(T)` — 
- `cloud(M=1000, N=100, seed=SEED_CLOUD)` — 
- `cloud_stats()` — 
- `fig63(n=6, N=1000, seed=SEED_FIG63)` — 
- `fig63_ratios()` — 
- `fig64(N=1000, seed=SEED_FIG64)` — 
- `fig64_stats(theta=math.pi / 6)` — 
- `invariants(T)` — 
- `det_from_traces(T)` — 
- `asph(vals)` — 
- `asph_book16(vals, N)` — 書上 (6.16) 括號內照字面：分母 (d−1) Σ_{k=1}^{N} R_k²（沒平方；k > d 的 R_k 當 0）。
- `A_d(d)` — 
- `A_from_moments(d, t11sq, t1122, t12sq)` — 
- `A_S21(d, book=False)` — (6:S2-1) 用補充 2 的比值 ⟨T11²⟩/⟨T11T22⟩ = 9/5、⟨T11²⟩/⟨T12²⟩ = 9/2。book=True 照書上（分母少乘 ⟨T11T22⟩/⟨T12²⟩）。
- `trace_a(N)` — 
- `trace_a2(N)` — 
- `diag_a2(N)` — 
- `moments_exact(N, d, kind='gauss')` — 
- `moments_enum(N, d)` — 暴力列舉全部 (2d)^N 條格子漫步（精確 Fraction）：(⟨T11²⟩, ⟨T11T22⟩, ⟨T12²⟩)。
- `A_exact(N, d, kind='gauss')` — 
- `M_p(p)` — 
- `Aprime(d)` — 
- `Aprime_laplace(d)` — ⟨TrT²/(TrT)²⟩ = ∫_0^∞ t⟨TrT² e^{−t TrT}⟩dt；高斯：傾斜後每個模態的變異數變小，Z = (y/sinh y)^{d/2}，y = √(2t/d)。
- `Aprime_closed(d)` — 
- `sim3d(M=20000, N=100, seed=SEED_FIG65)` — 
- `sim3d_stats()` — 
- `asph_hist(a, bins=20)` — 
- `staircase(N)` — 
- `staircase_com(N)` — 
- `eigs_exact(N)` — 
- `eigs_asym(N, k)` — 
- `sine_mode(N, k)` — 
- `staircase_asph(N)` — 
- `orth_fraction(N, d)` — 
- `orth_approx(N, d)` — 
- `gauss_moment(m)` — 
- `wick_pairings(m)` — 
- `wick_count(m)` — 
- `closed_cov(N)` — 
- `sim_gauss_moments(N=100, d=3, M=20000, seed=SEED_GAUSS)` — 
- `resolvent_weight(eps=0.001)` — 
- `sum_identity(n, M=200000)` — 
- `shift_sum_exact(i, N=4000)` — 
- `shift_1d(d)` — 
- `mean_1d(k, d, N=N_D)` — 
- `var_1d(k, d, N=N_D)` — 
- `eig_sim(d, N=N_D, M=None, seed=None)` — 
- `xi_identity_error(N=12, d=5, seed=1)` — 
- `onedee_table(d)` — 
- `gamma_density(x, alpha, d)` — 
- `fixman_K(s, d, book=False)` — 
- `fixman_K_product(s, d, n_max=200000, book=False)` — 
- `fixman_density(x, d, book=False, h=0.002)` — (6.74)：P(x) = (1/2π)∫e^{isx}K(s)ds = (1/π)Re∫_0^∞ e^{isx}K(s)ds，s = u²。
- `fixman_mean_var(d, N=None)` — 
- `fixman_moments_numeric(d, xmax=1.6, n=801)` — 
- `A_saw(d)` — 
- `f4(x)` — 
- `summary()` — 
- `check()` — 

## `combkit.walk_var`

walk_var.py —— 《Elements of the Random Walk》第 5 章 Variations on the random walk 的數學

常數：`PI`、`_XP`、`P0`

- `talbot(F, t: float, M: int=32)` — 數值 Laplace 反轉 f(t)（Abate–Valkó 固定 Talbot；F 在圍道內只要奇點都在左邊、λ 大時解析）。
- `biased_pmf(N: int, p)` — (5.1)：{x = 2n − N: C(N,n) p^n (1−p)^{N−n}}。
- `biased_dp(N: int, p)` — (5.2) 逐步遞迴：P(N+1; y) = p P(N; y−1) + (1−p) P(N; y+1)。
- `moments(P: dict)` — 
- `biased_mean(N, p, l=1)` — 
- `biased_var(N, p, l=1)` — 
- `biased_gauss(x, N, p, l=1.0)` — (5.7)：exp[−(x − (2p−1)Nl)²/(8Nl²p(1−p))]/√(8πNp(1−p)l²)。
- `taylor_residual(p: float, book: bool=False, l: float=1.0, s: float=30.0, y: float=7.0)` — (5.3)：|pP(y−l) + (1−p)P(y+l) − P(y) − [c ∂P + (l²/2)∂²P]|，P = e^{−y²/2s²}；c = l(1−2p)（更正）或 l(2p−1)（書上）。
- `diffusion_D(p, l=1.0, kind='careful')` — Ex. 5.1：'naive' l²/2（直接把 P(N+1) − P(N) 換成 ∂P/∂N）；'careful' 加上 −½∂²P/∂N² ≈ −(v²/2)∂²P/∂y² → 2p(1−p)l²。
- `width_ratio(p)` — 
- `shift_pmf(N: int, b: float)` — (5.10) 一維：每步先走 b、再 ±1（公正）→ [(x, P)]，x = k + Nb，P = C(N,(k+N)/2)/2^N。
- `shift_dp(N: int, b: float)` — 
- `shift_fourier(N: int, b: int, x: int)` — (5.10)：(1/2π)∫_{−π}^{π} χ(q)^N e^{iqx} dq，χ = e^{−iqb} cos q（(5.9)；b 整數，均分求積精確）。
- `persistent_dp(N: int, Q)` — (5.11)(5.12)：{x: (L, R)}，L/R = 最後一步向左／向右；N = 0 時 L = R = ½ 在 0。
- `persistent_pmf(N: int, Q)` — 
- `persistent_np(N: int, Q: float)` — numpy 版：(xs, L, R)，xs = −N..N。
- `persistent_var(N: int, Q)` — 
- `persistent_var_formula(N: int, Q)` — 
- `persistent_var_asym(N, Q)` — 
- `persistence_length(Q)` — 
- `T_matrix(q: float, Q: float)` — 
- `T_total(N: int, q: float, Q: float)` — 
- `eigs(q: float, Q: float)` — 
- `eig_residual(lam, q, Q, book=False)` — 
- `lambda1_small(q, Q)` — 
- `total_fourier(N: int, q: float, Q: float, part: str='LR')` — 
- `total_eig(N: int, q: float, Q: float)` — 
- `recursion2_residual(N: int, q: float, Q: float, book: bool=False)` — 
- `gauss_persistent(x, N, Q)` — 
- `tele_cont(x, t, v=1.0, gam=1.0)` — (5.61) 更正後的連續部分（向量）：(γe^{−γt}/2v)[I0(γs) + t I1(γs)/s]，s = √(t² − x²/v²)。
- `tele_book(x, t, v=1.0, gam=1.0)` — 照書上 (5.60)(5.61)：前面是 ½ 而不是 1/(2v)。
- `tele_spike(t, gam=1.0)` — 
- `tele_mass(t, v=1.0, gam=1.0, book=False, n=2000)` — 
- `tele_var(t, v=1.0, gam=1.0)` — 
- `tele_var_numeric(t, v=1.0, gam=1.0, n=2000)` — 
- `tele_fourier(q, t, v=1.0, gam=1.0)` — σ(q, t) = e^{−γt}[cosh ωt + γ sinh(ωt)/ω]，ω² = γ² − q²v²（(5.52) 的 Laplace 反轉）。
- `tele_fourier_numeric(q, t, v=1.0, gam=1.0, n=4000)` — 
- `tele_laplace(q, lam, v=1.0, gam=1.0)` — 
- `tele_D(v=1.0, gam=1.0)` — 
- `tele_lp(v=1.0, gam=1.0)` — 
- `tele_gauss(x, t, v=1.0, gam=1.0)` — 
- `tele_lattice(t, dt, v=1.0, gam=1.0)` — Ex. 5.3 的格子版：步長 l = vδt、Q = 1 − γδt、N = t/δt 步 → (xs, P)（xs 是實際座標）。
- `angle_avg_2d(lam, q, v=1.0, gam=1.0, n=512)` — 
- `angle_avg_2d_closed(lam, q, v=1.0, gam=1.0)` — 
- `angle_avg_3d(lam, q, v=1.0, gam=1.0, n=2000)` — 
- `angle_avg_3d_closed(lam, q, v=1.0, gam=1.0)` — 
- `sm_residual(n=64, lam=0.7, q=1.3, v=1.0, gam=1.0)` — Ex. 5.7：A = diag(λ+γ+ivq cos θ_k)，|0⟩ = (1,…,1)/√n；(5.75) 乘 (5.73) 與單位矩陣的最大差。
- `rt2_cont(r, t, v=1.0, gam=1.0)` — 二維密度的連續部分（更正）：(γ/2πv) e^{−γt + γs/v}/s，s = √(v²t² − r²)（向量；r ≥ vt 為 0）。
- `rt2_pole(r, t, v=1.0, gam=1.0)` — 
- `rt2_cut(r, t, v=1.0, gam=1.0)` — 
- `rt2_book_pole(r, t, v=1.0, gam=1.0)` — 照書上 (5.84)：2γ√(π/2) e^{−γt} cosh(γ√(v²t²−r²))/(v√(v²t²−r²))（書上的歸一化是 2πρ）。
- `rt2_book_cut(r, t, v=1.0, gam=1.0)` — 照書上 (5.89) 的連續部分：−2π γ e^{−γt} e^{−γs/v}/(v s)。
- `rt2_radial(r, t, v=1.0, gam=1.0)` — 
- `rt2_mass(t, v=1.0, gam=1.0, n=400, book=False)` — ∫ 2πrρ dr ＋ 尖峰（r = vt sin φ）。book=True：用書上 (5.84)+(5.89)（2πρ 的歸一化：∫ r(…) dr）。
- `rt2_hankel(q, t, v=1.0, gam=1.0, n=400)` — 
- `rt2_laplace(q, lam, v=1.0, gam=1.0)` — (5.78)：1/(√((λ+γ)²+q²v²) − γ)。
- `hankel_pair(r, n=3000, qmax=14.0)` — (5.88) 的核對：f(r') = e^{−r'²/2} → F(q) = ∫ r'J0(qr')f dr'（數值）→ ∫ q J0(qr) F(q) dq（數值）；應回到 f(r)。
- `rt3_pole_book(r, t, v=1.0, gam=1.0)` — (5.101)（書上歸一化 2π²ρ）：(γ/v)² ∫_0^{π/2} x³/(r sin²x) e^{−γt} e^{γxt cot x} sin(γxr/v) dx。
- `rt3_5106(r, t, v=1.0, gam=1.0, book=True)` — (5.106)：−e^{−γt}∫du/sin²u (γ/v)² u³ sin(γur/v)/r e^{∓γut cot u}（書上 −；更正 +，即 −(5.101)）。
- `rt3_cut_book(r, t, v=1.0, gam=1.0, n=600)` — (5.103) 第一項（x = tanh y）：e^{−γt}(1/r)(γ/v)² Re∫ sech²y (π/2 − iy)³ exp[−i(r − vt tanh y)(γ/v)(π/2 − iy)] dy。
- `rt3_once_book(r, t, v=1.0, gam=1.0)` — (5.103) 第三項：e^{−γt}(γ/v)(1/(rvt))(π/2) ln((1+r/vt)/(1−r/vt))（只轉一次彎的人；書上歸一化）。
- `rt3_density(r, t, v=1.0, gam=1.0)` — 三維密度的連續部分（標準歸一化）＝ [(5.101) + (5.103)]/(2π²)；r ≥ vt 為 0。
- `rt3_radial(r, t, v=1.0, gam=1.0)` — 徑向機率密度 4πr²ρ（∫_0^{vt} dr ＝ 1 − e^{−γt}）；陣列逐點算。
- `rt3_r2sigma_book(r, t, v=1.0, gam=1.0)` — 書上 Fig. 5.6 的縱軸 r²Σ（Σ = 2π²ρ）＝ (π/2)·4πr²ρ。
- `rt3_once_radial(r, t, v=1.0, gam=1.0)` — 
- `rt3_mass(t, v=1.0, gam=1.0, n=300)` — 
- `rt3_fourier(q, t, v=1.0, gam=1.0, n=300)` — 
- `rt3_laplace(q, lam, v=1.0, gam=1.0)` — (5.91)–(5.93)：g/(1 − γg)，g = arctan(vq/(λ+γ))/(vq)。
- `rt_mc(dim: int, t: float, M: int=100000, seed: int=5, v: float=1.0, gam: float=1.0)` — 固定種子的粒子模擬：一維＝電報過程（以速率 γ 反轉）；二三維＝以速率 γ 換成隨機新方向。回傳一維位置或二三維距離。
- `mc_fraction(vals, a, b, absval=True)` — 
- `radial_fraction(dim: int, a: float, b: float, t: float, v=1.0, gam=1.0, n=400)` — 公式：a ≤ r < b 的比例（連續部分；一維是 |x|）。
- `hazard(t, kind='poisson', p0=P0)` — 
- `survival(t, kind='poisson', p0=P0)` — 
- `survival_numeric(t, kind='poisson', p0=P0, n=400)` — 
- `wait_density(t, kind='poisson', p0=P0)` — 
- `mean_wait(kind='poisson', p0=P0)` — B(0) = ∫a dt = 平均等待時間。
- `AB(lam, kind='poisson', p0=P0)` — 
- `AB_numeric(lam, kind='poisson', p0=P0, T=60.0, n=6000)` — 
- `steps_prob(m: int, t: float, kind='poisson', p0=P0)` — 
- `steps_prob_laplace(m: int, t: float, kind='poisson', p0=P0)` — 
- `steps_prob_conv(m: int, t: float, kind='poisson', p0=P0, n=3000)` — (5.110)(5.111)：f_0 = a，f_{k+1}(τ) = ∫_0^τ a(s)q(s) f_k(τ−s) ds（梯形法數值卷積）。
- `mean_steps(t, kind='poisson', p0=P0)` — 
- `ctrw_lattice(x: int, t: float, kind='poisson', p0=P0, mmax: int=120)` — 
- `ctrw_bessel(x: int, t: float, p0=P0)` — Ex. 5.8：1/(λ + p0(1−cos q)) 的反轉 = e^{−p0t} I_|x|(p0 t)（walk_bc.ine_list）。
- `ctrw_var(t, kind='poisson', p0=P0)` — 
- `ctrw_var_numeric(t, kind='poisson', p0=P0, X=40)` — 
- `ctrw_gauss(x, t, kind='poisson', p0=P0, C=0.5)` — (5.131)–(5.133) ÷ (2π)（一維）：√(B0/4πCt) e^{−B0x²/4Ct}，B0 = B(0)。
- `ctrw_transform(lam, q, kind='poisson', p0=P0)` — 
- `ex58(lam, q, p0=P0)` — 
- `mortal_total(t, p0=1.0, p1=0.1)` — 
- `mortal_series(t, p0=1.0, p1=0.1, mmax=200)` — Σ_m (還活著且走了 m 步) = Σ_m (p0t)^m e^{−(p0+p1)t}/m!（Ex. 5.9，第二種）。
- `mortal_laplace(t, p0=1.0, p1=0.1)` — Talbot 反轉 Σ A^m B = B/(1 − A)，A = p0/(p0+p1+λ)、B = 1/(p0+p1+λ)（第三種）。
- `check()` — 

## `combkit.young_extras`

young_extras.py —— youngart（Fulton《Young Tableaux》全書 60 張 4K 藝術桌布）用到、可以重複使用的計算（Python＋numpy）

常數：`TAU`

- `patience(w)` — 
- `slide_path(T)` — 洞從 (0,0) 出發往右或往下（取較小者，相等取下面），走到外角為止。
- `bump_routes(w)` — 
- `shadow_lines(w)` — 回傳折線串列，每條是 (列, 行) 座標的陣列；邊界用 0.5 與 n+0.5。
- `four_corners(w)` — 
- `rsk_rows(w)` — 
- `spectral_xy(n, edges)` — 
- `pipe_route(D, m, i, cs)` — 
- `ruling(fam, th, t)` — 
- `line_meets_hyperboloid(p0, d)` — 
- `ellipse_curve(a, b, phi, c=(0.0, 0.0))` — 
- `line_curve(p0, ang)` — 
- `cubic_curve(c=(0.0, 0.0), s=1.6, k=0.9)` — 
- `meet(A, B, N=60000)` — 
- `check()` — 

## `combkit.draw.styles`

styles.py —— 全專案共用的視覺規則

常數：`CJK_FONT`、`POINT_COLOR`、`CIRCLE_COLOR`、`TEXT_COLOR`、`MUTED_COLOR`、`TYPE1_COLOR`、`TYPE2_COLOR`、`TYPE3_COLOR`、`TYPE4_COLOR`、`GENERIC_COLOR`、`CROSSING_COLOR`、`HIGHLIGHT_COLOR`、`ROOT_COLOR`、`OK_COLOR`、`D_COLOR`、`REGION_COLORS`、`SIZE_TO_COLOR`、`DOT_RADIUS`、`BLOCK_STROKE`、`BLOCK_FILL_OPACITY`、`CIRCLE_STROKE`

- `block_color(size: int)` — 依群組的點數回傳顏色。
- `zh(text: str, font_size: float=32, color=TEXT_COLOR, **kwargs)` — 中文（或中英混合）文字。
- `mixed(*parts, font_size: float=32, color=TEXT_COLOR, buff: float=0.12, t2c=None)` — 把中文與數學式排成一行。
- `boxed(mob: Mobject, color=YELLOW, buff: float=0.25, stroke_width: float=3)` — 在物件外面加一個框，回傳 VGroup(框, 物件)。

## `combkit.draw.scene_base`

scene_base.py —— 所有場景共用的基底類別

- **class `ChuScene`** — 
  - `.setup()` — 
  - `.show_header(text: str, animate: bool=True)` — 
  - `.make_caption(*parts, color=WHITE, t2c=None)` — 
  - `.say(*parts, wait: float=0.0, run_time: float=0.6, color=WHITE, t2c=None)` — 在底部顯示字幕；已有字幕時會換掉舊的。
  - `.clear_caption(run_time: float=0.4)` — 
  - `.clear_stage(run_time: float=0.7)` — 淡出標題以外的所有東西（字幕也會淡出）。
  - `.fade_all(run_time: float=0.8)` — 把畫面上所有東西淡出。

## `combkit.draw.circle_points`

circle_points.py —— 圓、圓周上的點、弦

- `point_angle(i: int, n: int, start_angle: float=PI / 2, clockwise: bool=True)` — 第 i 個點（共 n 個）在圓上的角度。
- `circle_coords(n: int, radius: float=2.5, center=ORIGIN, start_angle: float=PI / 2, clockwise: bool=True)` — 回傳 n 個等距圓周點的座標，shape = (n, 3)。
- `make_circle_points(n: int, radius: float=2.5, center=ORIGIN, dot_radius: float=DOT_RADIUS, circle_color=CIRCLE_COLOR, circle_stroke: float=CIRCLE_STROKE, point_color=POINT_COLOR)` — 建立一個圓與圓周上的 n 個點。
- `chord_between(coords, i: int, j: int, color=TYPE1_COLOR, stroke_width: float=BLOCK_STROKE)` — 連接第 i 與第 j 個點的弦。
- `label_points(coords, center=ORIGIN, labels=None, buff: float=0.32, font_size: float=26, color=GREY_A, start: int=1)` — 在每個點的外側放上標籤（預設 1, 2, 3, ...）。
- `chords_intersection(coords, chord_a, chord_b)` — 兩條弦（各為一對點編號）在圓內的交點；不相交則回傳 None。

## `combkit.draw.polygon_groups`

polygon_groups.py —— 群組（block）的畫法，以及「一整個分法」的組裝

- `make_block(indices, coords, color=None, fill_opacity: float=BLOCK_FILL_OPACITY, stroke_width: float=BLOCK_STROKE)` — 依點編號畫出一個群組。
- `make_blocks(blocks, coords, **kwargs)` — 一次畫出多個群組。blocks 是「點編號串列」的串列。
- **class `Configuration`** — 一個完整的分法：圓 + 群組 + 圓周點。
- `make_configuration(n_points: int, blocks, **kwargs)` — 建立一個分法（圓 + 群組 + 點）。blocks 例如 [(0, 1, 2), (3, 4)]。
- `dot_color_anims(dots, indices, color)` — 把指定的點變成某個顏色（回傳動畫串列，放進 self.play）。
- `mark_crossing(block: VMobject)` — 把一個群組塗成「交叉」的紅色。
- `highlight_block(block: VMobject, scale_factor: float=1.08, color=None, **kwargs)` — 讓一個群組閃一下（放大再縮回）。回傳動畫。
- `circular_segment(i: int, j: int, n: int, radius: float, center=ORIGIN, color=TEAL, opacity: float=0.22, n_samples: int=48)` — 弓形區域：從第 i 個點沿圓周順時針走到第 j 個點的圓弧，加上弦 i-j 圍成的區域。
- `root_regions(root_block, n: int, radius: float, center=ORIGIN, colors=None, opacity: float=0.22)` — 根群組（依順時針排序的點編號）切出的各個區域。

## `combkit.draw.root_figure`

root_figure.py —— Scene 05 之後反覆用到的「根分解」示意圖，以及群組小圖示

常數：`ROOT_N_POINTS`、`ROOT_BLOCK`、`ROOT_SUB_BLOCKS`、`ROOT_ALT_R2`、`ROOT_R_VECTORS`、`ROOT_N_VEC`、`ROOT_V_VEC`

- **class `RootFigure`** — 屬性：
- `block_icon(size: int, radius: float=0.28, color=None, stroke_width: float=4, dot_radius: float=None)` — 不放在圓上的「群組圖示」：正 size 邊形（size=2 時是一條線段）加上頂點。

## `combkit.draw.lattice_paths`

lattice_paths.py —— Raney 公式證明用的「路徑」圖形

常數：`UP_STEP_COLOR`

- `letter_color(letter: int, v_vec)` — 
- `weight_tex(w: int)` — 
- `word_tiles(word, v_vec, size: float=0.5, buff: float=0.08, font_size: float=26)` — 把一串步畫成一排方塊，方塊上寫高度變化（+1、-1、-2）。
- **class `PathMob`** — 一條路徑。origin 是「第 0 個點」的位置，start_height 是它的高度（只影響 .heights）。
  - `.vertex(j: int)` — 
  - `.y_of(h: float)` — 
  - `.x_of(j: float)` — 
- `level_line(path: PathMob, h: float, j0: float, j1: float, color=GREY_B, stroke_width: float=2, dashed: bool=True)` — 在路徑的座標系裡，從第 j0 格到第 j1 格、高度 h 的水平線。
- `path_thumbnail(word, v_vec, unit: float=0.17, origin=ORIGIN, mark_bad: bool=True, stroke_width: float=3.2, dot_radius: float=0.03)` — 小圖：路徑 + 高度 0 的虛線；壞路徑在第一個高度 <= 0 的點畫紅點。

## `combkit.draw.catalan_draw`

catalan_draw.py —— 卡塔蘭基礎模組（場景 00a–00m）的圖形

常數：`CATALAN_COLOR`、`LEFT_COLOR`、`RIGHT_COLOR`、`EDGE_COLOR`、`LEAF_COLOR`、`PART_COLORS`、`DYCK_STEPS`、`LATTICE_STEPS`、`FAMILIES`、`FAMILY_NAMES`

- `part_color(part, default=CATALAN_COLOR)` — 
- `polygon_vertices(m: int, radius: float=1.0, center=ORIGIN)` — 正 m 邊形的頂點。底邊（第 0 與第 m-1 點）水平放在最下面；
- **class `TriangulationMob`** — 凸 m 邊形 + 對角線 + 頂點。子物件順序：[fills, outline, diagonals, dots]。
  - `.vertex(k: int)` — 
  - `.region(indices, color, opacity: float=0.35)` — 
  - `.edge(a: int, b: int, color=ROOT_COLOR, width: float=7)` — 
- `triangulation_split(tree, m: int)` — 回傳 (k, left_indices, right_indices)：含底邊的三角形是 (0, k, m-1)。
- `triangulation_item(tree, radius: float=1.0, center=ORIGIN, colored: bool=False, **kwargs)` — 由樹畫出凸 (n+2) 邊形的三角剖分；colored=True 時把根三角形、左右子多邊形上色。
- **class `BinaryTreeMob`** — 滿二元樹：內部節點（實心點）、葉子（小方塊或字母）、邊。
  - `.color_parts()` — 
- `tree_item(tree, scale: float=1.0, colored: bool=False, **kwargs)` — 
- `paren_mob(tree, factors='abcdefghij', font_size: float=40, colored: bool=False, color=WHITE)` — 回傳 Tex，附加屬性：
- `catalan_product_tex(k: int, j: int, tail: str='', font_size: float=34)` — 「C_k C_j」＋ tail（例如 "=2"），C_k 用左邊的青色、C_j 用右邊的粉紅。
- `paren_item(tree, scale: float=1.0, colored: bool=False, factors='abcd')` — 
- `paren_depth_stages(tokens)` — 把括號依「由內而外」分批：回傳 [[(開括號位置, 閉括號位置), ...], ...]，
- `non_paren_glyphs(mob: Tex)` — 括號以外的符號（因子與乘號），依閱讀順序。
- `paren_build_animations(plain: Tex, target: Tex, run_time: float=0.7)` — 從沒有括號的「a·b·c·d」長出某一種括號法：
- **class `StepPathMob`** — 一條由「步」組成的路徑。子物件：[steps, dots]。
  - `.vertex(j: int)` — 
- `dyck_path_mob(word: str, unit: float=0.5, origin=ORIGIN, **kwargs)` — 
- `lattice_path_mob(word: str, unit: float=0.5, origin=ORIGIN, **kwargs)` — 
- `lattice_grid(nx: int, ny: int, unit: float=0.5, origin=ORIGIN, color=GREY_D, stroke_width: float=1.5)` — 格子線：x = 0..nx、y = 0..ny。
- `lattice_diagonal(n: int, unit: float=0.5, origin=ORIGIN, color=GREY_A, stroke_width: float=2, offset: float=0.0)` — 對角線 y = x + offset（從 x=0 畫到 y=n 或 x=n）。
- `lattice_thumbnail(word: str, n: int, unit: float=0.25, origin=ORIGIN, valid=None, stroke_width: float=3.5)` — 小圖：格子 + 對角線 + 路徑。valid=None 時路徑用一般藍色。屬性 .path、.grid、.diag。
- `dyck_item(tree, scale: float=1.0, colored: bool=False, **kwargs)` — 
- `lattice_item(tree, scale: float=1.0, colored: bool=False)` — 
- `matching_item(tree, radius: float=0.4, center=ORIGIN, colored: bool=False, **kwargs)` — 由樹畫出圓周 2n 點的不交叉完美配對；colored=True 時含第 0 點的弦是根。
- `catalan_item(family: str, tree, scale: float=1.0, colored: bool=False)` — 畫出同一棵樹在某一種物件下的樣子，中心放在原點（大小約 0.9 × scale）。

## `combkit.draw.mit_draw`

mit_draw.py —— MIT 18.212 L1–4 補充場景（M01–M14）的共用圖形

常數：`NEST_COLOR`、`HOOK_COLOR`、`CORNER_COLOR`、`CELL_STROKE`

- `arc_curve(p, q, color=CATALAN_COLOR, stroke_width: float=4, height_ratio: float=1.0, below: bool=False)` — p、q（同一條水平線上）之間的半圓弧，從 p 畫到 q（ShowCreation 的方向）。
- **class `ArcDiagramMob`** — 直線上 n_points 個點（0 起算，由左到右，間距 gap，整體中心 center）與上方的弧。
  - `.pairs()` — 目前每條弧的端點（由 .arcs 算出，所以 .copy() 之後仍正確）。
  - `.point(i: int)` — 
  - `.make_arc(i: int, j: int, color=None, stroke_width=None, height_ratio=None)` — 
  - `.add_arc(i: int, j: int, color=None)` — 
  - `.color_arcs(colors)` — 
- **class `YoungDiagramMob`** — Young 圖（英式）：格子 (i, j) = 第 i 列、第 j 行（0 起算），第 0 列在最上面、每列靠左。
  - `.cell(i: int, j: int)` — 
  - `.cell_center(i: int, j: int)` — 
  - `.cell_group(cells)` — 
  - `.fill_cells(cells, color, opacity: float=0.5)` — 
  - `.set_entries(table, font_size=None, color=WHITE)` — table[i][j] 是數字、字串或 Mobject；None 表示空白。會換掉原本的 entries。
  - `.entry(i: int, j: int)` — 
  - `.hook_highlight(i: int, j: int, color=HOOK_COLOR, opacity: float=0.35, cell_color=None, cell_opacity: float=0.7)` — 鉤（自己、同列右邊、同行下面）的覆蓋方塊；第 0 個是格子本身（顏色 cell_color，預設同 color）。
- **class `PlaneTreeMob`** — 平面樹（counting 第 6 節的 tuple 形式）。根在 root_pos、往下長；同一個頂點的小孩由左到右。
  - `.labels()` — {頂點: 文字 Mobject}（由 .label_group 算出，所以 .copy() 之後仍正確）。
  - `.vertex(v: int)` — 
  - `.edge_to(v: int)` — 
- `box_row(values, size: float=0.55, buff: float=0.08, direction=RIGHT, colors=None, color=GREY_A, fill_opacity: float=0.15, font_size=None, text_color=WHITE, stroke_width: float=2.5, center=ORIGIN)` — 一排方塊，每個方塊 = VGroup(方框, 文字)，有 .value、.square、.label。
- `perm_row(perm, highlight=(), highlight_color=YELLOW, **kwargs)` — 排列（或任何數列）的方塊列；highlight 裡的位置（0 起算）用 highlight_color 標出。

## `combkit.draw.fulton_draw`

fulton_draw.py —— Fulton《Young Tableaux》第 2–4 章場景（F09–F19）共用的畫圖工具（ManimGL）

常數：`HOLE_COLOR`、`MU_COLOR`、`RULE_COLORS`、`RULE_TEX`

- `assert_in_frame(*mobs)` — 物件都在畫面的安全範圍內（不碰到標題與字幕）。
- `assert_apart(a, b, gap=0.05)` — 兩個物件的外框不重疊（至少一個方向分開 gap）。
- `tuple_str(lam)` — 
- `labeled(tex: str, mob: Mobject, font_size=36, buff=0.25)` — VGroup(公式, 物件)，左右排好；公式的 Tex 用 tex。
- `tableau_mob(T, cs: float=0.55, fs=None, top_left=None, center=ORIGIN, **kwargs)` — 填好數字的 Young 圖；每個數字物件記下 .tab_value。空表傳回只有一個點的 VGroup 不適用，請自行處理。
- **class `SkewMob`** — 斜表或有洞的表：S 的每一列可以有 None。列首連續的 None 是 μ 的格子（灰色，show_mu=False 時不畫），
  - `.cell(i, j)` — 
  - `.entry(i, j)` — 
  - `.cell_center(i, j)` — 
- `word_tex(w, colors=None, font_size=36, gap='\\ ')` — 字的文字版（例如 5\,6\,4…）；colors = {位置: 顏色}。字母依序是 tex[k]（每個字母都是一位數時）。
- **class `WordMob`** — 一排字母方塊（同 mit_draw.box_row 的樣式）。位置（slot）固定，方塊會移動：
  - `.letters()` — 
  - `.slot(k)` — 
  - `.bracket(i, n=3, label=None, color=YELLOW, font_size=26)` — 
  - `.swap_anims(a, b, run_time=0.6)` — 
  - `.apply_move(scene, i, rule, run_time=0.6, hold=0.25, show_label=True)` — 一步基本 Knuth 變換：K'／K'⁻¹ 交換三個字母中的後兩個，K''／K''⁻¹ 交換前兩個。回傳新字。
  - `.color_slots(colors: dict, opacity=0.3)` — 
- **class `TwoRowArray`** — 雙列陣列（上列 top、下列 bottom），左右有圓括號。.col(k) = VGroup(上, 下)，.cols。
  - `.col(k)` — 
- **class `BallMatrixMob`** — 矩陣球圖：m×n 的格子，(i, j) 裡畫 a(i,j) 個球（由西北往東南斜排），numbering = tableaux.ball_numbering(A)
  - `.ball(i, j, t=0)` — 
  - `.cell_square(i, j)` — 
  - `.cell_center(i, j)` — 

## `combkit.draw.rep_draw`

rep_draw.py —— Fulton 第 7 章場景（F28–F41）共用的畫圖元件（ManimGL）

常數：`ROW_COLOR`、`COL_COLOR`、`OP_COLOR`、`POS_COLOR`、`NEG_COLOR`

- `perm_tex(p)` — 
- `oneline_tex(p)` — 
- **class `TabloidMob`** — 列 tabloid {T}：數字排成列，每列上下各一條橫線（Fulton p.85 的畫法）。
  - `.entry(i, j)` — 
- **class `ColumnTabloidMob`** — 帶方向的行 tabloid [T]（Fulton p.95）：每行左右各一條直線；sign = ±1 畫在前面。
- `signed_sum(terms, font_size=40, buff=0.2, dot=False)` — terms = [(係數, Mobject), ...]；回傳 VGroup，.terms[k] = VGroup(係數記號, 物件)。
- **class `MatrixMob`** — 整數矩陣（rows）；括號依高度拉長。.entry(i, j)、.col(j)、.row(i)。
  - `.entry(i, j)` — 
  - `.col(j)` — 
  - `.row(i)` — 

## `combkit.draw.formulas`

formulas.py —— 常用公式的 Tex 物件

常數：`GF_T2C`

- `catalan_formula(font_size: float=48)` — C_n = 1/(n+1) * binom(2n, n)
- `vec_str(values)` — 
- `create_vector_formula(n_values, v_values, font_size: float=44)` — 例如 \vec n=(1,1), \vec v=(2,3)。
- `create_total_points_formula(n_values, v_values, font_size: float=44)` — N = n_1 x v_1 + ... = 總點數（帶入數字）。
- `create_multinomial_formula(n_values, v_values, font_size: float=44, steps: bool=True)` — 把數字代入論文式 (2)：
- `general_formula(font_size: float=44)` — 論文式 (2) 的一般形式。
- `gf_definition(font_size: float=44)` — 論文式 (3)：生成函數 D 的定義。
- `gf_equation(font_size: float=48)` — 論文式 (6)：D = 1 + x_1 D^{v_1} + ... + x_k D^{v_k}。
- `recurrence_formula(font_size: float=40)` — 論文式 (5)：遞迴關係。
- `size_colored(text: str, size: int, font_size: float=40)` — 用「size 點群組」的顏色畫一段 Tex。

## `combkit.draw.art_kit`

藝術桌布共用工具（walkart；combkit 的 combkit/draw/art_kit.py 是同一份）：

常數：`BG_DIR`、`PAL`、`ROMAN`、`SERIF`、`SANS`、`STEPS4`、`D4`、`GROW`、`CAP`

- `rng(seed)` — 
- `hexrgb(h)` — 
- `lerp_hex(a, b, t)` — 
- `ramp(colors, n)` — n 個沿 colors 漸變的顏色。
- `background(pal, center=(0.5, 0.5), spread=1.0, grain=0.005, seed=0, w=3840, h=2160)` — 徑向漸層＋暗角＋細微顆粒的 4K PNG（快取在 _bg/），回傳鋪滿畫面的 ImageMobject。
- `separator()` — 
- `to3(pts)` — 
- `path(pts, color=WHITE, width=2.0, opacity=1.0, smooth=False)` — 
- `gradient_path(pts, colors, width=2.0, opacity=1.0, smooth=False)` — 
- `glow(mob, color=None, width=None, layers=((7.0, 0.05), (3.5, 0.12), (1.0, 1.0)))` — 同一個形狀疊三層：粗而淡 → 細而亮。
- `dots(pts, color=WHITE, radius=0.03, opacity=1.0, glow_factor=0.0)` — 
- `glow_dots(pts, color=WHITE, radius=0.12, opacity=1.0, glow_factor=2.0)` — 
- `fit(pts, width=None, height=None, center=ORIGIN)` — 把點縮放進 width×height 的框並置中（保持比例）。
- `contours(F, xs, ys, levels)` — F[iy, ix] 在格點 (xs, ys) 上的等高線：回傳 [(level, (n,2) 陣列), ...]。
- `grid(f, xr, yr, n=400)` — 
- `lattice_walk(n, seed)` — 
- `gauss_walk(n, seed, dim=2)` — 
- `walk1d(n, seed, k=1)` — 
- `pivot_saw(n, seed, attempts=None)` — 方格上的自迴避漫步（樞軸法，從直棒開始）。
- **class `Plate`** — 
  - `.P()` — 
  - `.C(i)` — 
  - `.add(*mobs)` — 這台 wgpu 渲染器的毛病：DotCloud 緊接在 VMobject 後面畫會出錯，中間夾一張透明小圖就好。
  - `.on_add(m)` — 每加入一個物件（不含夾在中間的透明小圖）就呼叫一次；動畫版用它記下加入順序。
  - `.construct()` — 
  - `.art()` — 
  - `.caption()` — 
- `render_png(cls, res, out)` — 建構 Plate（或任何 Scene）並把最後畫面存成 RGB PNG。
- `contact_sheet(files, out, cols=4, tw=480, th=270)` — 把多張圖縮小拼成一張，左上角標檔名前兩個字（編號）。
- `win(t0, t1, T, f=smooth)` — 把整段動畫的 alpha（0–1，對應 T 秒）換成只在 [t0, t1] 秒之間走完的 rate function。
- `reveal(m, rf, T)` — 
- `grow_points(m, keys, op, rf, T, w=0.03)` — 點雲依 keys（0–1）一顆顆亮起。
- `along(P, u)` — 
- `travelers(scene, polys, color, T, t0=0.0, dur=None, stagger=0.0, rf=None, radius=0.14, glow_factor=2.5)` — 幾顆光點各沿一條折線走（依點的編號均勻走，和 ShowCreation 同步）；還沒出發或已到的點隱藏。
- `keep(m, fn, T)` — 記下 m 一家每個物件的點，之後每格用 fn(原點, alpha) 重設。
- `spin(items, angle=None, T=T_B)` — 
- `mirror(items, axis, T=T_B)` — 每個 (物件, 中心) 沿平面上的 axis 方向翻 180° 再翻回來（3D 翻轉＝鏡射；youngart 20、30）。
- `twinkle(m, seed=0, amp=0.75, T=T_B)` — 點雲逐點閃爍；頭尾不變，接在「長出來」後面不會跳。
- `pulse(m, seed=0, amp=0.45, T=T_B)` — 點雲逐點呼吸（半徑）。
- `keyed_reveal(m, keys, w=0.03)` — 把 m 的「長出來」換成依 key（0–1）一個個出現：keys 是陣列時對 m 的子物件，純量時對 m 本身。
- `animate_plate(base, cfg=None)` — 把 Plate 子類別包成動畫版（最後一格＝桌布）。cfg：together、A(scene)、B(scene)、drift、
- `render_video(cls, res, fps, outdir, stem, crf=16)` — 一個行程裡直接渲染成 mp4（不用 manimgl 指令）；呼叫前同樣要讓 manimlib 以「寫檔」初始化。

## `combkit.draw.lattice_draw`

lattice_draw.py —— 《Elements of the Random Walk》第 2 章（格子上的生成函數）的畫圖元件（ManimGL）。數學在 walk_gf.py（先算、再畫）。

常數：`SITE_COLOR`、`NEIGHBOR_COLOR`、`TRAIL_COLOR`、`POS_COLOR`、`NEG_COLOR`、`ZERO_COLOR`、`GRID_COLOR`、`START_COLOR`、`END_COLOR`、`DEPTH`

- `lat_point(x: float, y: float, cell: float=0.5, origin=ORIGIN)` — 
- `lattice_lines(n: int, cell: float=0.5, origin=ORIGIN, color=GRID_COLOR, stroke_width: float=1.5)` — 
- `lattice_dots(n: int, cell: float=0.5, origin=ORIGIN, radius: float=0.035, color=GREY_C)` — 
- `count_grid(values: dict, cell: float=0.5, origin=ORIGIN, font_size: float=22, color=WHITE, zero_dots: bool=False, n: int | None=None)` — values：{(x, y): 整數}。數字的自訂屬性 .lat_xy = (x, y)。zero_dots=True 時，−n..n 裡值為 0 的格點畫小灰點。
- `grid_number(g: VGroup, x: int, y: int)` — 
- `neighbor_star(cell: float=0.9, origin=ORIGIN, labels: bool=True, font_size: float=26)` — 
- `trail(points, cell: float=0.2, origin=ORIGIN, color=TRAIL_COLOR, stroke_width: float=3, ends: bool=True, end_radius: float=0.07)` — 
- `heat_color(v: float, vmax: float=1.0)` — 
- `heatmap(fn, x_range=(-PI, PI), y_range=(-PI, PI), n: int=30, width: float=4.0, height: float | None=None, origin=ORIGIN, vmax: float=1.0, opacity: float=1.0)` — 色塊 (i, j)：i 是列（由上往下）、j 是行（由左往右）；中心值 fn(x_j, y_i)。
- `heat_cells(values: dict, cell: float=0.2, origin=ORIGIN, vmax: float=1.0)` — 
- `bz_frame(width: float=4.0, origin=ORIGIN, labels: bool=True, color=GREY_B, font_size: float=24)` — 
- `roots_of_unity(n: int, radius: float=1.2, origin=ORIGIN, m: int=1, color=TEAL)` — 
- `root_chain(n: int, m: int=1, length: float=0.8, origin=ORIGIN, color=TEAL)` — 
- `iso_point(x: float, y: float, z: float, unit: float=1.0, origin=ORIGIN)` — 
- `lattice3d(kind: str='sc', unit: float=1.0, origin=ORIGIN, color=TEAL, cube_color=GREY_D)` — 
- `branch_cut(axes: Axes, x0: float, x1: float, color=RED, amp: float=0.06, teeth: int=14, stroke_width: float=3)` — 

## `combkit.draw.osc_draw`

osc_draw.py —— 對稱群與耦合振盪子的畫圖元件（ManimGL）。數學在 symosc.py（先算、再畫）。

常數：`C_I`、`C_J`、`C_K`、`MASS_COLOR`、`SPRING_COLOR`

- `tri_points(radius: float=1.5, center=ORIGIN)` — 
- `spring(p, q, coils: int=7, amp: float=0.08, color=SPRING_COLOR, stroke_width: float=2.5)` — 
- `osc_triangle(radius: float=1.5, center=ORIGIN, labels: bool=True, font_size: float=26)` — 
- `masses_updater(tracker_fn, radius: float=1.5, center=ORIGIN, amp: float=0.45, color=MASS_COLOR)` — 
- `radial_bars(values, radius: float=1.5, center=ORIGIN, scale: float=0.35, color=MASS_COLOR, width: float=7)` — 
- `phasor(z: complex, at, scale: float=0.45, color=YELLOW, zero_dot: bool=True)` — 
- `phasors(zs, points, scale: float=0.45, color=YELLOW, zero_dot: bool=True)` — 
- `level_diagram(levels, x: float=0.0, y0: float=-2.0, yscale: float=1.0, width: float=1.2, color=WHITE, font_size: float=26, label_side=RIGHT)` — 
- `q8_positions(r_in: float=0.9, r_out: float=1.9, center=ORIGIN)` — 
- `q8_graph(pos: dict, mul, font_size: float=26, node_radius: float=0.22)` — mul(a, b) 是 Q8 的乘法（symosc.q8_mul）。

## `combkit.draw.raney_draw`

raney_draw.py —— 整數序列的 Raney 引理（動畫 R01–R06）用到的圖形

常數：`UP_COLOR`、`BIG_UP_COLOR`、`WIN_COLOR`、`USELESS_COLOR`、`LEAF_COLOR`、`UNARY_COLOR`

- `step_color(x: int)` — 
- **class `IntPathMob`** — 一條路徑。origin 是第 0 點的位置（高度 start_height）。
  - `.vertex(j: int)` — 第 j 點目前的座標（整條路徑移動或縮放後也正確）。
  - `.get_unit()` — 
  - `.y_of(h: float)` — 
- `level_line(path: IntPathMob, h: float, j0: float, j1: float, color=GREY_C, stroke_width: float=2, dashed: bool=True)` — 
- `step_tiles(seq, size: float=0.55, buff: float=0.08, font_size: float=30)` — 一排方塊，每格寫一個步的值（Tex），邊框依值上色。
- **class `CardMob`** — 摸牌遊戲的牌。w ≥ 1：特殊牌（寫 +w，顏色同 w 點群組）；w = 0：無用牌；w = None：勝利牌（星號）。
- `draws_badge(k: int, radius: float=0.24, font_size: float=26)` — 「還能摸 k 張」的小圓牌：k ≥ 1 綠、k = 0 紅。
- `symbol_color(j: int)` — 
- `word_tiles(word, size: float=0.5, buff: float=0.06, font_size: float=28, sep_after=())` — 字（或字串）的方塊列；sep_after：在這些位置之後多空一點（用來把字列切開）。
- **class `WordTreeMob`** — 把字（前序記法）畫成平面樹：葉子依序等距，內部節點在第一個與最後一個小孩的中間；節點依元數上色。
  - `.node(k: int)` — 
- `forest_mob(words, gap: float=0.6, **kwargs)` — 字列（一串字）畫成森林，由左到右。
- `rank_walk(word, unit: float=0.35, origin=ORIGIN, **kwargs)` — Raney 的秩路徑：由左往右讀，符號 j 走 j − 1；線段依符號上色。
- `chu_walk(word, unit: float=0.35, origin=ORIGIN, **kwargs)` — Chu 講義的路徑：把字串倒過來讀，符號 j 走 1 − j；線段依符號上色。
- `caption_zh(text: str, font_size: float=26, color=GREY_A)` — 
- `readings_row(seq, starts, center=ORIGIN, spacing: float=2.5, unit: float=0.22, label_gap: float=0.25, mark: str='good', font_size: float=26)` — 把「從第 i 點開始讀」的路徑排成一列（i 取自 starts），各條路徑的第 0 點在同一高度，基線（高度 0）畫虛線。

## `combkit.draw.schubert_draw`

schubert_draw.py —— Fulton 第 9–10 章（場景 F46–F55）的圖形（ManimGL）

常數：`STAR_COLOR`、`ONE_COLOR`、`ZERO_COLOR`、`BOX_COLOR`

- **class `CellMatrixMob`** — 
  - `.cell(i, j)` — 
  - `.stars()` — 
  - `.ones()` — 
- **class `BoxDiagramMob`** — 
- **class `HasseMob`** — 
- **class `PipeDreamMob`** — 格子 (i, j)（1 起算，i 列 j 行，i+j ≤ m）。十字畫兩條線，其餘畫一對彎管（左上、右下的四分之一圓）。
- `poly_tex(f, **kwargs)` — 
- `word_tex(w)` — 

## `combkit.draw.sym_draw`

sym_draw.py —— Harter PSDS 第 1–8 章共用的畫圖元件（ManimGL，`from manimlib import *`）。數學在 sym_mat.py（先算、再畫）。

常數：`KET_COLOR`、`BRA_COLOR`、`P1_COLOR`、`P2_COLOR`、`BEAM_COLOR`、`BOX_COLOR`、`MASS_COLOR`、`SPRING_COLOR`、`BOB_COLOR`、`ROD_COLOR`、`WEAK_COLOR`、`STRONG_COLOR`、`PHASOR_COLOR`、`MIRROR_COLOR`、`N_COLOR`、`H_COLOR`、`LABEL_COLORS`、`SIX_ORDER`、`SIX_COORD_DEG`、`SIX_COORD_MASS`、`SIX_COORD_TEX`、`DIAG_COLORS`、`POS_COLOR`、`NEG_COLOR`、`IPOS_COLOR`、`INEG_COLOR`、`X_COLOR`、`F_COLOR`、`OCT_COLORS`、`OCT_SITES`、`OCT_DIAGONALS`、`XF6_COORDS`、`BODY_COLORS`、`LAB_COLOR`、`CRANK_COLOR`、`PROB_COLOR`、`A_COLOR`、`B_COLOR`、`STRESS_COLOR`、`STRAIN_COLOR`、`EFIELD_COLOR`、`BFIELD_COLOR`、`STIM_COLOR`、`RESP_COLOR`、`P_COLOR`、`E_COLOR`、`J1_COLOR`、`J2_COLOR`、`J3_COLOR`、`POL_COLOR`、`STOKES_COLOR`、`POT_G_COLOR`、`POT_U_COLOR`、`TRAJ_COLOR`、`FRONT_COLOR`、`ENV_COLOR`、`PACKET_COLOR`、`PACKET_ENV_COLOR`、`SPIN_COLOR`、`DRESSED_COLOR`、`BARE_COLOR`

- `entry_tex(x)` — 
- `mat_tex(rows, colors=None, font_size: float=34, h_buff: float=0.45, v_buff: float=0.28, paren_color=WHITE)` — 
- `mat_entry(m: VGroup, i: int, j: int)` — 
- `mat_row(m: VGroup, i: int)` — 
- `mat_col(m: VGroup, j: int)` — 
- `labeled_mat(lhs_tex: str, rows, font_size: float=34, **kw)` — 
- `analyzer(angle: float=0, width: float=1.6, height: float=1.0, label: str | None=None, recombine: bool=False, color=BOX_COLOR, font_size: float=26)` — 
- `port(an: VGroup, side: str, k: int=0)` — 
- `beam(start, end, intensity: float=1.0, color=BEAM_COLOR, max_width: float=14)` — 
- `split_beams(an: VGroup, p1: float, length: float=1.4, spread: float=0.35, **kw)` — 
- `vec_plane(unit: float=1.0, xr: float=3, yr: float=3)` — 
- `vec_arrow(plane, v, color=KET_COLOR, label_tex: str | None=None, font_size: float=30)` — 
- `eigen_line(plane, v, color=GREY_B, extent: float=3.0)` — 
- `circle_image(plane, M, n: int=72, alpha: float=1.0, color=WHITE, stroke_width: float=3)` — 
- `contour_around(axes, center: float, r: float, color=TEAL, stroke_width: float=3)` — 
- **class `SpringChain`** — 
  - `.wall_pts()` — 兩面牆的中點（隨整條鏈移動）。
  - `.eq()` — 各質量的平衡位置（隨整條鏈移動）。
  - `.set_displacements(xs)` — 
- `spring(a, b, coils: int=6, amp: float=0.1, color=SPRING_COLOR, width: float=2.5)` — a → b 的一根彈簧（任意方向）；要動就在 updater 裡 .set_points_as_corners(spring_points(a, b, …))。
- `spring_points(a, b, coils: int=6, amp: float=0.1)` — 
- **class `PendulumRow`** — 一排單擺（側視）。第 i 擺的偏角 = unit·x_i（弧度）；相鄰擺桿在 spring_height（擺長的比例）處以彈簧相連。
  - `.pivots()` — 各支點座標（隨整排移動）。
  - `.bob_position(i: int, x: float)` — 
  - `.set_displacements(xs)` — 
- **class `PendulumRing`** — 俯視的單擺環。第 p 個錘在角度 start_angle + s·2πp/n（clockwise=True 時 s = −1）；位移沿切向或徑向。
  - `.rest()` — 各錘的平衡位置（隨整個環移動）。
  - `.set_displacements(xs)` — 
- **class `PhasorClock`** — 相量鐘面：指針指到 center + radius·(Re z, Im z)。e^{−iωt} 的相量順時針轉。
  - `.set_value(z)` — 
- `band_dots(axes, ks, values, color=YELLOW, radius: float=0.05)` — 能帶／色散的點 (k, 值)（≤ 60 個）。
- `band_curve(axes, f, k0: float=-PI, k1: float=PI, n: int=60, color=YELLOW, stroke_width: float=3)` — y = f(k) 的曲線（n 點，一條 VMobject）。
- `level_lines(axes, values, x0: float, x1: float, color=YELLOW, stroke_width: float=3)` — 能階的水平短線（y = 各值，x 從 x0 到 x1）。
- `quad_contours(plane, F, levels=(0.5, 1.0, 2.0), n: int=60, color=TEAL, stroke_width: float=2.5)` — ½xᵀ𝓕x = level 的等位橢圓（𝓕 對稱正定 2×2，可給 sympy Matrix）。
- `mirror_lines(radius: float=2.0, angles=(120, 60, 0), labels=('\\sigma_1', '\\sigma_2', '\\sigma_3'), color=MIRROR_COLOR, font_size: float=26, center=ORIGIN)` — 固定在紙上的鏡線（Fig. 3.1.2；角度＝ sym_nonabel.MIRROR_DEG）：過 center 的虛線＋標籤（放在右半邊那一端）。
- **class `BobTriangle`** — 俯視的三個擺錘（Fig. 3.1.2）：名字 a、b、c，平衡位置在 0°、120°、240°（虛線小圓）；
  - `.point(p)` — state 座標 → 畫面座標（隨整組移動）。
  - `.set_state(state)` — 
  - `.anim_to(state)` — 移到新的 state 的動畫列表（直線移動；標籤字母不翻轉）。
- **class `SpringTriangle`** — Fig. 3.4.1 的六自由度彈簧質量模型（俯視）：三個質點在 0°、120°、240°（半徑 radius），三邊彈簧 k₁，
  - `.rest()` — 
  - `.disp(x6)` — 6 個正則座標 → 三個質點的位移向量（已乘 unit）。
  - `.set_displacements(x6)` — 
  - `.coord_arrows(length: float=0.55, font_size: float=18, color=KET_COLOR)` — 六個座標箭頭 |g⟩ 與標籤 ⟨g|x⟩（Fig. 3.4.1）。
  - `.displacement_arrows(x6, scale: float=1.0, color=KET_COLOR)` — 每個質點的位移箭頭（長度 = scale·unit·|位移|；太短的略過）。
- `sym_table(col_labels, row_labels, rows, cell_w: float=0.62, row_h: float=0.42, font_size: float=24, corner: str='', first_w: float | None=None, label_color=YELLOW)` — 群乘法表、類代數表、特徵標表、相關表：表頭列＝ col_labels、第一行＝ row_labels（TeX），格子是 TeX 字串／int／sympy（entry_tex）。
- `perm_grid(M, cell: float=0.3, on_color=YELLOW, off_color=GREY_E, row_labels=None, col_labels=None, font_size: float=18, show_values: bool=False)` — 矩陣的格子圖（正則表示 𝓡(g)、塊對角的 T†𝓡T）：非零格填 on_color。→ VGroup(cells, labels, values)，
- **class `View3D`** — 2D 畫面上的 3D 正交投影（不用 3D 相機：這台 wgpu 渲染器的 Surface／DotCloud 緊接在 VMobject 後面畫會出錯）。
  - `.matrix()` — 
  - `.project(v)` — 
  - `.depth(v)` — 
- **class `SphereView`** — Hamilton 轉的球（Fig. 3.1.5–3.1.7、3.2.1）：經緯線（前面亮、後面淡）＋可加的大圓弧、轉軸、點，全部用 View3D 投影。
  - `.project(v)` — 
  - `.add_arc(N1, N2, color=YELLOW, width: float=6, tip: bool=True)` — 
  - `.add_axis(axis, color=GREY_B, length: float=1.3)` — 
  - `.add_point(v, color=WHITE, radius: float=0.06)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
- `polar_curve(f, radius: float=1.5, amp: float=0.3, center=ORIGIN, n: int=60, color=YELLOW, width: float=3)` — 環上的波（Fig. 3.6.2、3.6.5）：r(φ) = radius + amp·f(φ) 的閉曲線（n 點）。
- **class `NH3Model`** — 球棍 NH₃（Fig. 3.7.1），用 View3D 投影到畫面（不用 3D 相機）：N 在上（N_COLOR）、三個 H 在下（H_COLOR），
  - `.positions3(x12)` — [H1, H2, H3, N] 的 3D 位置（分子座標）。
  - `.set_displacements(x12)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
- `irrep_color(label: str)` — 不可約表示標記 → 顏色（"T_{1u}"、"E_g"、"A_1" 都可以；g/u 不分）。
- `rot_matrix(axis, deg)` — 繞 axis 逆時針轉 deg 度的 3×3 矩陣（Rodrigues，浮點；動畫用 rot_matrix(axis, t·deg)）。
- **class `CubeView`** — 投影的立方體（Fig. 4.1.1–4.1.2、4.2.1–4.2.2）：12 條稜、四條上色的體對角線（DIAG_COLORS，端點標 1–4，編號＝ sym_oct.DIAGONALS）、
  - `.redraw()` — 
  - `.set_rotation(Rm)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
  - `.add_axis(axis, color=GREY_B, length: float=1.6, label_tex: str | None=None, font_size: float=22)` — 
  - `.add_arrow(v, color=KET_COLOR, label_tex: str | None=None, body: bool=False, font_size: float=22)` — 
- **class `OctaSites`** — 八面體的六個頂點（Fig. 4.3.1–4.3.2、4.3.4–4.3.7）：|1⟩…|6⟩ 在 +z, −z, +x, −x, +y, −y（＝ sym_oct.SITES），用 View3D 投影。
  - `.point(i: int)` — 
  - `.redraw()` — 
  - `.set_amplitudes(v, min_frac: float=0.25, value_texs=None)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
- `level_diagram(columns, links=(), col_w: float=1.0, gap: float=1.4, unit: float=0.5, font_size: float=22, deg_gap: float=0.07, colors=None, headers=None, header_size: float=20, stroke_width: float=4)` — 能階分裂／相關圖（Fig. 4.2.3、4.3.2、4.3.5、4.4.3）：columns ＝ 每一欄的 [(能量, 標記 TeX, 簡併度), …]；簡併度 d 畫成 d 條間隔 deg_gap
- **class `XF6Model`** — 球棍 XF6（Fig. 4.4.1–4.4.4、4.4.7），用 View3D 投影：中心 X（X_COLOR）、六個 F（F_COLOR，原子 0–5 在 OCT_SITES·bond）、
  - `.positions3(x21)` — [F0, …, F5, X] 的 3D 位置。
  - `.set_displacements(x21)` — 
  - `.arrows(x21, scale: float=1.0, color=KET_COLOR, min_len: float=0.04)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
- **class `TriadBall`** — 投影的球＋本體三腳架（Fig. 5.1.1、5.3.1–5.3.5、5.5.10–5.5.11）：外輪廓（正交投影下永遠是圓）、跟著物體轉的赤道與本初子午線
  - `.project(v)` — 
  - `.redraw()` — 
  - `.set_rotation(Rm)` — 
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
  - `.add_crank(axis, color=CRANK_COLOR, label_tex: str | None=None, length: float=1.6)` — 
  - `.tip(k: int)` — 
  - `.trail(k: int, Rs, color=None, width: float=3)` — 
- **class `VectorCone`** — 角動量向量錐（Fig. 5.4.4、5.5.3、5.5.6）：頂點在原點、軸沿 R_y(β)ẑ、斜高 √(j(j+1))·unit、軸上的高 m·unit（cos θ = m/√(j(j+1))）。
  - `.set_beta(beta: float)` — 
- `cone_family(j, unit: float=0.5, view: View3D | None=None, colors=None)` — Fig. 5.4.4：m = j…−j 的一整族圓錐（同一個頂點）＋ z 軸虛線。colors：{m: 顏色}；預設由暖到冷。→ VGroup，.cones[k] 對應 m = j − k。
- **class `ProbBars`** — 機率／振幅長條（Fig. 5.5.2–5.5.5）：n 根長條、底線、下面的標記。ProbBars(values, labels=None, width=0.32, gap=0.1,
  - `.set_values(values)` — 
  - `.top(i: int)` — 
- `ladder_column(levels, width: float=1.0, unit: float=0.5, font_size: float=20, stroke_width: float=3, left_buff: float=0.1, right_buff: float=0.1)` — 轉子／殼層能階梯（Fig. 5.5.7、5.6.7、5.8.6）：levels = [(能量, 左邊 TeX 或 None, 右邊 TeX 或 None, 顏色), …]；
- `weight_diagram(weights, unit: float=0.6, roots=None, axes: bool=True, axis_labels=('H_1', 'H_2'), extent: float=2.2, dot_radius: float=0.07, color=YELLOW, root_color=TEAL, font_size: float=20, ring_gap: float=0.07)` — 根圖／權圖（Fig. 5.8.2–5.8.5）：weights = {(x, y): 重數} 或 [(x, y), …]（重數 > 1 加同心圓圈）；roots = [(x, y, 標記 TeX 或 None), …]
- `mnemonic_wheel(labels, radius: float=1.5, start_deg: float=180, direction: int=-1, font_size: float=24, colors=None, ring_color=GREY_B, angles=None)` — 助記輪（Fig. 5.6.8、5.6.9、5.7.3、5.7.6）：圓環＋沿圓周的 irrep 標記（TeX）。第 k 個標記在角度 start_deg + direction·k·360/n
- `wheel_marker(wheel, k: int, text: str | None=None, inner: bool=True, color=YELLOW, font_size: float=18)` — 助記輪上第 k 格的箭頭標（inner=True 在圈內指向該格、False 在圈外）；text 以 $ 開頭結尾的當 TeX，其餘用 zh。
- **class `PairOcta`** — 兩個粒子在八面體上（Fig. 6.1.1）：OctaSites 加上粒子 A（大的空心光暈，A_COLOR）與 B（小實心點，B_COLOR）。
  - `.set_pair(i: int, j: int)` — 
  - `.redraw()` — 
- **class `PairGrid`** — 6×6 的兩粒子格子（(6.1.1)、(6.1.8) 的字典序；列 = 粒子 A 的 i、行 = 粒子 B 的 j）：PairGrid(n=6, cell=0.42, font_size=16,
  - `.color_cells(M, palette: dict, opacity: float=0.75)` — 
- **class `StrainCube`** — 投影的立方體（Fig. 6.3.2、6.3.4、6.3.6）：StrainCube(size=1.0, view=None, color=GREY_A, show_ref=True)。
  - `.vertex(v)` — 
  - `.face_center(f: int)` — 
  - `.redraw()` — 
  - `.set_strain(S)` — 
  - `.face_arrows(T, scale: float=0.5, color=STRESS_COLOR, faces=(0, 1, 2, 3, 4, 5), thickness: float=3)` — 
- `field_glyph(kind: str='out', radius: float=0.12, color=BFIELD_COLOR)` — 場與電荷的小記號（Fig. 6.3.10–6.3.11）：kind = "out"（⊙，指出紙面）、"in"（⊗）、"+"、"-"（帶號的電荷圓）。
- `current_loop(center=ORIGIN, rx: float=0.5, ry: float=0.18, sense: int=1, color=BFIELD_COLOR, width: float=3)` — 斜看的電流環（橢圓＋前緣一個箭頭；sense = +1 前緣向右＝從上往下看逆時針，B 向上；−1 反向）→ VGroup(環, 箭頭)。
- `trace_curve(axes, ts, ys, color=RESP_COLOR, stroke_width: float=2.5)` — 時間曲線（Fig. 6.5.1）：plot_axes 上的折線（一條 VMobject；n ≤ 1000 點，靜態畫、ShowCreation，不放進 updater）。
- **class `SmithChart`** — Smith 圖（Fig. 6.5.3）：SmithChart(S=1, gamma=0.1, unit=0.25, deltas=(), extent=None)：實軸、虛軸、Γ 圓（直徑 S/Γ 沿虛軸）、
  - `.point(delta: float)` — 
  - `.get_origin()` — 
  - `.phasor(delta: float, color=YELLOW)` — 
- `hat_contours(center=ORIGIN, unit: float=2.0, j: float=1.0, k: float=1.0, r: float=0.0, levels=(-0.2, -0.1, 0.0, 0.2), n: int=60, colors=None, stroke_width: float=2.5)` — 墨西哥帽（Fig. 6.7.1）下層 ε−(x1, x2) = kρ² − (j²ρ² + 2jrρ³cos3φ + r²ρ⁴)^{1/2} 的等高線（(6.7.8)；自己算，不 import sym_prod）：
- `hat_profile(axes, j: float=1.0, k: float=1.0, r: float=0.0, phi_deg: float=0.0, rho_max: float=1.4, n: int=60, colors=(GREY_A, TEAL), stroke_width: float=3)` — Fig. 6.7.3 的截面：沿 φ（右半，x > 0）與 φ + 180°（左半，x < 0）的 ε+（上層）、ε−（下層）→ VGroup(上層, 下層)（各一條 2n 點的線）。
- **class `SpinTop`** — 小陀螺（Fig. 7.1.1）：球、自旋箭頭（與向上夾 angle 度）、繞箭頭的進動虛線橢圓、可選的標籤（球下方）。
  - `.set_angle(angle: float)` — 
- **class `SpinPair`** — 兩個 SpinTop（質子大、電子小；Fig. 7.1.1 與 H701 的四個乘積態）。.set_ms(m1, m2)：+½ ↑、−½ ↓、0 → 水平（疊加態的示意）。
  - `.set_ms(m1: float, m2: float)` — 
- `vector_triangle(j1, j2, j3, unit: float=0.5, start=ORIGIN, angle: float=0.0, colors=(J1_COLOR, J2_COLOR, J3_COLOR), labels: bool=True, font_size: float=24)` — Fig. 7.1.5：|j1| + |j2| ＝ |j3| 的三角形（邊長 √(j(j+1))·unit）；j3 沿 angle 度，j1 在 j3 的左側（逆時針）。
- **class `CouplingCone`** — Fig. 7.1.6：總角動量 J（長 √(J(J+1))，與 z 軸夾 Θ，cos Θ = M/|J|，在 xz 平面）、j1 繞 J 的錐（夾角 θ13 由餘弦定律）、
  - `.j1_vector(phase: float)` — 
  - `.set_phase(phase: float)` — 
- **class `RESurface`** — RE 曲面的投影（不用真 3D）：每個能量一組等高線（單位向量折線 × 該組半徑），前面亮、後面淡；可加淡色經緯網與本體三軸。
  - `.mesh_lines(radius_fn, n_lat: int=5, n_lon: int=8, n: int=48)` — 淡色經緯網：radius_fn(n×3 單位向量) → 半徑陣列；回傳 [n×3 點]（已乘半徑）。
  - `.set_view(theta: float | None=None, phi: float | None=None)` — 
- **class `PoincareSphere`** — Poincaré–Stokes 球（Fig. 7.5.1(b2)）：SphereView ＋ A、B、C 三軸（A→x、B→y、C→z）與標籤；
  - `.add_stokes(S, color=STOKES_COLOR)` — 
  - `.set_stokes(g: VGroup, S)` — 
- `polarization_ellipse(a: float, b: float, phi: float=0.0, center=ORIGIN, unit: float=1.0, color=POL_COLOR, n: int=60, arrow: bool=True, width: float=3)` — (7.5.20a,b)：E(Φ) = (a cos Φ cos φ − b sin Φ sin φ, a cos Φ sin φ + b sin Φ cos φ)·unit；b > 0 時 Φ 增加是逆時針（右旋）。
- `ellipse_from_stokes(S, center=ORIGIN, unit: float=1.0, color=POL_COLOR, **kw)` — Stokes (S_A, S_B, S_C)（(7.5.25) 圓欄）→ ψ = ½asin(S_C/I)、φ = ½atan2(S_B, S_A)、a = √I cos ψ、b = √I sin ψ。
- `projected_curve(points3, view: View3D, color=TEAL, width: float=2.5, front: float=1.0, back: float=0.3)` — 3D 折線（例如振子相空間 (x1, x2, p1) 的軌跡，Fig. 7.5.9–7.5.10）投影；前（depth ≥ 0）亮、後淡。
- `box_frame(X: float, Y: float, center=ORIGIN, unit: float=1.0, angle: float=0.0, color=GREY_B, width: float=1.6)` — 2X×2Y 的虛線盒（Fig. 7.5.5 的線偏振盒；Fig. 7.5.6 的傾斜盒用 angle，度）。
- `potential_curves(axes, xs, ys_list, colors=(POT_G_COLOR, POT_U_COLOR), labels=None, label_x=None, font_size: float=22, width: float=3)` — plot_axes 上的位能曲線（Fig. 7.6.2）：每條一個 VMobject（建議 ≤ 60 點），標籤放在 x = label_x（預設最右）。→ VGroup(curves, labels)。
- `min_marker(axes, x: float, y: float, color=YELLOW, label_tex: str | None=None, font_size: float=20, drop: bool=True)` — 極小點：點＋（可選）落到 y = 0 的虛線＋標籤（點的右下方）。
- `polyline(axes, pts, color=TRAJ_COLOR, width: float=2.5, dashed: bool=False, num_dashes: int=40)` — plot_axes 上的一條折線（pts：(n, 2) 資料座標；建議 ≤ 60 點、靜態）；dashed=True 畫成虛線。
- `line_family(axes, lines, x0: float, x1: float, color=GREY_B, width: float=1.6, highlight=None, hcolor=YELLOW)` — 直線族 Y = sX + b（lines = [(s, b)]，sym_dyn.line_family 的輸出；Fig. 8.1.2）在 x ∈ [x0, x1] 的線段；highlight = 要標亮的序號。
- `trajectory_fronts(axes, trajs, fronts, traj_color=TRAJ_COLOR, front_color=FRONT_COLOR, traj_width: float=1.6, front_width: float=2.5)` — Fig. 8.1.6：軌跡族（trajs：[(n, 2)]，sym_dyn.trajectory）＋ 等 S_H 波前（fronts：[(n, 2)] 或 {level: [段…]}，sym_dyn.fig816_curves）。
- **class `PhasePhasor`** — 相空間的相量（Fig. 8.2.1）：圓心 center、半徑 radius 的虛線軌道（|α| = amp 對到 radius）、兩軸、從圓心到 α 的箭頭、
  - `.point(z)` — 
  - `.set_alpha(z)` — 
- `snapshot_stack(xs, rows, unit_x: float=0.45, unit_y: float=0.5, gap: float=0.45, top=ORIGIN, re_color=PACKET_COLOR, env_color=PACKET_ENV_COLOR, width: float=2.0, envelope: bool=True, baseline: bool=True)` — 波包快照列（Fig. 8.2.2、8.3.2(a)）：rows = [複數陣列（在 xs 上）]，第 k 列畫在 y = top_y − k·gap：實線 Re χ、虛線 |χ|（envelope）、
- **class `CrankSphere`** — 曲柄向量 ω 與自旋向量 S 的 3D 圖（Fig. 8.5.1、7.5.3）：SphereView（球面代表 |v| = ½，即自旋向量的長度）＋ 三軸與標籤；
  - `.vec_point(v)` — 
  - `.add_vector(v, color=SPIN_COLOR, label_tex: str | None=None, font_size: float=22)` — 
  - `.set_vector(g: VGroup, v)` — 
  - `.add_trail(points3, color=SPIN_COLOR, width: float=2.0, front: float=1.0, back: float=0.3)` — 
  - `.cone_rim(axis, half_angle: float, length: float=0.5, n: int=48)` — 
- `dressed_ladder(blocks, unit: float=1.0, base=ORIGIN, bare_w: float=0.7, dressed_w: float=0.8, gap: float=0.55, font_size: float=18, labels: bool=True)` — Jaynes–Cummings 穿衣能階梯（Fig. 8.7.1）：blocks = [(N, 下穿衣, 上穿衣, 裸 |N⟩|↓⟩, 裸 |N−1⟩|↑⟩)]（sym_dyn.jc_ladder 的輸出；

## `combkit.draw.walk_draw`

walk_draw.py —— 《Elements of the Random Walk》第 1 章的畫圖元件（ManimGL）。數學在 walk1d.py（先算、再畫）。

常數：`HEAD_COLOR`、`TAIL_COLOR`、`WALK_COLOR`、`BAR_COLOR`、`GAUSS_COLOR`、`EXACT_COLOR`、`WRONG_COLOR`、`RIGHT_COLOR`、`MUTED`

- `coin(face: str='H', radius: float=0.3)` — 
- `coin_row(coins: str, radius: float=0.3, buff: float=0.12)` — 
- `number_line(N: int, unit: float=0.5, label_every: int=2, font_size: float=22, color=MUTED, labels: bool=True)` — 
- `line_point(nl: VGroup, x: float)` — 數線 nl 上座標 x 的點（nl 可以移動）。
- `walker_dot(color=WALK_COLOR, radius: float=0.12)` — 
- `path_point(k: float, x: float, ux: float, uy: float, origin=ORIGIN)` — 
- `walk_axes(N: int, ymax: int, ux: float=0.5, uy: float=0.4, origin=ORIGIN, x_every: int | None=None, y_every: int | None=None, font_size: float=20, x_label: str='步數', y_label: str='位置', color=MUTED)` — 
- `walk_path(pos, ux: float=0.5, uy: float=0.4, origin=ORIGIN, color=WALK_COLOR, stroke_width: float=4, dots: bool=True, dot_radius: float=0.05)` — 
- `fan(pos_list, ux: float=0.05, uy: float=0.2, origin=ORIGIN, color=BLUE_B, stroke_width: float=1.0, opacity: float=0.35)` — 
- `sqrt_envelope(Nmax: int, ux: float, uy: float, origin=ORIGIN, c: float=1.0, color=GAUSS_COLOR, stroke_width: float=3, samples: int=120)` — 
- `bar_chart(values: dict, ux: float=0.4, uy: float=1.0, origin=ORIGIN, color=BAR_COLOR, width: float=0.8, fill_opacity: float=0.75, labels: bool=False, label_fmt=str, font_size: float=20, label_color=WHITE)` — 
- `bar_of(chart: VGroup, x)` — 
- `curve(fn, x0: float, x1: float, ux: float=0.4, uy: float=1.0, origin=ORIGIN, color=GAUSS_COLOR, stroke_width: float=3, samples: int=200, dashed: bool=False, num_dashes: int=40)` — 
- `plot_axes(x_range, y_range, width: float=6.0, height: float=3.5, x_ticks=(), y_ticks=(), x_tick_labels=None, y_tick_labels=None, font_size: float=20, color=MUTED, x_label: str | None=None, y_label: str | None=None)` — x_range=(a, b, step)。刻度數字用 Tex 字串（x_tick_labels 可給對應的 TeX，否則用 str(值)）。
- `pascal_triangle(rows, dx: float=0.32, dy: float=0.5, origin=ORIGIN, font_size: float=24, color=WHITE)` — 
- `pascal_entry(tri: VGroup, N: int, y: int)` — 
- `delta_grid(N: int, cell: float=0.55, origin=ORIGIN, diag_color=HEAD_COLOR, off_color=GREY_D)` — 
- `grid_cell(g: VGroup, j: int, k: int)` — 
- `tex_table(rows, col_widths, row_h: float=0.5, font_size: float=26, header_color=YELLOW, line_color=GREY_C)` — rows：list of list，元素是 TeX 字串或 Mobject（中文表頭用 zh(...)）。回傳 VGroup(格子們, 橫線)，左上角在 ORIGIN。
- `table_cell(t: VGroup, i: int, j: int, ncols: int)` — 
- `erratum_box(where: str, wrong_tex: str, right_tex: str, width: float=5.6, font_size: float=30)` — 
- `contour_circle(axes: Axes, r: float, color=TEAL, stroke_width: float=3)` — 
- `pole_marker(axes: Axes, z: complex, label_tex: str='', color=RED, font_size: float=26)` — 

## `combkit.draw.young_art`

young_art.py —— Young 表的藝術桌布工具（來自 youngart：Fulton《Young Tableaux》全書 60 張 4K 桌布，codes/art/ykit.py）

常數：`PARTS`

- **class `YPlate`** — 與 art_kit.Plate 相同，只有右上角的標籤用 PARTS（章或附錄）。
  - `.caption()` — 
- `shape_of(T)` — 
- `cell_center(i, j, cs, origin)` — 
- `cell_squares(cells, cs=0.3, origin=(0, 0), colors=WHITE, opacity=0.85, stroke=None, sw=1.0, gap=0.08, stroke_opacity=0.9)` — 
- `tab_squares(T, cs=0.3, origin=(0, 0), cmap=None, vmax=None, vmin=1, **kw)` — 
- `outline_pts(lam, cs=0.3, origin=(0, 0))` — 
- `perm_pts(w, size=1.0, center=ORIGIN)` — 
- `rot3(P, yaw=0.6, pitch=0.35)` — 
- `proj(P, yaw=0.6, pitch=0.35, scale=1.0, center=ORIGIN)` — 
- `turn3d(items, view, angle=0.7, T=T_B)` — 動畫：items = [(物件, 3D 折線), …]（物件是用 proj 畫出的 path／glow 群組或點雲），view = (yaw, pitch, scale, center)；

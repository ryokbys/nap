# Extended Lagrangian Shell Model (XL-SM) 実装計画

## 背景と目的

現在の adiabatic 実装ではshell粒子に微小質量を与えてvelocity-Verletで積分する。
shell振動周期 τ_s = 2π√(m_s/k2s) が小さいため dt ≤ 0.16 fs 程度が必要で非常に遅い。

Nomura et al., CPC 192 (2015) の extended Lagrangian (XL) 法を shell model に適用すると：
- 補助変数 θ（shell位置の時間可逆アナログ）が1回の勾配ステップの初期推定値を提供
- long-time エネルギー保存性を維持しつつ coreに合わせた大きな dt（0.5–2 fs）が使用可能
- 2回のforce評価/stepが必要だが dt×10倍で実質5倍の高速化が期待される

## アルゴリズム（1MDステップ）

パラメータ（自動設定）：
- ω² = 2/dt²  （K=2、Nikolassonの推奨値）
- κ = 1/max(sm_k2s)  （最大バネ定数の逆数；調和ポテンシャルで厳密なNewtonステップ）

```
1. 第1キック（coreのみ）:   va_core += aa_core * fa2v_core * dt
2. core位置更新:             ra_core += (hi*va_core)*dt
3. θ予測ステップ（fractional座標）:
     θ_dot_half = θ̇ + 0.5*dt * θ̈
     θ(t+dt)    = θ(t) + dt * θ_dot_half
     ra_shell   = θ(t+dt)            ← shell をθ位置に置く
4. force評価 #1: get_force() at (r_c(t+dt), θ)
     → F_core と F_shell が aa(:,i) に格納される
5. 勾配ステップ（shellのみ）:
     δra = hi*(κ * aa_shell)         ← Cartesian力 → fractional変位
     ra_shell = θ(t+dt) + δra
6. θ修正ステップ:
     θ̈(t+dt) = ω² * δra
     θ̇(t+dt) = θ_dot_half + 0.5*dt * θ̈(t+dt)
7. force評価 #2: get_force() at (r_c(t+dt), r_s*(t+dt))  ← 正確なF_core
8. 第2キック（coreのみ）:   va_core += aa_core * fa2v_core * dt
```

## 変更ファイルと変更内容

### 1. `pmd/force_ShellModel.F90`（主要変更）
新モジュール変数：
- `logical :: use_xl = .false.`
- `real(rp),allocatable,save:: xl_theta(:,:), xl_thdot(:,:), xl_thacc(:,:)`  (3,namax)
- `real(rp):: xl_omega2, xl_kappa`

新サブルーチン：
- `xl_init(namax, ra, tag_isp, dt)`: θ配列確保・初期化、ω²・κ自動設定
- `xl_predict(natm, tag_isp, ra, dt)`: θ予測ステップ→ra_shell=θにセット
- `xl_gradient_step(natm, tag_isp, ra, aa, hi, dt)`: 勾配ステップ + θ̈・θ̇修正

### 2. `pmd/pmd_core.F90`
- `use` に `use ShellModel, only: use_xl, is_shell_sp, xl_init, xl_predict, xl_gradient_step` 追加
- VVループ直前（初期化フェーズ）に `if(use_xl) call xl_init(...)` 追加
- 第1キック・位置更新ループに shell skipを追加（`if(use_xl .and. is_shell_sp(is)) cycle`）
- 位置更新直後に `if(use_xl) call xl_predict(...)` 追加
- get_force()呼び出し直後に XL gradient step と 2回目のget_force()を追加
- 第2キックループに shell skipを追加

### 3. `pmd/read_input.F90`
- `shell_dynamics` キーワード追加（`'xl'` → `use_xl_shell = .true.`）

### 4. `pmd/mod_pmdvars.F90`
- `logical:: use_xl_shell = .false.` を追加

### 5. `pmd/pmd_main.F90`（重要）
- `bcast_params` サブルーチン末尾に以下を追加：
  ```fortran
  call mpi_bcast(use_xl_shell, 1, mpi_logical, 0, mpicomm, ierr)
  ```
  （read_inputはrank 0のみが読むため、他のMPIランクへのブロードキャストが必須）

### 6. `pmd/force_common.F90`
- init_force 内の ShellModel 初期化ブロックに `use_xl = use_xl_shell` を追加

### 7. `tmp/tmp_shell-model/BaTiO3/in.pmd`
- `shell_dynamics xl` を追加
- `time_interval` を 1.0 fs に変更（adiabaticの10倍）
- shell speciesを fix（XLではVV不要）

## 検証方法

1. **BaTiO3 NVE**：XLモード(dt=1.0 fs, 1000 steps)でエネルギードリフト < 0.01 eV/atom/ps
2. **adiabatic比較**：同じ物理時間(100 fs)で adiabatic(dt=0.1 fs) vs XL(dt=1.0 fs)の wall-clock時間比較
3. **単純spring テスト**：O+O_s 2原子系でθがshell位置を追随することを確認

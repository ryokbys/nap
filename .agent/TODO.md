# TODO: Extended Lagrangian Shell Model 実装

## 実装ステータス一覧

### Phase 1–4: 基盤実装（完了）

- [x] `pmd/mod_pmdvars.F90`: `use_xl_shell` フラグ追加
- [x] `pmd/force_ShellModel.F90`: `use_xl`, `xl_theta/thdot/thacc` 配列追加
  - xl_theta/thdot/thacc は `(3,ntot)` で確保し、`tag_itot(i)` でインデックス
    （bamove によるarray再配置後も valid になるよう修正済み）
- [x] `pmd/force_ShellModel.F90`: `xl_init(namax,natm,ra,tag_isp,dt)` 実装
- [x] `pmd/force_ShellModel.F90`: `xl_predict(natm,tag_isp,ra,dt)` 実装
- [x] `pmd/force_ShellModel.F90`: `xl_gradient_step(natm,tag_isp,ra,aa,hi,dt)` 実装
- [x] `pmd/read_input.F90`: `shell_dynamics` キーワード追加
- [x] `pmd/pmd_main.F90`: `bcast_params` に `use_xl_shell` の MPI_Bcast 追加
- [x] `pmd/force_common.F90`: `use_xl = use_xl_shell` 転写追加
- [x] `pmd/pmd_core.F90`: `use ShellModel` + VVループへの xl_init/predict/gradient_step 組み込み
- [x] `tmp/tmp_shell-model/BaTiO3/in.pmd`: `shell_dynamics xl`、dt=1.0 fs、100 steps

### Phase 5: デバッグ（作業中、ブロック中）

- [x] コンパイル通過（`cd pmd && make pmd`）
- [x] バグ修正: `bamove()` 後に xl_theta が無効化される問題
  → `xl_theta/thdot/thacc` のインデックスを `i`（配列位置）から `tag_itot(i)`（グローバルID）に変更
- [ ] **[現在のブロッカー] XL エネルギー発散バグの修正**
  - 症状：BaTiO3 NVE、dt=1.0 fs、T=300K で istp=2 (t=2 fs) にepot が +271 eV 急騰
  - istp=1 は完璧にエネルギー保存（Δetot < 0.01 eV）
  - istp=2 から catastrophic divergence が始まる
  - 詳細は `log_260525.md` を参照
- [ ] エネルギー保存性確認（Δetot/atom/ps < 0.01 eV が目標）
- [ ] adiabatic (dt=0.1 fs) との wall-clock 時間比較

## 現在調査すべき事項（次のエージェントへ）

1. **istp=2 のエネルギー急騰の原因特定**（最優先）
   - `xl_gradient_step` が istp=2 でシェルを間違った位置に移動している疑い
   - または `xl_theta` の更新ロジックに累積誤差がある可能性
   - 詳細は `log_260525.md` の「現在の課題」セクションを参照

2. **修正候補（試す順）**
   - (A) 診断出力を増やして step2 での xl_theta, ra_shell, epot_#1, epot_#2 を確認
   - (B) multi-step gradient iteration：シェルをN回反復して平衡位置に収束させる
   - (C) xl_step_factor をさらに小さくして安定性確認（現在 0.1）
   - (D) そもそも xl_gradient_step の formula を見直し

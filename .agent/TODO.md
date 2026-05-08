# TODO: 単精度・倍精度切り替え実装

## Step 1: 精度定義モジュールの新設

- [x] `pmd/mod_precision.F90` を新規作成し，`rp` kind定数を定義する

## Step 2: configure.ac の修正

- [x] `configure.ac` に `--enable-single` オプションを追加し，`-D__SINGLE__` フラグを渡す

## Step 3: 浮動小数点変数宣言の置換

- [x] 全ソースファイルに `use mod_precision` を追加（未記載ファイルのみ）
- [x] `real(8)` → `real(rp)` に置換（約2,422箇所）
- [x] `real*8` → `real(rp)` に置換（24箇所）
- [x] `implicit real*8(a-h,o-z)` → `implicit real(rp)(a-h,o-z)` に置換（5箇所）

## Step 4: 物理定数リテラルの修正

- [x] `.h` ヘッダファイルの `real(8),parameter::` → `real(rp),parameter::` に置換
- [x] 倍精度指数リテラル（`1.234d0`）→ `real(1.234d0, rp)` に変換
- [x] ヘッダインクルード前に `use mod_precision` が記述されていることを確認

## Step 5: MPI型定数の置換

- [x] MPI使用ファイルに `use pmdmpi` を追加（未記載ファイルのみ）
- [x] `mpi_real8` → `mpi_real_rp` に置換
- [x] `MPI_REAL8` → `mpi_real_rp` に置換
- [x] `mpi_double_precision` → `mpi_real_rp` に置換
- [x] `MPI_DOUBLE_PRECISION` → `mpi_real_rp` に置換

## Step 6: ビルド・テスト

- [x] `autoconf` を実行してビルドシステムを再生成
- [x] 倍精度モード（デフォルト）でビルドし，既存テストが通ることを確認（既存バグは変化なし）
- [x] 単精度モード（`--enable-single`）でビルドし，`rp=4` となることを確認
- [ ] 両モードで簡単なMDシミュレーションを実行し，結果が定性的に一致することを確認

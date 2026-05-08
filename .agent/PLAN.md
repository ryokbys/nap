# 単精度・倍精度切り替え実装プラン

## 概要

`configure` 時に `--enable-single` を指定した場合に，pmd全体で単精度浮動小数点を使用するよう改変する．

---

## Step 1: 精度定義モジュールの新設

`pmd/mod_precision.F90` を新規作成し，kind定数 `rp`（real precision）を定義する．

```fortran
module mod_precision
#ifdef __SINGLE__
  integer,parameter:: rp = 4
#else
  integer,parameter:: rp = 8
#endif
end module mod_precision
```

---

## Step 2: `configure.ac` に `--enable-single` オプションを追加

`configure.ac` に以下を追加し，`-D__SINGLE__` コンパイラフラグを渡す．

```autoconf
AC_ARG_ENABLE([single],
  [AS_HELP_STRING([--enable-single], [Use single precision floating point])],
  [if test "x$enableval" = "xyes"; then
    FCFLAGS="$FCFLAGS -D__SINGLE__"
  fi])
```

---

## Step 3: 浮動小数点変数宣言の置換

全 `.F90` / `.F` ファイルで以下を行う．

1. `use mod_precision` を追加（未記載のファイルのみ）
2. `real(8)` → `real(rp)` に置換（2,422箇所）
3. `real*8` → `real(rp)` に置換（24箇所）
4. `implicit real*8(a-h,o-z)` → `implicit real(rp)(a-h,o-z)` に置換（5箇所）

---

## Step 4: 物理定数リテラルの修正

`.h` ヘッダファイル（`params_unit.h`, `params_au.h`, `params_LJ.h` 等）の定数定義を修正する．

- `real(8),parameter::` → `real(rp),parameter::` に置換
- 倍精度指数リテラル（`1.234d0`）→ `real(1.234d0, rp)` に変換

なお，ヘッダファイルは `use mod_precision` でなく `include` で取り込まれるため，
`rp` が参照できるよう `use mod_precision` がヘッダ読み込み前に記述されていることを確認する．

---

## Step 5: MPI型定数の置換

### 既存インフラ（変更不要）

`mod_pmdmpi.F90` に既に以下の定義がある：

```fortran
#ifdef __SINGLE__
  integer,parameter:: mpi_real_rp = mpi_real4
#else
  integer,parameter:: mpi_real_rp = mpi_real8
#endif
```

### 置換作業

全 `.F90` / `.F` ファイルで以下を行う（556箇所・58ファイル）．

1. `use mod_pmdmpi` を追加（未記載のファイルのみ）
2. `mpi_real8` → `mpi_real_rp` に置換
3. `MPI_REAL8` → `mpi_real_rp` に置換
4. `mpi_double_precision` → `mpi_real_rp` に置換
5. `MPI_DOUBLE_PRECISION` → `mpi_real_rp` に置換

---

## Step 6: ビルド・テスト

1. `autoconf` / `./configure --enable-single` を実行してビルドシステムを再生成
2. 倍精度モード（デフォルト）でビルドし，既存のテストが通ることを確認
3. 単精度モード（`--enable-single`）でビルドし，コンパイルエラーがないことを確認
4. 簡単なMDシミュレーションを両モードで実行し，結果が定性的に一致することを確認

---

## ファイル変更サマリ

| 対象 | 変更内容 |
|------|----------|
| `pmd/mod_precision.F90` | **新規作成**：`rp` kind定数の定義 |
| `configure.ac` | `--enable-single` オプションの追加 |
| `pmd/*.F90`, `pmd/*.F` | `real(8)` → `real(rp)`，`use mod_precision` 追加 |
| `pmd/*.h` | `real(8)` → `real(rp)`，リテラル修正 |
| `pmd/*.F90`（MPI使用箇所） | `mpi_real8` 等 → `mpi_real_rp`，`use mod_pmdmpi` 追加 |

---

## 注意事項

- `rp` はコンパイル時定数のため，ランタイムの分岐は一切不要
- `mpi_real_rp` の定義は `mod_pmdmpi.F90` に既存のため変更不要
- ヘッダファイル内で `rp` を参照する際は，インクルード先の `use mod_precision` の順序に注意
- 単精度では精度が落ちるため，収束判定の閾値等を緩める必要が生じる可能性がある

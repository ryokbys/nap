# Plan: debug resize_iarr2 heap corruption (260526)

## 問題
`force_ShellModel` + extended Lagrangian ダイナミクス実行中に，`realloc_namax_related` 内の
`resize_iarr2(tag_igrp, ...)` 呼び出し時にヒープ破壊でクラッシュする．

## 原因特定
`mod_util.F90` の `resize_iarr2` において，配列コピー時に `size(iarr)` を使っているが，
これは2次元配列の全要素数（`ngrpmax * namax = 4 * 1296 = 5184`）を返す．
意図は第2次元サイズ（`namax = 1296`）であり，`size(iarr,2)` を使うべき．

`resize_darr2` / `resize_darr3` は正しく `size(darr,2)` / `size(darr,3)` を使っており，
整数版(`resize_iarr2`, `resize_iarr3`)のみが間違っていた．

## 修正方針
`mod_util.F90` の `resize_iarr2` と `resize_iarr3` で `size(iarr)` を
正しい次元指定 `size(iarr,2)`, `size(iarr,3)` に変更する．

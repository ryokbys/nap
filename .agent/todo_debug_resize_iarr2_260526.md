# TODO: debug resize_iarr2 heap corruption (260526)

- [x] out.pmd を読んでエラー内容を確認
- [x] バックトレースから `resize_iarr2` → `bacopy_` → `realloc_namax_related` の呼び出し経路を特定
- [x] `mod_util.F90` の `resize_iarr2` のバグ（`size(iarr)` を全要素数として使用）を確認
- [x] `resize_darr2` との比較で誤りを確認（darr版は `size(darr,2)` を正しく使用）
- [x] `resize_iarr2`: `size(iarr)` → `size(iarr,2)` に修正
- [x] `resize_iarr3`: `size(iarr)` → `size(iarr,3)` に修正
- [ ] ビルドして動作確認（ユーザーに委ねる）
- [ ] commit（ユーザーの許可を得てから）

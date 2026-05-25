# AI AGENTへの指導書

napパッケージをAI Agentが改変する際のルールを記す．
napパッケージには，pmdという分子動力学プログラムと，nappyという前後処理のpythonパッケージが含まれている．
ルートディレクトリに存在する次のディレクトリは無視して良い．
- `JOSS_paper`
- `mkconf`
- `neb`
- `not_used`
- `qmcl`


## プログラミング言語に関するルール

### pmd

- `pmd/`内には，Fortranプログラムpmdに関するソースコードを置く．
- Fortranプログラムは基本的にはFortran90に準拠する．ただし，2002の機能も必要であれば用いても良い．
- インデントはスペース２文字とする．タブは仕様しない．
- １行をできるだけ78文字とし，長い行は行末に `&` をつけて折り返す．

### nappy

- `nappy/`以下にはPythonパッケージnappyが格納されている．
- 基本的にはpython3系（3.9以上）でコードを書く．


## gitに関するルール

- ai-devブランチで作業する．
- git commitをする際には許可を得る．
- できるだけ作業の塊（１つの機能追加や１つのバグ修正）ごとにcommitする．多くの作業をいくつもこなしてからcommitしない．

## 作業報告

- 何かしらの実装や改変を行う際，作業プランを`.agent/plan_NAME_YYMMDD.md`に，その作業のtodoリストを`.agent/todo_NAME_YYMMDD.md`に，その作業ログを`.agent/log_NAME_YYMMDD.md`に保存する．ただし，NAMEは作業を簡潔に表す文言で，YYMMDDは作業を始めた日付．

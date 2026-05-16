{TARGET_DIR}　＝ {../}（この`DOC-GUIDE.md`ファイルを起点とする．）

# pmdプログラムの構造をJSONとHTMLファイルに書き出す

このファイルは，Fortranプログラムであるpmdの構造をJSONとHTMLファイルに書き出すための指導書．
目的は，実際にFortranファイルを読まずともJSONファイルを読めば，大まかな構造を把握することができるようにする．
HTMLファイルはそれをネットワークダイアグラム表示して，人間が分かりやすくするためのもの．

## 出力ファイル
- `./pmd_structure.json` -- pmdプログラムの構造（subroutine/functionの間の関係）
- `./diagram.html` -- `fetch('./pmd_structure.json')`で読み込み描画．**file:// 直開きでも動くよう同 JSON を `<script type="application/json" id="workflow-data">` に埋め込んでフォールバック．**

## HTMLのビジュアル仕様（ダーク + ウォームゴールド）
**CSS 変数:**
```css
--bg:#0c1117; --bg-2:#11171f; --panel:#161d27; --line:#2a3340;
--text:#d6dde6; --muted:#8a94a3;
--edge:#2c3441; --edge-active:#ffd166; --edge-active-glow:rgba(255,209,102,0.55);
```
グループごとに `stroke` 色と、彩度を落とした近似 `fill` を持たせる（例: stroke `#4cc9f0` / fill `#11242e`）。

**CSS Grid 全体レイアウト:**
```
grid-template-columns: 320px 1fr;
grid-template-rows:    auto 1fr 6px var(--anno-h, 300px);
grid-template-areas: "header header" "sidebar canvas" "sidebar resizer" "sidebar annotations";
```
- header（タイトル + パス + hint）
- 左 sidebar 320px（フローボタン + legend + reset）
- canvas（SVG）
- **6px の resizer**（縦ドラッグで `--anno-h` 更新、`localStorage` に保存）
- 下部 annotations（番号付きステップリスト）

**SVG レイヤー順序（最重要・Z 問題の根源）:**
必ず以下の順で `<g>` を生成すること。バッジを最上層にしないと Node に隠れてホバーが Node に奪われる事故が起きる。
```
defs → col-layer → edge-layer → node-layer → badge-layer
```
badge-layer は node-layer の **後** に追加して常に最前面に置く。

**ノード:** 角丸 rect (rx=8) + drop-shadow + 1.5px stroke、内側に title (12px bold) と subtitle (10.5px muted)。`in-flow` 状態は stroke 2.5px + glow、`dim` 状態は opacity 0.18。

**エッジ:** SVG bezier。エンドポイントは「同カラムなら上下、別カラムなら左右」。ベース `--edge`、アクティブ `--edge-active` + glow + 太矢印マーカー。アクティブ用 marker を別途 `<defs>` に。

**往復エッジの分離（lane offset・必須）:**
- `A→B` と `B→A` が両方存在するペアを検出 → 片方に `+lane`、もう片方に `-lane`（lane ≒ 38px）を割り当てる
- 縦エッジなら cp1/cp2 を x 方向に lane だけずらす、横エッジなら y 方向にずらす
- これをやらないと往復エッジが同じ曲線上に重なってバッジ同士がスタックする
- 割り当ては `fromId < toId` のような決定的な lex 比較で行うとデバッグ容易

**カラム見出し:** viewBox 上端に letter-spacing 広め、中央揃えで `EXTERNAL` 等を配置。`<line stroke-dasharray="3,5">` で divider。

**バッジ（最重要・かぶり対策をここで実装）:**
- **1 エッジ = 1 バッジ**。同じ from→to に複数ステップが乗っても 1 個に集約
- 形状は pill（rounded rect, rx=h/2）。中身を `getBBox()` で測って幅を伸縮
- ラベル: 1 ステップ→`"3"` / 2〜5 ステップ→`"3 · 5 · 7"` / 6+ ステップ→`"3–21 (4)"`
- **配置オフセット（広めに用意・15 個推奨）:**
  `[0.50, 0.44, 0.56, 0.38, 0.62, 0.32, 0.68, 0.26, 0.74, 0.20, 0.80, 0.15, 0.85, 0.10, 0.90]`
- **二段階の衝突判定（必須・片方だけでは隙間に潜む）:**
  1. 既配置 pill との楕円距離（w/h で正規化）が 1.05 以上か
  2. **Node の bbox に重なっていないか**（margin 4px）。Node 内に落ちる候補には大きなスラックペナルティ（-3.5 など）を与える
  - 全候補がダメなら、最良スコアの t を採用したうえで `getPointAtLength` で path 上を ±方向に微小ステップ歩いて Node bbox から脱出する **escape ループ**を最後の砦に入れる
- 配置順は「マルチステップ→短いパス」を先に処理して midpoint を確保、最後にステップ順で DOM 描画
- midpoint から離れて置かれたバッジには **破線テザー** を中点へ引く
- **DOM 構造（3 重ネスト）:**
  ```html
  <g class="edge-step-badge" transform="translate(x,y)">
    <ellipse class="edge-step-halo" .../>              <!-- focus 時の拡散リング -->
    <rect class="edge-step-hit" pointer-events="all"/>  <!-- 透明 hit-area (±8px padding) -->
    <g class="edge-step-scale">                         <!-- pulse 用の内側スケール -->
      <rect class="edge-step-pill" .../>                <!-- 可視 pill -->
      <text>...</text>
    </g>
  </g>
  ```
  - スケール用 CSS transform は内側 `.edge-step-scale` に当てて外側 translate と競合させない
  - `.edge-step-pill` にクラスを切って CSS セレクタを `:hover .edge-step-pill { fill: ... }` のように**可視 pill だけにスコープ**。`.edge-step-badge:hover rect` のようにタグ単位で書くと hit-area まで色が乗る
  - 透明 hit-area は `fill: rgba(0,0,0,0.001)` + `pointer-events: all` で確実に拾う

**イベントハンドラ:** SVG の `<g>` には `mouseenter` ではなく **`mouseover` / `mouseout`（bubble する版）** を使う。`mouseout` は `outer.contains(ev.relatedTarget)` で同一バッジ内の移動を除外する。

**フォーカス（PIKAPIKA・5 層演出）:**
注釈クリックやバッジクリックで focusedStep が立ったら、以下を **同時に**走らせる：
1. **バッジ pulse**: `.edge-step-scale` に `pulse` keyframe（scale 1.0 → 1.42、`drop-shadow` を 3 段重ね＝白コア + ゴールド + ハロー、0.9s ease-in-out cubic-bezier 強め、infinite）
2. **バッジ色反転**: focus 時に pill の fill を `--edge-active`、stroke を `#fff8d6`、text を `#1a1100` に
3. **halo リング**: `<ellipse class="edge-step-halo">` を `halo-burst` keyframe で 0.55x → 3.0x にスケールしながら opacity 0.85 → 0 で fade。1.3s infinite
4. **エッジ脈動**: 該当 `<path>` に `.focused-edge` クラス → `edge-pulse` keyframe で stroke-width 2.8 ↔ 4.2、色 gold ↔ 白、`drop-shadow` 強化。1.1s
5. **両端ノード**: `.endpoint` クラス → stroke 3px に上げ、`node-glow` keyframe で drop-shadow が gold-glow ↔ 純白でパルス。1.1s
6. **注釈 li フラッシュ**: `.focused` クラス + `anno-flash` keyframe（1.2s で背景が強い金 → 落ち着いた金へワンショット）。**同じ step を再クリックしたら再生し直す**ため、`li.classList.remove → 強制 reflow (li.offsetWidth 参照) → add` の順で再付与する

各 keyframe は `transform-origin: center; transform-box: fill-box;` を併用して、SVG 要素の自分自身の中心でスケールするようにする。

**注釈パネル:** CSS counter で番号付き `<ol>`、各 li は `from → to`（等幅）+ `passes`（背景パネル + 左 2px のゴールドボーダー）+ `note`（muted）。クリックで該当ステップへ pan + 上記 PIKAPIKA + `scrollIntoView`。

**ツールチップ:** 固定位置の `<div id="tooltip">`。背景 `rgba(14,20,28,0.97)` + ゴールド枠 + backdrop-filter blur。バッジホバーで出す：
- 1 ステップ: flow 名+ステップ番号、`from → to` チップ、passes コードブロック、note
- マルチステップ: 区切り線でステップを連結表示
- カーソル右下 +16px、画面外なら反対側に反転
- **フロー選択中はバッジ専用にし、Node ホバーでは出さない**（`if (activeFlowId) return`）


## インタラクション
- ホイール = カーソル位置をピボットにズーム（SVG viewBox を JS で操作）
- ドラッグ = パン（mousedown→mousemove→mouseup）
- バッジクリック = 単体ならフォーカス、**マルチステップは循環フォーカス**（クリックするたびにステップ番号を巡回）
- 注釈 li クリック = フォーカス + 5 層 PIKAPIKA + scrollIntoView
- 数字キー `1`〜`9` `0` でフロー選択、`F`=全体表示、`+`/`-`=ズーム、`Esc`=フォーカス解除 → もう一度押すとフロー全体解除
- フロー切り替え時にフィットアニメ（該当ノード群の bbox に padding 60px でフィット）

## やってはいけないこと

- 想像のsubroutine/function/connectionを作る

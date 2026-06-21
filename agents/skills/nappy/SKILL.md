---
name: nappy
description: >-
  nappy post-processing toolkit for MD trajectory files produced by pmd.
  Sub-skills handle RDF (radial distribution function, g(r)) and MSD
  (mean square displacement, diffusion coefficient).
  Use whenever the user asks to analyse MD trajectories or compute physical
  quantities from extxyz or pmd output files.
  Trigger on: "RDF", "MSD", "g(r)", "radial distribution function",
  "mean square displacement", "diffusion", "trajectory analysis",
  "RDFを計算", "RDFをプロット", "動径分布関数",
  "MSDを計算", "MSDをプロット", "平均二乗変位", "拡散のMSD", "MSDから拡散係数".
---

# nappy Skill

nappy は pmd MD シミュレーション結果の Python 後処理パッケージ。
このスキルはトラジェクトリ解析タスクを担う。ユーザーの要求に応じて
下記のサブスキルファイルを読み込み、その指示に従って実行する。

## サブスキル一覧

| タスク | 読み込むファイル |
|--------|----------------|
| 動径分布関数 (RDF / g(r)) | `skill://nappy/rdf.md` |
| 平均二乗変位 (MSD / 拡散係数) | `skill://nappy/msd.md` |

該当するサブスキルを特定し、`read` ツールで読み込んでから Step 1 以降の指示に従うこと。

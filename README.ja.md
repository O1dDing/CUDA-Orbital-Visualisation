# Chemical Orbital Visualiser (COV)

[English](README.md) · [简体中文](README.zh-CN.md) · **日本語** · [Français](README.fr.md)

軌道のエネルギーと電子占有数、原子軌道と分子軌道のつながりを調べられます。

COV は軌道のエネルギーと占有数を操作可能なエネルギー準位図に表示し、図とデータを書き出せます。Gaussian の FCHK/FCH ファイルと Molden ファイルを読み込みます。

v0.4 プレビュー版では NBO 解析にも対応しました。同じ計算の NBO 出力を使って、原子軌道や局在化軌道が分子軌道にどのように寄与するかを確認し、準位図上で軌道間のつながりをたどり、電荷や結合に関する情報を調べられます。

軌道を選ぶと、その形状を 3D で表示できます。描画には NVIDIA CUDA を使用します。画面表示は英語、簡体字中国語、日本語、フランス語に対応しています。

## 主な機能

- **エネルギー準位図を読む。** 軌道のエネルギーと電子占有数を確認し、エネルギーの単位を切り替え、図を PNG または SVG で書き出せます。
- **軌道間のつながりをたどる — v0.4 プレビュー版。** 原子軌道や局在化軌道の分子軌道への寄与を表示します。つながりを選ぶと、その寄与を詳しく確認できます。
- **電荷と結合を調べる — v0.4 プレビュー版。** NBO ファイルにデータが含まれていれば、NPA 電荷、スピン分布、Wiberg 結合指数、供与体–受容体相互作用を確認できます。
- **図とデータを書き出す。** エネルギー準位図を PNG または SVG で保存し、対応するデータを CSV または JSON で書き出せます。
- **軌道を探す。** HOMO や LUMO に移動し、軌道一覧を検索できます。占有軌道、空軌道、内殻軌道、価電子軌道だけを表示することもできます。
- **軌道を 3D で見る。** 軌道を選び、等値面を調整し、分子の表示を回転・拡大縮小できます。

## ダウンロード

| バージョン | 主な内容 | Windows 用ダウンロード |
|---|---|---|
| [安定版 v0.3.0](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.3.0) | 軌道のエネルギーと占有数、エネルギー準位図、3D 表示 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.3.0/CUDA-Orbital-Visualisation-v0.3.0-Windows-sm120.zip) |
| [プレビュー版 v0.4.0-pre.1](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/tag/v0.4.0-pre.1) | NBO 解析、軌道の組成とつながりを追加 | [ZIP](https://github.com/O1dDing/Chemical-Orbital-Visualiser/releases/download/v0.4.0-pre.1/CUDA-Orbital-Visualisation-v0.4.0-pre.1-Windows-sm120.zip) |

現在配布中のファイルには旧製品名が表示されます。

現在の Windows 用ダウンロードは NVIDIA GeForce RTX 50 シリーズ GPU 向けです。ほかの NVIDIA アーキテクチャでは、対応するビルドが必要です。詳しくは[ソースからのビルド](docs/BUILD.md)をご覧ください。互換性のある NVIDIA ドライバーをインストールしてください。CUDA Toolkit が必要なのはソースからビルドする場合だけです。

## はじめに

1. Windows 用 ZIP をダウンロードして展開し、`cov.exe` を起動します。
2. Gaussian の FCHK/FCH ファイルか互換性のある Molden ファイルを開くか、ウィンドウにドラッグします。
3. 軌道一覧とエネルギー準位図でエネルギーと占有数を確認します。準位を選んで対応する軌道を表示するか、図を書き出します。
4. v0.4 プレビュー版では、計算フォルダーを開くと波動関数と NBO ファイルをまとめて読み込めます。フォルダーに複数の計算がある場合は、開く計算を選んでください。

Gaussian の `.chk` ファイルは、開く前に Gaussian の `formchk` で `.fchk` に変換してください。NBO 機能には、波動関数と同じ分子構造および電子状態に対応するファイルを使ってください。COV は計算結果を読み込むアプリであり、Gaussian や NBO の計算は実行しません。

## 入力ファイルと動作要件

- **波動関数：** Gaussian の `.fchk` / `.fch`、または互換性のある `.molden` / `.mol` / `.input` ファイル。Gaussian の `.chk` は事前に変換が必要です。
- **NBO ファイル — v0.4 プレビュー版：** レポート、アーカイブ、軌道行列はそれぞれ解析の異なる部分に必要です。各表示に必要なファイルは[NBO の結果を使う（英語）](docs/AOMO_NBO.md)をご覧ください。
- **分子の大きさ：** 入力ファイルあたり最大 100 原子。
- **グラフィックス：** NVIDIA GPU、互換性のあるドライバー、OpenGL 2.1 以降。GPU アーキテクチャに対応したビルドを使ってください。

COV は既存の計算結果を読み込みます。自動で付けられる軌道や結合のラベルには誤りがあり得るため、解釈する際は元の出力を確認してください。

## ドキュメント

- [COV の使い方](docs/UI.md)
- [NBO の結果を使う（英語）](docs/AOMO_NBO.md) — v0.4 プレビュー版
- [ソースからのビルド](docs/BUILD.md)
- [安定版のリリースノート](docs/releases/v0.3.0.md) · [プレビュー版のリリースノート](docs/releases/v0.4.0-pre.1.md)

## ライセンス

[Apache License 2.0](LICENSE)。同梱ライブラリのライセンスは[サードパーティーに関するお知らせ](THIRD_PARTY_NOTICES.md)に記載しています。

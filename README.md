Bifurcation analysis program for difference equations (descrete dynamical systems).

See aw02m/bifurcation_theory for the calcuration algorythm.

# Changelog
Aug/22/2021 : PD，G追跡に対応(initial commit)  
Aug/27/2021 : NS追跡に対応，数値微分が選択可能  
Oct/21/2021 : fixとbifの統合，json出力の自動化，数値微分が一時的に使用不可  
Oct/22/2021 : 二階変分の数値微分が復活，ppをC++(QT+QCustomPlot)に変更  
Oct/25/2021 : 変分方程式の記述を"CによるカオスCG"形式に変更(構造が`aw02m/autonomous_bif`とほぼ同一化)  
Nov/11/2021 : NS条件をBialternate Productを用いた形式に変更，虚数の取り扱いが不要のため高速化，なぜかわからんが収束速度が前の虚部分離アルゴリズムよりもかなり早い！G，PDは独立したモードで実装  
Jun/10/2022 : `bif`のMacOSでのCmakeに対応．`brew`経由で`eigen`と`nlohmann-json`を導入すれば動きます．QtのCmake対応はマダ  
Jun/21/2022 : `bif`と`pp`の`CMakeFile.txt`とツリー構造を変更．Qmakeを廃止したためすべてのプロジェクトはCmakeでコンパイル可能です．

# descrete_bif
離散力学系(差分方程式)の分岐解析ツールです．  
本プログラムは2部に分かれています:

1. pp (Phase portrait; 相平面描画ツール)
2. bif (Bifurcation; 分岐集合追跡 + 固定点追跡)

CMakeをプロジェクトのビルドに用いるため，適宜インストールをしておいてください．

## pp概要
`pp`は相平面をリアルタイムに描画します．
固定点計算のための近似値取得に使用してください．  

### 動作環境
* Qt6 : クロスプラットフォーム描画ライブラリ．Qt5, Qt4は非対応です．
* Eigen-3.4-rc1 : 線形代数ライブラリです．現在のstableに実装されていない関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`# pamac build eigen-git`で，MacOSなひとは`$brew install eigen`でインストールされます．
* nlohmann-3.10.1 : jsonライブラリ．

### pp入力ファイル要素
* "xrange", "yrange" : x,y描画区間．
* "max_plot" : 一回の描画に計算するステップ数．
* "x0" : プログラム実行時の初期値．
* "params" : 力学系のパラメタ．
* "dparams" : プログラム内でパラメタを変更する際の刻み幅．
* "explode" : 座標が発散したと判定するまでのノルム値．
* "out_path" : CSV形式での軌道の履歴出力パスを指定します．
* "json_out_path" : 出力されるjsonのパスを指定します．

### 使用法
```
cd pp
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```
* ←→ : 変更するパラメタを選択します．
* ↑↓ : 選択したパラメタを変更します．
* Space : Poincare断面上の固定点，パラメタ，周期τを出力します．
* x : x軸に対応する変数を変更します．
* y : y軸に対応する変数を変更します．
* w : 最新の状態座標及びパラメタをjsonに書き出します．
* q : プログラムを終了します．

`cmake`で`qcustomplot.h`が見つからないエラーが出たら，`build/main_autogen/ui_mainwindow.h`の`#include"../../qcustomplot.h"`を`#include"qcustomplot.h"`に書き換えてください．
`CMakeFile.txt`の更新により対応予定.

## bif概要
`bif`は分岐点の近似値をもとに分岐集合を計算します．
本プログラムはPeriod-doubling; I(PD), Tangent; G, Neimark-Sacker; NSの局所分岐計算が可能です．
Newton法の目的関数はNeimark-Sacker分岐条件ですが，I,GはNS条件に含まれます．
後述の固定点計算モードを指定することで固定点の追跡も可能です．

### 動作環境
* Eigen-3.4-rc1
* nlohmann-3.10.1

### bif入力ファイル概要
* "mode" : 0:Fixed, 1:Tangent, 2:Period-Doubling, 3:Neimark-Sacker．本バージョンのNS条件はGとPDの条件を含まないため，不正な収束が発生しにくくなりました．
* "x0" : ppにてjsonファイルを出力した際に更新されます．この値が固定点計算の初期値として使用されます．
* "period" : 固定点の周期を指定します．ppにて周期を確認して正確な整数値を与えてください．
* "inc_param" : 増加させるパラメタを指定します．"params"のインデックスで指定してください．
* "var_param" : 変数として扱うパラメタを指定します．
* "delta_inc" : パラメタ変分量です．経験的に0.1~0.0001の範囲で設定すると良いです．あまりに収束しない場合は増分パラメタの変更もしくは刻みを小さくしてください．
* "inc_iter" : 指定した回数固定点を計算します．(指定した計算回数に到達もしくは発散しないとプログラムは固定点を計算し続けます．)
* "max_iter" : Newton法の最大ステップを指定します．Newton法は通常数回の繰り返しで収束するため，10~32の整数値を指定します．
* "eps" : Newton法の収束判定を与えます．誤差が指定数値以下になった場合に計算が終了します．10^-8くらい．
* "numerical_diff" : `bool`値を指定します．`true`の場合2階変分は数値微分により計算され，`false`の場合は代数的に求められます．
* "dif_strip" : 数値微分を用いる際の微小量を指定します．`numerical_diff`がfalseの場合は無視されます．

### 使用法
```
cd bif
mkdir build
cd build
cmake ../cmake-tree
make
./main ../input/<your_input_file>.json
```
力学系の写像及びその微分は`sys_func.cpp`に記述し，それ以外のC++ファイルは変更しないでください．
固定点計算の場合は`T`と`dTdx`の定義のみ，数値微分を用いる場合は`T`, `dTdx`, `dTdlambda`の定義のみでOK．
なお，`Eigen`配列へのアクセスは`[i]`ではなく`(i)`を用いてください．こっちのほうが早いらしいです．
計算に成功すると固定点座標，パラメタ値，特性定数，特性定数のノルム・偏角が出力されます．
分岐計算の最終点はカレントディレクトリにout.jsonとして書き出されます．
収束しなかった場合は増分パラメタと変分パラメタを入れ替えるなどして追跡継続を試してください．
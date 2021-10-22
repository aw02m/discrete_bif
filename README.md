Bifurcation analysis program for difference equations (descrete dynamical systems).

See aw02m/bifurcation_theory for the calcuration algorythm.

# Changelog
Aug/22/2021 : PD，G追跡に対応(initial commit)  
Aug/27/2021 : NS追跡に対応，数値微分が選択可能  
Oct/21/2021 : fixとbifの統合，json出力の自動化，数値微分が一時的に使用不可  
Oct/22/2021 : 二階変分の数値微分が復活，ppをC++(QT+QCustomPlot)に変更 

# descrete_bif
離散力学系(差分方程式)の分岐解析ツールです．  
本プログラムは2部に分かれています:

1. pp (Phase portrait; 相平面描画ツール)
2. bif (Bifurcation; 分岐集合追跡 + 固定点追跡)

## pp概要
ppは相平面をリアルタイムに描画します．
固定点計算のための近似値取得に使用してください．  

### 動作環境
* Eigen-3.4-rc1 : 線形代数ライブラリです．現在のstableに実装されていない関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`# pamac build eigen-git`でインストールされます．
* nlohmann-3.10.1 : jsonライブラリ．
* Qt5 : クロスプラットフォーム描画ライブラリ．

### pp入力ファイル要素
* "xrange", "yrange" : x,y描画区間．
* "max_plot" : 一回の描画に計算するステップ数．2000以上は重いかも．
* "x0" : プログラム実行時の初期値．
* "params" : 力学系のパラメタ．
* "dparams" : プログラム内でパラメタを変更する際の刻み幅．
* "explode" : 座標が発散したと判定するまでのノルム値．
* "out_path" : CSV形式での軌道の履歴出力パスを指定します．
* "json_out_path" : 出力されるjsonのパスを指定します．

### 使用法
`mkdir build`で`cd build`の後`qmake ../`で`Makefile`を作成し，`make`します．`qmake`は`qt5`にバンドルされています．
* ←→ : 変更するパラメタを選択します．
* ↑↓ : 選択したパラメタを変更します．
* Space : Poincare断面上の固定点，パラメタ，周期τを出力します．
* x : x軸に対応する変数を変更します．
* y : y軸に対応する変数を変更します．
* w : 最新の状態座標及びパラメタをjsonに書き出します．
* q : プログラムを終了します．

## bif概要
`bif`は分岐点の近似値をもとに分岐集合を計算します．
本プログラムはPeriod-doubling; I(PD), Tangent; G, Neimark-Sacker; NSの局所分岐計算が可能です．
Newton法の目的関数はNeimark-Sacker分岐条件ですが，I,GはNS条件に含まれます．
後述の固定点計算モードを指定することで固定点の追跡も可能です．

### 動作環境
* Eigen-3.4-rc1 : 線形代数ライブラリです．現在のstableに実装されていない関数を使用しているため，gitリポジトリの最新バージョンを利用してください．Arch Linuxなひとは`# pamac build eigen-git`でインストールされます．
* nlohmann-3.10.1 : jsonライブラリ．
* cmake : MakeFileの自動生成に使います．

### bif入力ファイル概要
* "fix_mode" : trueの場合固定点追跡を行います．
* "x0" : ppにてjsonファイルを出力した際に更新されます．この値が固定点計算の初期値として使用されます．
* "theta" : 分岐点特性定数の偏角を与えます．
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
`./main [input json file]`にて実行．コンパイルは`mkdir build`の後`bif/build`ディレクトリの中で`cmake ../`するとMakefileが自動で作成されますので，その後`make`で実行ファイルが生成されます．
力学系の写像及びその微分は`sys_func.cpp`に記述し，それ以外のC++ファイルは変更しないでください．
固定点計算の場合は`T`と`dTdx`の定義のみ，数値微分を用いる場合は`T`, `dTdx`, `dTdlambda`の定義のみでOK．
なお，`Eigen`配列へのアクセスは`[i]`ではなく`(i)`を用いてください．こっちのほうが早いらしいです．
計算に成功すると固定点座標，パラメタ値，特性定数，特性定数のノルム・偏角が出力されます．
分岐計算の最終点はカレントディレクトリにout.jsonとして書き出されます．
収束しなかった場合は増分パラメタと変分パラメタを入れ替えるなどして追跡継続を試してください．

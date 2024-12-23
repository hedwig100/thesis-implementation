# Thesis

このレポジトリは「効率的なバックトラック判定による同種写像ハッシュ関数の高速化」に関するプログラムがおいてあるレポジトリです.

## CGLハッシュ関数の計算コストの評価

後にあるように必要なパッケージをインストールした後, `./do run count-cgl`としてください. 

NOTE: DPB-CGLは衝突耐性を**持っていません**, 比較のために載せているだけです.

```
...
Consolidate compiler generated dependencies of target count-cgl
[100%] Built target count-cgl
+ ./build/src/count-cgl
Size of p (bit):256
n:240
P: 94482146789024102024466833852227567633874701385423723785552891614122842521599
e: 240
  name: ModularPolynomial-CGL
  Count:   M_COUNTER: 279.57
  S_COUNTER: 758.615
  A_COUNTER: 24.445
  I_COUNTER: 0.985
  name: Yoshida-Takashima-CGL
  Count:   M_COUNTER: 273.095
  S_COUNTER: 762.45
  A_COUNTER: 12.095
  I_COUNTER: 0.995
  name: Hashimoto-Nuida-CGL
  Count:   M_COUNTER: 383.41
  S_COUNTER: 510.88
  A_COUNTER: 7.835
  I_COUNTER: 0.505
  name: Radical-CGL
  Count:   M_COUNTER: 388.73
  S_COUNTER: 510.925
  A_COUNTER: 9.165
  I_COUNTER: 0.505
  name: DPB-CGL
  Count:   M_COUNTER: 165.076
  S_COUNTER: 3.43542
  A_COUNTER: 176.323
  I_COUNTER: 0.004375
  name: DPB-CGL-prevent-backtrack
  Count:   M_COUNTER: 205.565
  S_COUNTER: 6.41
  A_COUNTER: 215.793
  I_COUNTER: 0.00833333
  name: My-CGL
  Count:   M_COUNTER: 199.82
  S_COUNTER: 5.24104
  A_COUNTER: 211.639
  I_COUNTER: 0.00666667
...
```

## バックトラック判定をしていない場合に衝突を見つける.

バックトラック判定を行わない場合に衝突を発券する.
`./do run find_collision`を実行してください.
これは失敗することがありますが, 成功した場合は衝突を見つけられています.

```
...
Consolidate compiler generated dependencies of target find_collision
[100%] Built target find_collision
+ ./build/src/find_collision
Found collision!!
P:112200089154621740081528141189677554510592709903292430697596358600685233635327
n: 240
m0: 101001001100110110001110100011000001110100111100111101000111000001111011101000000111101100110110011101001000101101100101011011111111000111101010010011111100010010010001001111011000101111011011111101111011000101110101101100110001110011110010
m0': 001110010110110001000110100101011111101110011100110010000011011111101110001101000000110111101010111010011010000010010100011011111111000111101010010011111100010010010001001111011000101111011011111101111011000101110101101100110001110011110010
m1': 110111000011011001111001011010110110101011110010000111101011011010000101010000111110001111001000111000101110000000110100011101000101100110100011110001100111000110111001111111110000101001101111101000010101011111001100110001101011110001101111
```

## Requirements [C++]

多倍長整数の計算ライブラリとして`libboost-dev`が必要です.

```
sudo apt-get install libboost-dev
```

## 実行 [C++]

まずビルドする

```
./do build
```

次にテストを走らせる.

```
./do test
```

最後に`src/<exec>.cpp`を実行する.

```
./do run <exec>
```

## Requirements [Python]

poetryをインストールする.

## 実行 [Python]

まず, ライブラリをインストールする.

```
poetry install
```

つぎにスクリプトを動かす.

```
./do python strategy_optimize
```

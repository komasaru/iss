ビルド方法
==========

`make`

（やり直す場合は、 `make clean` をしてから）

準備
====

* `tle_iss_nasa.txt`, `Leap_Second.dat`, `eop.csv` を最新のものにしておく。
  + `tle_iss_nasa.txt` の取得は「[Python - TLE（2行軌道要素形式）の取得(NASA)！](https://www.mk-mode.com/octopress/2018/08/14/python-tle-getting-from-nasa "Python - TLE（2行軌道要素形式）の取得(NASA)！")」を参考に。
  + `Leap_Second.dat` は「[こちら](https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat "Leap_Second.dat")」を取得したものをそのまま配置する。
  + `eop.csv` の取得は「[Ruby, Python - EOP（地球姿勢パラメータ）CSV 生成！](https://www.mk-mode.com/octopress/2018/08/29/ruby-python-eop-getting-from-iers "Ruby, Python - EOP（地球姿勢パラメータ）CSV 生成！")」
* `Leap_Second.dat` は `file` ディレクトリ配下に、 `tle_iss_nasa.txt`, `eop.csv` は `data` ディレクトリ配下に配置する。

実行方法
========

`./iss_sgp4_json [YYYYMMDD[HHMMSS[MMM]]]`

* JST（日本標準時）は「年・月・日・時・分・秒・ミリ秒」を最大17桁で指定する。
* JST（日本標準時）を指定しない場合は、システム日時を JST とみなす。
* 実行が正常に終了すれば、 `data` ディレクトリ配下に `iss.json` が作成される。


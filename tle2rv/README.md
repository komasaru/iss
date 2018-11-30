ビルド方法
==========

`make`

（やり直す場合は、 `make clean` をしてから）

準備
====

* `tle_iss_nasa.txt` を最新のものにしておく。
  （`tle_iss_nasa.txt` の取得は「[Python - TLE（2行軌道要素形式）の取得(NASA)！](https://www.mk-mode.com/octopress/2018/08/14/python-tle-getting-from-nasa "Python - TLE（2行軌道要素形式）の取得(NASA)！")」を参考に）

実行方法
========

`./tle2rv [YYYYMMDD[HHMMSS[MMM]]]`

* UT1（世界時1）は「年・月・日・時・分・秒・ミリ秒」を最大17桁で指定する。
* UT1（世界時1）を指定しない場合は、システム日時を UT1 とみなす。


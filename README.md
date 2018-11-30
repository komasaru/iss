プログラム一覧
==============

それぞれのディレクトリ内で完結するようにしている。（モジュールの共有化は図っていない）  
（詳細説明はそれぞれのディレクトリ内の `README.md` を参照）

* `tle2rv`  
  TLE から任意の時刻(UT1)の ISS の位置・速度を計算する。  
  （但し、座標系は TEME(True Equator, Mean Equinox; 真赤道面平均春分点）
* `iss_sgp4`  
  TLE から任意の時刻(JST)の ISS の位置・速度を計算する。  
  （但し、座標系は ECEF: Earth Centered, Earth Fixed; 地球中心・地球固定直交座標系; いわゆる、緯度／経度／高度）
* `iss_sgp4_json`  
  TLE から任意の時刻(JST)から10秒間隔2日分の ISS の位置・速度を計算を計算し、 JSON 出力する。  
  （但し、座標系は ECEF: Earth Centered, Earth Fixed; 地球中心・地球固定直交座標系; いわゆる、緯度／経度／高度）


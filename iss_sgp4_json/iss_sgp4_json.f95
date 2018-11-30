!*******************************************************************************
! Getting of ISS position.
!  : 指定日時から48時間分の ISS 位置／速度を計算し、 JSON 出力する。
!
!   date          name            version
!   2018.11.28    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
! ---
! 引数 : JST（日本標準時）
!          YYYYMMDD[HHMMSS[MMM]]
!          無指定なら現在(システム日時)を JST とみなす。
!          （上記最後の MMM はミリ秒）
! ---
! MEMO:
!   TEME: True Equator, Mean Equinox; 真赤道面平均春分点
!    PEF: Pseudo Earth Fixed; 擬地球固定座標系
!   ECEF: Earth Centered, Earth Fixed; 地球中心・地球固定直交座標系
!*******************************************************************************
!
program iss_sgp4_json
  use const,       only : SP, DP, JST_UTC, UID, F_TLE, SEC_D, TT_TAI, &
                        & DAY_JC, UID_J, F_JSON, FMT_DT_0, FMT_DT_1, &
                        & FMT_DT_3, getgravconst
  use propagation, only : propagate
  use io
  use model
  use ext
  use eph
  implicit none
  integer(SP), parameter :: DAY = 2   ! 計算日数(日)
  integer(SP), parameter :: SEC = 10  ! 計算間隔(秒)
  type(t_time)  :: jst, jst_w, utc, utc_w, ut1, tai, tt
  type(t_eop)   :: eop(DAY + 1)
  type(t_sat)   :: satrec
  type(t_cst)   :: gravconst
  character(69) :: tle(2)
  character(10) :: date_w, date_sav
  character(19) :: date_jst, date_utc
  integer(SP)   :: dat, i, ios, idx
  real(DP)      :: jd_ut1, jd_tt, t_tt, gmst, om, gmst_g, r(3), v(3)
  real(DP)      :: mtx_z(3, 3), mtx_pm(3, 3)
  real(DP)      :: r_pef(3), r_ecef(3), om_e(3), v_pef(3), v_ecef(3)
  real(DP)      :: lat, lon, ht, vel_ecef

  ! コマンドライン引数（現在日時(JST)）取得
  call get_arg(jst)
  if (jst%year == 0) stop

  ! 各種時刻換算
  call add_day(jst, -JST_UTC / 24.0_DP, utc)
  do i = 1, DAY + 1
    call add_day(utc, real(i - 1, DP), utc_w)
    call get_eop(utc_w, eop(i))
  end do
  call get_dat(utc, dat)

  ! TLE 読み込み
  call get_tle(utc, tle)

  ! 書き込み用 JSON ファイル OPEN
  open (unit   = UID_J,       &
      & iostat = ios,         &
      & file   = F_JSON,      &
      & action = "write",     &
      & form   = "formatted", &
      !& status = "new")
      & status = "replace")
  if (ios /= 0) then
    print '("[ERROR:", I0 ,"] Failed to open file: ", A)', ios, F_JSON
    stop
  end if

  ! 主処理
  idx = 0
  write (UID_J, '(A)') "{"
  write (UID_J, '("  ""tle1"": """, A, """,")') tle(1)
  write (UID_J, '("  ""tle2"": """, A, """,")') tle(2)
  write (UID_J, '("  ""counts"": ", I0, ",")') int(SEC_D * DAY / SEC)
  write (UID_J, '(A)') "  ""data"": ["
  do i = 0, int(SEC_D * DAY / SEC) - 1  ! 10秒間隔
    call add_day(jst, i * SEC / SEC_D, jst_w)
    call add_day(jst_w, -JST_UTC / 24.0_DP, utc_w)
    write (date_w, '(I4, "-", I0.2, "-", I0.2)') &
      & utc_w%year, utc_w%month, utc_w%day
    write (date_jst, FMT_DT_3) &
      & jst_w%year, jst_w%month, jst_w%day, &
      & jst_w%hour, jst_w%minute, jst_w%second
    write (date_utc, FMT_DT_3) &
      & utc_w%year, utc_w%month, utc_w%day, &
      & utc_w%hour, utc_w%minute, utc_w%second
    if (date_w /= date_sav) idx = idx + 1

    call add_day(utc_w, eop(1)%dut1 / SEC_D, ut1)
    call get_dat(utc_w, dat)
    call add_day(utc_w, dat / SEC_D, tai)
    call add_day(tai, TT_TAI / SEC_D, tt)
    call gc2jd(ut1, jd_ut1)
    call gc2jd(tt, jd_tt)
    t_tt = (jd_tt - 2451545.0_DP) / DAY_JC

    ! WGS84 用定数の取得
    gravconst = getgravconst("wgs84")
    ! ISS 初期位置・速度の取得
    call twoline2rv(tle, gravconst, satrec)
    ! 指定 UT1 の ISS 位置・速度の取得
    call propagate(gravconst, satrec, ut1, r, v)
    ! Error 時、終了
    !   1 - mean elements, ecc >= 1.0 or ecc < -0.001 or a < 0.95 er
    !   2 - mean motion less than 0.0
    !   3 - pert elements, ecc < 0.0  or  ecc > 1.0
    !   4 - semi-latus rectum < 0.0
    !   5 - epoch elements are sub-orbital
    !   6 - satellite has decayed
    if (satrec%error /= 0) then
      print '(A, I0, A)', "ERROR! [", satrec%error, "]"
      stop
    end if

    ! GMST（グリニッジ平均恒星時）計算
    gmst = calc_gmst(jd_ut1)
    ! Ω（月の平均昇交点黄経）計算（IAU1980章動理論）
    om = calc_om(t_tt)
    ! GMST に運動項を適用（1997年より新しい場合）
    gmst_g = apply_kinematic(gmst, jd_ut1, om)
    ! GMST 回転行列（z軸を中心とした回転）
    mtx_z = r_z(gmst_g)
    ! 極運動(Polar Motion)回転行列
    mtx_pm = r_pm(eop(idx)%pm_x, eop(idx)%pm_y, t_tt)
    ! PEF 座標の計算（GMST 回転行列の適用）
    r_pef = matmul(r, mtx_z)
    ! ECEF 座標（位置）の計算（極運動(Polar Motion)の適用）
    r_ecef = matmul(r_pef, mtx_pm)
    ! Ω_earth値の計算
    om_e = calc_om_e(eop(idx)%lod)
    ! PEF 座標（速度）の計算（GMST 回転行列の適用）
    v_pef = matmul(v, mtx_z) - v_cross(om_e, r_pef)
    ! ECEF 座標（速度）の計算（極運動(Polar Motion)の適用）
    v_ecef = matmul(v_pef, mtx_pm)
    ! ECEF 座標 => BLH(Beta, Lambda, Height) 変換
    call ecef2blh(r_ecef, lat, lon, ht)
    vel_ecef = sqrt(sum(v_ecef * v_ecef))

    write (UID_J, '(A)') "    {"
    write (UID_J, '("      ""jst"": """,     A,    """,")') date_jst
    write (UID_J, '("      ""utc"": """,     A,    """,")') date_utc
    write (UID_J, '("      ""latitude"": ",  F17.12, ",")')      lat
    write (UID_J, '("      ""longitude"": ", F17.12, ",")')      lon
    write (UID_J, '("      ""height"": ",    F17.12, ",")')       ht
    write (UID_J, '("      ""velocity"": ",  F14.12)'     ) vel_ecef
    if (i == int(SEC_D * DAY / 10) - 1) then
      write (UID_J, '(A)') "    }"
    else
      write (UID_J, '(A)') "    },"
    end if

    date_sav = date_w
  end do
  write (UID_J, '(A)') "  ]"
  write (UID_J, '(A)') "}"

  ! JSON ファイル CLOSE
  close(UID_J)

  stop
contains
  ! コマンドライン引数取得
  ! * YYYYMMDD[HHMMSS[MMM]] 形式
  ! * 17桁超入力された場合は、18桁目以降の部分は切り捨てる
  ! * コマンドライン引数がなければ、システム日時を JST とする
  ! * 日時の整合性チェックは行わない
  !
  ! :param(out) type(t_time) jst
  subroutine get_arg(jst)
    implicit none
    type(t_time), intent(out) :: jst
    character(17) :: gc
    integer(SP)   :: dt(8)
    integer(SP)   :: len_gc

    if (iargc() == 0) then
      call date_and_time(values=dt)
      jst = t_time(dt(1), dt(2), dt(3), dt(5), dt(6), dt(7), dt(8))
    else
      call getarg(1, gc)
      len_gc = len(trim(gc))
      if (len_gc /= 8 .and. len_gc /= 14 .and. len_gc /= 17) then
        print *, "Format: YYYYMMDD[HHMMSS[MMM]"
        return
      end if
      read (gc, FMT_DT_0) jst
    end if
  end subroutine get_arg

  ! TLE 取得
  !
  ! :param(in)  type(t_time)     ut1: UT1
  ! :param(out) character(69) tle(2): TLE(2行)
  subroutine get_tle(ut1, tle)
    implicit none
    type(t_time),  intent(in)  :: ut1
    character(69), intent(out) :: tle(2)
    character(69) :: buf, buf_p(2)
    integer(SP)   :: ios
    type(t_time)  :: utc  ! TLE 内の日時
    integer(SP)   :: i, l, y
    real(DP)      :: d
    character(17) :: s_ut1, s_utc

    write (s_ut1, FMT_DT_1) ut1

    ! ファイル OPEN
    open (unit   = UID,         &
        & iostat = ios,         &
        & file   = F_TLE,       &
        & action = "read",      &
        & form   = "formatted", &
        & status = "old")
    if (ios /= 0) then
      print '("[ERROR:", I0 ,"] Failed to open file: ", A)', ios, F_TLE
      stop
    end if

    buf_p = ""
    do
      read (UID, '(A)', iostat = ios) buf
      if (ios < 0) then
        exit
      else if (ios > 0) then
        print '("[ERROR:", I0 ,"] Failed to read file: ", A)', ios, F_TLE
      end if
      if (len(trim(buf)) == 0) cycle
      if (buf(1:1) /= "1" .and. buf(1:1) /= "2") cycle
      if (buf(1:1) == "1") then
        read (buf(19:20), '(I2)') y
        y = 2000 + y
        read (buf(21:32), '(F12.8)') d
        call add_day(t_time(y, 1, 1, 0, 0, 0, 0), d, utc)
        write (s_utc, FMT_DT_1) utc
        if (s_utc > s_ut1) then
          tle = buf_p
          exit
        end if
      end if
      read (buf(1:1), '(I1)') l
      buf_p(l) = buf
    end do

    ! 上記の処理で該当レコードが得られなかった場合は、最初の2行
    if (len(trim(tle(1))) == 0) then
      rewind(UID)
      do i = 1, 2
        read (UID, '(A)', iostat = ios) tle(i)
      end do
    end if

    ! ファイル CLOSE
    close(UID)
  end subroutine get_tle
end program iss_sgp4_json


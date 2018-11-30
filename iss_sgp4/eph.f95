!*******************************************************************************
! Modules about ephemeris
!
!   date          name            version
!   2018.11.27    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module eph
  use const, only : SP, DP, DAY_JC, PI, PI2, PI_180, SEC_D, &
                  & WGS84_A, WGS84_B, WGS84_E2, WGS84_ED2
  implicit none
  private
  public :: calc_gmst, calc_om, apply_kinematic, r_z, r_pm, &
          & calc_om_e, ecef2blh
contains
  ! =============================
  ! Private subroutines/functions
  ! =============================

  ! ============================
  ! Public subroutines/functions
  ! ============================

  ! GMST（グリニッジ平均恒星時）計算
  ! * IAU1982理論(by David Vallado)によるもの
  !     GMST = 18h 41m 50.54841s
  !          + 8640184.812866s T + 0.093104s T^2 - 0.0000062s T^3 
  !     (但し、 T = 2000年1月1日12時(UT1)からのユリウス世紀単位)
  !
  ! :param(in) real(8) jd_ut1: Julian Day (UT1)
  ! :return    real(8)   gmst: グリニッジ平均恒星時(単位:radian)
  function calc_gmst(jd_ut1) result (gmst)
    implicit none
    real(DP), intent(in) :: jd_ut1
    real(DP) :: gmst
    real(DP) :: t_ut1

    t_ut1 = (jd_ut1 - 2451545.0_DP) / DAY_JC
    gmst =  67310.54841_DP &
       & + (876600.0_DP * 3600.0_DP + 8640184.812866_DP &
       & + (0.093104_DP    &
       & -  6.2e-6_DP * t_ut1) * t_ut1) * t_ut1
    gmst = mod(gmst * PI_180 / 240.0_DP, PI2)
    if (gmst < 0.0_DP) gmst = gmst + PI2
  end function calc_gmst

  ! Ω（月の平均昇交点黄経）計算（IAU1980章動理論）
  ! * Ω = 125°02′40″.280
  !      - ((5 * 360)° + 134°08′10″.539) * T
  !      + 7″.455 * T2
  !      + 0″.008 * T3
  !
  ! :param(in) real(8) t_tt: ユリウス世紀数(TT)
  ! :return    real(8)   om: Ω（月の平均昇交点黄経）
  function calc_om(t_tt) result(om)
    implicit none
    real(DP), intent(in) :: t_tt
    real(DP) :: om

    om =  125.04452222_DP   &
     & + ((-6962890.5390_DP &
     & +  (7.455_DP         &
     & +   0.008_DP * t_tt) * t_tt) * t_tt) / 3600.0_DP
    om = mod(om, 360.0_DP)
    if (om < 0.0_DP) om = om + 360.0_DP
    om = om * PI_180
  end function calc_om

  ! GMST に運動項を適用（1997年より新しい場合）
  ! * gmst_g = gmst \
  !          + 0.00264 * PI / (3600 * 180) * sin(om) \
  !          + 0.000063 * PI / (3600 * 180) * sin(2.0 * om)
  !
  ! :param(in) real(8)   gmst: GMST（グリニッジ平均恒星時）（適用前）
  ! :param(in) real(8) jd_ut1: ユリウス日(UT1)
  ! :param(in) real(8)     om: Ω（月の平均昇交点黄経）
  ! :return    real(8) gmst_g: GMST（グリニッジ平均恒星時）（適用後）
  function apply_kinematic(gmst, jd_ut1, om) result(gmst_g)
    implicit none
    real(DP), intent(in) :: gmst, jd_ut1, om
    real(DP) :: gmst_g

    if (jd_ut1 > 2450449.5_DP) then
      gmst_g = gmst &
           & + 0.00264_DP * PI / (3600 * 180) * sin(om) &
           & + 0.000063_DP * PI / (3600 * 180) * sin(2 * om)
    else
      gmst_g = gmst
    end if
    gmst_g = mod(gmst_g, PI2)
  end function apply_kinematic

  ! 回転行列生成(z軸中心)
  !   ( cos(θ)  -sin(θ)  0 )
  !   ( sin(θ)  +cos(θ)  0 )
  !   (    0         0     1 )
  !
  ! :param(in)  real(8)       theta: 角度(Unit: rad)
  ! :return     real(8) r_mtx(3, 3): 回転行列
  function r_z(theta) result(r_mtx)
    implicit none
    real(DP), intent(in) :: theta
    real(DP) :: r_mtx(3, 3)
    real(DP) :: s, c

    s = sin(theta)
    c = cos(theta)
    r_mtx = reshape((/ &
      &      c,      s, 0.0_DP, &
      &     -s,      c, 0.0_DP, &
      & 0.0_DP, 0.0_DP, 1.0_DP  &
    /), (/3, 3/))
  end function r_z

  ! 極運動(Polar Motion)回転行列
  !
  ! :param(in) real(8)        pm_x: Polar Motion X
  ! :param(in) real(8)        pm_y: Polar Motion Y
  ! :param(in) real(8)        t_tt: ユリウス世紀数(TT)
  ! :return    real(8) r_mtx(3, 3): 回転行列
  function r_pm(pm_x, pm_y, t_tt) result(r_mtx)
    implicit none
    real(DP), intent(in) :: pm_x, pm_y, t_tt
    real(DP) :: r_mtx(3, 3)
    real(DP) :: pm_x_r, pm_y_r
    real(DP) :: conv, c_xp, s_xp, c_yp, s_yp, sp, s_sp, c_sp

    pm_x_r = pm_x * PI / (180 * 60 * 60 * 1000)
    pm_y_r = pm_y * PI / (180 * 60 * 60 * 1000)
    conv = PI / (3600 * 180)
    c_xp = cos(pm_x_r)
    s_xp = sin(pm_x_r)
    c_yp = cos(pm_y_r)
    s_yp = sin(pm_y_r)
    ! approximate sp value in rad
    sp = -47.0e-6_DP * t_tt * conv
    s_sp = sin(sp)
    c_sp = cos(sp)
    r_mtx = reshape((/ &
      &  c_xp * c_sp,  c_xp * s_sp,  s_xp, &
      & -c_yp * s_sp + s_yp * s_xp * c_sp, &
      &  c_yp * c_sp + s_yp * s_xp * s_sp, &
      & -s_yp * c_xp,                      &
      & -s_yp * s_sp - c_yp * s_xp * c_sp, &
      &  s_yp * c_sp - c_yp * s_xp * s_sp, &
      &  c_yp * c_xp                       &
    /), (/3, 3/))
  end function r_pm

  ! Ω_earth値の計算
  !
  ! :param(in) real(8)  lod: Length Of Day
  ! :return    real(8) om_e: Ω_earth
  function calc_om_e(lod) result(om_e)
    implicit none
    real(DP), intent(in) :: lod
    real(DP) :: om_e(3)
    real(DP) :: theta_sa

    theta_sa = 7.29211514670698e-05_DP * (1.0_DP - lod / SEC_D)
    om_e = (/0.0_DP, 0.0_DP, theta_sa/)
  end function calc_om_e

  ! ECEF 座標 => BLH(Beta, Lambda, Height) 変換
  !      β = atan {(z + e'^2 * b * sin^3(θ)) / (p − e^2 * a * cos^3(θ))}
  !      λ = atan (y / x)
  !      h  = (p / cos(β)) − N
  !   但し、
  !       p = sqrt(x^2 + y^2)
  !      θ = atan(za / pb)
  !     e^2 = (a^2 - b^2) / a^2
  !    e'^2 = (a^2 - b^2) / b^2
  !
  ! :param(in)  real(8) r_ecef(3): ECEF 座標 (unit: km)
  ! :param(out) real(8)       lat: BLH (latitude) (°)
  ! :param(out) real(8)       lon: BLH (longitude)(°)
  ! :param(out) real(8)        ht: BLH (height)   (km)
  subroutine ecef2blh(r_ecef, lat, lon, ht)
    implicit none
    real(DP), intent(in)  :: r_ecef(3)
    real(DP), intent(out) :: lat, lon, ht
    real(DP) :: x, y, z, p, theta

    x = r_ecef(1) * 1.0e3_DP
    y = r_ecef(2) * 1.0e3_DP
    z = r_ecef(3) * 1.0e3_DP
    p = sqrt(x * x + y * y)
    theta = atan2(z * WGS84_A, p * WGS84_B) / PI_180
    lat = atan2( &
      & z + WGS84_ED2 * WGS84_B * sin(theta * PI_180)**3, &
      & p - WGS84_E2 * WGS84_A * cos(theta * PI_180)**3   &
    & ) / PI_180
    lon = atan2(y, x) / PI_180
    ht = p / cos(lat * PI_180) &
     & - WGS84_A / sqrt(1 - WGS84_E2 * sin(lat * PI_180)**2)
    ht = ht / 1.0e3_DP
  end subroutine ecef2blh
end module eph


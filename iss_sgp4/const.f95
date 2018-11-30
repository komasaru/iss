!*******************************************************************************
! 定数モジュール
!
!   date          name            version
!   2018.11.21    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module const
  implicit none

  ! SP: 単精度(4), DP: 倍精度(8)
  integer,      parameter :: SP = kind(1.0)
  integer(SP),  parameter :: DP = selected_real_kind(2 * precision(1.0_SP))
  integer(SP),  parameter :: UID     = 10                     ! TLE 用 Unit ID
  integer(SP),  parameter :: UID_E   = 11                     ! EOP 用 Unit ID
  integer(SP),  parameter :: UID_D   = 12                     ! DAT 用 Unit ID
  character(*), parameter :: F_TLE   = "tle_iss_nasa.txt"     ! TLE ファイル
  character(*), parameter :: F_EOP   = "eop.csv"              ! EOP ファイル
  character(*), parameter :: F_DAT   = "Leap_Second.dat"      ! DAT ファイル
  integer(SP),  parameter :: JST_UTC = 9                      ! JST - UTC (hour)
  integer(SP),  parameter :: DAY_JC  = 36525                  ! Days per Julian century
  real(DP),     parameter :: PI      = atan(1.0_DP) * 4.0_DP  ! 円周率
  real(DP),     parameter :: PI2     = PI * 2.0_DP            ! PI * 2
  real(DP),     parameter :: PI_180  = PI / 180.0_DP          ! PI / 180
  real(DP),     parameter :: DEG2RAD = PI / 180.0_DP          ! 0.0174532925199433
  real(DP),     parameter :: MIN_D   = 1440.0_DP              ! Minutes per day
  real(DP),     parameter :: SEC_D   = 86400.0_DP             ! Seconds per day
  real(DP),     parameter :: TT_TAI  = 32.184_DP              ! TT - TAI
  real(DP),     parameter :: XPDOTP  = MIN_D / (2.0_DP * PI)  ! 229.1831180523293
  integer(SP),  parameter :: DAYS(1:12) = &
   & (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)       ! Days per month
  ! WGS84 座標パラメータ
  real(DP), parameter :: WGS84_A     = 6378137.0_DP           ! a(地球楕円体長半径(赤道面平均半径))
  real(DP), parameter :: WGS84_ONE_F = 298.257223563_DP       ! 1 / f(地球楕円体扁平率=(a - b) / a)
  real(DP), parameter :: WGS84_B     = WGS84_A * (1.0_DP - 1.0_DP / WGS84_ONE_F)
                                                              ! b(地球楕円体短半径)
  real(DP), parameter :: WGS84_E2    = (1.0_DP / WGS84_ONE_F) &
                                   & * (2.0_DP - (1.0_DP / WGS84_ONE_F))
                                                              ! e^2 = 2 * f - f * f
                                                              !     = (a^2 - b^2) / a^2
  real(DP), parameter :: WGS84_ED2   = WGS84_E2 * WGS84_A * WGS84_A &
                                   & / (WGS84_B * WGS84_B)    ! e'^2= (a^2 - b^2) / b^2
  ! 各種フォーマット
  character(*), parameter :: FMT_DT_0 = '(I4I2I2I2I2I2I3)'  ! 日時取得用
  character(*), parameter :: FMT_DT_1 = &
    & '(I4, I0.2, I0.2, I0.2, I0.2, I0.2, I0.3)'
  character(*), parameter :: FMT_DT_2 = &
    & '(I4, "-", I0.2, "-", I0.2, " ", I0.2, ":", I0.2, ":", I0.2, ".", I0.3)'
  ! 構造型
  type, public :: t_cst
    real(DP) :: mu            = 0.0_DP
    real(DP) :: radiusearthkm = 0.0_DP
    real(DP) :: xke           = 0.0_DP
    real(DP) :: tumin         = 0.0_DP
    real(DP) :: j2            = 0.0_DP
    real(DP) :: j3            = 0.0_DP
    real(DP) :: j4            = 0.0_DP
    real(DP) :: j3oj2         = 0.0_DP
  end type t_cst
contains
  ! function getgravconst
  ! * this function gets constants for the propagator. note that mu is identified to
  !   facilitiate comparisons with newer models. the common useage is wgs72.
  !
  ! :param(in) character(*) whichconst: which set of constants to use
  !                                      wgs72old, wgs72, wgs84
  ! :return    type(t_cst)   gravconst: tumin      - minutes in one time unit
  !                                     mu         - earth gravitational parameter
  !                                     radiusearthkm - radius of the earth in km
  !                                     xke        - reciprocal of tumin
  !                                     j2, j3, j4 - un-normalized zonal harmonic values
  !                                     j3oj2      - j3 divided by j2
  function getgravconst(whichconst) result(gravconst)
    implicit none
    character(*), intent(in) :: whichconst
    type(t_cst) :: gravconst
    real(DP)    :: mu, radiusearthkm, xke, tumin, j2, j3, j4, j3oj2

    select case (trim(whichconst))
    case ("wgs72old")
      radiusearthkm = 6378.135_DP  ! km
      mu    = 398600.79964_DP      ! in km3 / s2
      xke   =  0.0743669161_DP
      tumin =  1.0_DP / xke
      j2    =  0.001082616_DP
      j3    = -0.00000253881_DP
      j4    = -0.00000165597_DP
      j3oj2 = j3 / j2
    case ("wgs72")
      radiusearthkm = 6378.135_DP  ! km
      mu    = 398600.8_DP          ! in km3 / s2
      xke   = 60.0_DP / sqrt( &
            & radiusearthkm * radiusearthkm * radiusearthkm / mu)
      tumin =  1.0_DP / xke
      j2    =  0.001082616_DP
      j3    = -0.00000253881_DP
      j4    = -0.00000165597_DP
      j3oj2 = j3 / j2
    case ("wgs84")
      radiusearthkm = 6378.137_DP  ! km
      mu    = 398600.5_DP          ! in km3 / s2
      xke   = 60.0_DP / sqrt( &
            & radiusearthkm * radiusearthkm * radiusearthkm / mu)
      tumin =  1.0_DP / xke
      j2    =  0.00108262998905_DP
      j3    = -0.00000253215306_DP
      j4    = -0.00000161098761_DP
      j3oj2 = j3 / j2
    end select
    gravconst = t_cst(mu, radiusearthkm, xke, tumin, j2, j3, j4, j3oj2)
  end function getgravconst
end module const


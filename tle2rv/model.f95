!*******************************************************************************
! Modules for a satellite model
! (An earth-orbiting satellite as represented by the SGP4 model)
!
!   date          name            version
!   2018.11.22    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module model
  use const, only : SP, DP, MIN_D, t_cst
  use ext,   only : t_time, jday
  implicit none
  type, public :: t_sat
    ! Unique satellite number given in the TLE file.
    integer(SP)  :: satnum     = 0
    ! Full four-digit year of this element set's epoch moment.
    integer(SP)  :: epochyr    = 0
    ! Fractional days into the year of the epoch moment.
    real(DP)     :: epochdays  = 0.0_DP
    ! Julian date of the epoch (computed from epochyr and epochdays).
    real(DP)     :: jdsatepoch = 0.0_DP
    ! First time derivative of the mean motion (ignored by SGP4).
    real(DP)     :: ndot       = 0.0_DP
    ! Second time derivative of the mean motion (ignored by SGP4).
    real(DP)     :: nddot      = 0.0_DP
    ! Ballistic drag coefficient B* in inverse earth radii.
    real(DP)     :: bstar      = 0.0_DP
    ! Inclination in radians.
    real(DP)     :: inclo      = 0.0_DP
    ! Right ascension of ascending node in radians.
    real(DP)     :: nodeo      = 0.0_DP
    ! Eccentricity.
    real(DP)     :: ecco       = 0.0_DP
    ! Argument of perigee in radians.
    real(DP)     :: argpo      = 0.0_DP
    ! Mean anomaly in radians.
    real(DP)     :: mo         = 0.0_DP
    ! Mean motion in radians per minute.
    real(DP)     :: no         = 0.0_DP
    !
    ! Near Earth
    character(1) :: method  = "n"
    character(1) :: opsmode = "i"
    character(1) :: init    = "y"
    integer(SP)  :: isimp   = 0
    real(DP)     :: aycof   = 0.0_DP
    real(DP)     :: con41   = 0.0_DP
    real(DP)     :: cc1     = 0.0_DP
    real(DP)     :: cc4     = 0.0_DP
    real(DP)     :: cc5     = 0.0_DP
    real(DP)     :: d2      = 0.0_DP
    real(DP)     :: d3      = 0.0_DP
    real(DP)     :: d4      = 0.0_DP
    real(DP)     :: delmo   = 0.0_DP
    real(DP)     :: eta     = 0.0_DP
    real(DP)     :: argpdot = 0.0_DP
    real(DP)     :: omgcof  = 0.0_DP
    real(DP)     :: sinmao  = 0.0_DP
    real(DP)     :: t       = 0.0_DP
    real(DP)     :: t2cof   = 0.0_DP
    real(DP)     :: t3cof   = 0.0_DP
    real(DP)     :: t4cof   = 0.0_DP
    real(DP)     :: t5cof   = 0.0_DP
    real(DP)     :: x1mth2  = 0.0_DP
    real(DP)     :: x7thm1  = 0.0_DP
    real(DP)     :: mdot    = 0.0_DP
    real(DP)     :: nodedot = 0.0_DP
    real(DP)     :: xlcof   = 0.0_DP
    real(DP)     :: xmcof   = 0.0_DP
    real(DP)     :: nodecf  = 0.0_DP
    !
    ! Deep space
    integer(SP)  :: irez  = 0
    real(DP)     :: d2201 = 0.0_DP
    real(DP)     :: d2211 = 0.0_DP
    real(DP)     :: d3210 = 0.0_DP
    real(DP)     :: d3222 = 0.0_DP
    real(DP)     :: d4410 = 0.0_DP
    real(DP)     :: d4422 = 0.0_DP
    real(DP)     :: d5220 = 0.0_DP
    real(DP)     :: d5232 = 0.0_DP
    real(DP)     :: d5421 = 0.0_DP
    real(DP)     :: d5433 = 0.0_DP
    real(DP)     :: dedt  = 0.0_DP
    real(DP)     :: del1  = 0.0_DP
    real(DP)     :: del2  = 0.0_DP
    real(DP)     :: del3  = 0.0_DP
    real(DP)     :: didt  = 0.0_DP
    real(DP)     :: dmdt  = 0.0_DP
    real(DP)     :: dnodt = 0.0_DP
    real(DP)     :: domdt = 0.0_DP
    real(DP)     :: e3    = 0.0_DP
    real(DP)     :: ee2   = 0.0_DP
    real(DP)     :: peo   = 0.0_DP
    real(DP)     :: pgho  = 0.0_DP
    real(DP)     :: pho   = 0.0_DP
    real(DP)     :: pinco = 0.0_DP
    real(DP)     :: plo   = 0.0_DP
    real(DP)     :: se2   = 0.0_DP
    real(DP)     :: se3   = 0.0_DP
    real(DP)     :: sgh2  = 0.0_DP
    real(DP)     :: sgh3  = 0.0_DP
    real(DP)     :: sgh4  = 0.0_DP
    real(DP)     :: sh2   = 0.0_DP
    real(DP)     :: sh3   = 0.0_DP
    real(DP)     :: si2   = 0.0_DP
    real(DP)     :: si3   = 0.0_DP
    real(DP)     :: sl2   = 0.0_DP
    real(DP)     :: sl3   = 0.0_DP
    real(DP)     :: sl4   = 0.0_DP
    real(DP)     :: gsto  = 0.0_DP
    real(DP)     :: xfact = 0.0_DP
    real(DP)     :: xgh2  = 0.0_DP
    real(DP)     :: xgh3  = 0.0_DP
    real(DP)     :: xgh4  = 0.0_DP
    real(DP)     :: xh2   = 0.0_DP
    real(DP)     :: xh3   = 0.0_DP
    real(DP)     :: xi2   = 0.0_DP
    real(DP)     :: xi3   = 0.0_DP
    real(DP)     :: xl2   = 0.0_DP
    real(DP)     :: xl3   = 0.0_DP
    real(DP)     :: xl4   = 0.0_DP
    real(DP)     :: xlamo = 0.0_DP
    real(DP)     :: zmol  = 0.0_DP
    real(DP)     :: zmos  = 0.0_DP
    real(DP)     :: atime = 0.0_DP
    real(DP)     :: xli   = 0.0_DP
    real(DP)     :: xni   = 0.0_DP
    ! error
    integer(SP)  :: error = 0
  end type t_sat
end module model


!*******************************************************************************
! Read the TLE earth satellite file format.
!
!   date          name            version
!   2018.11.22    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module io
  use const,       only : SP, DP, PI, DEG2RAD, XPDOTP, t_cst
  use model,       only : t_sat
  use ext,         only : days2mdhms, jday
  use propagation, only : sgp4init
  implicit none
  private
  public :: twoline2rv

contains
  ! ============================
  ! Public subroutines/functions
  ! ============================

  ! Subroutine twoline2rv
  ! * this function converts the two line element set character string data to
  !   variables and initializes the sgp4 variables. several intermediate varaibles
  !   and quantities are determined. note that the result is a structure so multiple
  !   satellites can be processed simultaneously without having to reinitialize. the
  !   verification mode is an important option that permits quick checks of any
  !   changes to the underlying technical theory. this option works using a
  !   modified tle file in which the start, stop, and delta time values are
  !   included at the end of the second line of data. this only works with the
  !   verification mode. the catalog mode simply propagates from -1440 to 1440 min
  !   from epoch and is useful when performing entire catalog runs.
  !
  ! :param(in)  character(*)   tle(2): first and second line of the TLE
  ! :param(in)  type(t_cst) gravconst: which set of constants to use  72, 84
  !                                    `wgs72`    - Standard WGS 72 model
  !                                    `wgs84`    - More recent WGS 84 model
  !                                    `wgs72old` - Legacy support for old SGP4 behavior
  ! :param(out) type(t_sat)    satrec: a Satellite imported from two lines of TLE
  ! :param(in)  logical, optional  afspc_mode: 
  !               Normally, computations are made using various recent improvements
  !               to the algorithm.  If you want to turn some of these off and go
  !               back into "afspc" mode, then set `afspc_mode` to `.true.`.
  subroutine twoline2rv(tle, gravconst, satrec, afspc_mode)
    implicit none
    character(*), intent(in)  :: tle(2)
    type(t_cst),  intent(in)  :: gravconst
    type(t_sat),  intent(out) :: satrec
    logical,      intent(in), optional  :: afspc_mode
    character(1) :: opsmode, nddot_1, bstar_1
    character(5) :: nddot_2, bstar_2
    character(7) :: nddot_a, bstar_a
    integer(SP)  :: satnum, two_digit_year, nexp, ibexp, ecco_i
    integer(SP)  :: year, month, day, hour, minute, sec_whole
    real(DP)     :: second, sec_fraction
    real(DP)     :: tumin, epochdays, ndot, nddot, bstar
    real(DP)     :: inclo, nodeo, ecco, argpo, mo, no

    opsmode = "i"
    if (present(afspc_mode) .and. (afspc_mode)) opsmode = "a"

    tumin = gravconst%tumin

    ! TLE reading
    ! 1 NNNNNC NNNNNAAA NNNNN.NNNNNNNN +.NNNNNNNN +NNNNN-N +NNNNN-N N NNNNN
    ! 2 NNNNN NNN.NNNN NNN.NNNN NNNNNNN NNN.NNNN NNN.NNNN NN.NNNNNNNNNNNNNN
    ! (First line)
    read (tle(1), '(XXI5XX(8X)XI2F12.8XF10.8XA1A5I2XA1A5I2)') &
      & satnum, two_digit_year, epochdays, ndot, &
      & nddot_1, nddot_2, nexp, bstar_1, bstar_2, ibexp
    nddot_a = nddot_1 // "." // nddot_2
    bstar_a = bstar_1 // "." // bstar_2
    read (nddot_a, '(F7.5)') nddot
    read (bstar_a, '(F7.5)') bstar
    satrec%satnum    = satnum
    satrec%epochdays = epochdays
    satrec%ndot      = ndot
    satrec%nddot     = nddot
    satrec%bstar     = bstar
    ! (Second line)
    read (tle(2), '((8X)F8.4XF8.4XI7XF8.4XF8.4XF11.8)') &
      & inclo, nodeo, ecco_i, argpo, mo, no
    ecco = ecco_i * 1.0e-7_DP
    satrec%inclo = inclo
    satrec%nodeo = nodeo
    satrec%ecco  = ecco
    satrec%argpo = argpo
    satrec%mo    = mo
    satrec%no    = no

    ! ---- find no, ndot, nddot ----
    satrec%no    = no / XPDOTP  ! rad/min
    satrec%nddot = nddot * (10.0_DP**nexp)
    satrec%bstar = bstar * (10.0_DP**ibexp)

    ! ---- convert to sgp4 units ----
    satrec%ndot  = ndot  / (XPDOTP * 1440.0_DP)  ! ? * minperday
    satrec%nddot = nddot / (XPDOTP * 1440.0_DP * 1440.0_DP)

    ! ---- find standard orbital elements ----
    satrec%inclo = inclo  * DEG2RAD
    satrec%nodeo = nodeo  * DEG2RAD
    satrec%argpo = argpo  * DEG2RAD
    satrec%mo    = mo     * DEG2RAD

    ! ----------------------------------------------------------------
    ! find sgp4epoch time of element set
    ! remember that sgp4 uses units of days from 0 jan 1950 (sgp4epoch)
    ! and minutes from the epoch (time)
    ! ----------------------------------------------------------------

    ! ---------------- temp fix for years from 1957-2056 -------------------
    ! --------- correct fix will occur when year is 4-digit in tle ---------
    if (two_digit_year < 57) then
      year = two_digit_year + 2000
    else
      year = two_digit_year + 1900
    end if

    call days2mdhms(year, satrec%epochdays, &
      & month, day, hour, minute, second)
    sec_whole         = int(second)
    sec_fraction      = second - real(sec_whole, DP)
    satrec%epochyr    = year
    satrec%jdsatepoch = jday(year, month, day, hour, minute, second)

    ! ---------------- initialize the orbit at sgp4epoch -------------------
    call sgp4init( &
      & gravconst, opsmode, satrec%satnum, &
      & satrec%jdsatepoch - 2433281.5_DP, &
      & satrec%bstar, satrec%ecco, satrec%argpo, satrec%inclo, &
      & satrec%mo, satrec%no, satrec%nodeo, satrec &
    & )
  end subroutine twoline2rv
end module io


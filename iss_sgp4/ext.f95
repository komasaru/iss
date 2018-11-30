!*******************************************************************************
! Utility routines
!
!   date          name            version
!   2018.11.22    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module ext
  use const, only : SP, DP, DAYS, UID_E, F_EOP, UID_D, F_DAT, FMT_DT_2
  implicit none
  private
  public :: jday, days2mdhms, add_day, gc2jd, get_eop, get_dat, &
          & v_cross, date_fmt
  type, public :: t_time
    integer(SP) :: year    = 0
    integer(SP) :: month   = 0
    integer(SP) :: day     = 0
    integer(SP) :: hour    = 0
    integer(SP) :: minute  = 0
    integer(SP) :: second  = 0
    integer(SP) :: msecond = 0
  end type t_time
  type, public :: t_eop
    character(10) :: date     = "0000-00-00"
    real(DP)      :: mjd      = 0.0_DP
    character(1)  :: flag_pm  = ""
    real(DP)      :: pm_x     = 0.0_DP
    real(DP)      :: pm_x_e   = 0.0_DP
    real(DP)      :: pm_y     = 0.0_DP
    real(DP)      :: pm_y_e   = 0.0_DP
    character(1)  :: flag_dut = ""
    real(DP)      :: dut1     = 0.0_DP
    real(DP)      :: dut1_e   = 0.0_DP
    real(DP)      :: lod      = 0.0_DP
    real(DP)      :: lod_e    = 0.0_DP
    character(1)  :: flag_nut = ""
    real(DP)      :: nut_x    = 0.0_DP
    real(DP)      :: nut_x_e  = 0.0_DP
    real(DP)      :: nut_y    = 0.0_DP
    real(DP)      :: nut_y_e  = 0.0_DP
  end type t_eop

contains
  ! =============================
  ! Private subroutines/functions
  ! =============================

  ! うるう年判定
  !
  ! :param(in) integer(4)   ye: 年
  ! :return    logical is_leap: うるう年(T), うるう年でない(F)
  logical function is_leap(ye)
    implicit none
    integer(SP), intent(in) :: ye

    if (mod(ye, 400) == 0) then
      is_leap = .true.
    else
      if (mod(ye, 4) == 0 .and. mod(ye, 100) /= 0) then
        is_leap = .true.
      else
        is_leap = .false.
      end if
    end if
  end function is_leap

  ! JD(Julian Day) -> GC(Gregoria Calendar)
  !
  ! :param(in)  real(8)      jd: Julian Day
  ! :param(out) type(t_time) gc: Gregoria Calendar
  subroutine jd2gc(jd, gc)
    implicit none
    real(DP),     intent(in)  :: jd
    type(t_time), intent(out) :: gc
    integer(SP) :: i, n, a, b, tm(0:6)
    real(DP)    :: tm_f, tm_w

    ! n, a, b 計算
    n = int(jd - 1721119.5_DP)  ! = -2400000.5 + 678881.0
    a = 4 * n + 3 + 4 * floor((3.0_DP / 4.0_DP) &
      & * (floor(4 * (n + 1) / 146097.0_DP) + 1))
    b = 5 * floor(mod(a, 1461) / 4.0_DP) + 2

    ! 年・月・日 計算
    tm(0) = floor(a / 1461.0_DP)
    tm(1) = floor(b / 153.0_DP)
    tm(2) = floor(mod(b, 153) / 5.0_DP)
    tm_w  = floor((tm(1) + 2) / 12.0_DP)
    tm(0) = tm(0) + tm_w
    tm(1) = tm(1) + 2 - tm_w * 12 + 1
    tm(2) = tm(2) + 1

    ! 時・分・秒・ミリ秒 計算
    tm_f = 86400.0_DP * (jd - .5_DP - int(jd - .5_DP))
    tm(3) = int(tm_f / 3600.0_DP)
    tm(4) = int((tm_f - 3600 * tm(3)) / 60.0_DP)
    tm_w = tm_f - 3600 * tm(3) - 60 * tm(4)
    tm(5) = int(tm_w)
    tm(6) = nint((tm_w - tm(5)) * 1.0e3_DP)
    ! ミリ秒四捨五入で 1000 になった場合
    if (tm(6) > 999) then
      tm(5) = tm(5) + 1
      tm(6) = tm(6) - 1000
      if (tm(5) > 59) then
        tm(4) = tm(4) + 1
        tm(5) = tm(5) - 60
        if (tm(4) > 59) then
          tm(3) = tm(3) + 1
          tm(4) = tm(4) - 60
        end if
      end if
    end if
    gc = t_time(tm(0), tm(1), tm(2), tm(3), tm(4), tm(5), tm(6))
  end subroutine jd2gc

  ! ============================
  ! Public subroutines/functions
  ! ============================

  ! Subroutine jday
  ! *  this procedure finds the julian date given the year, month, day, and time.
  !    the julian date is defined by each elapsed day since noon, jan 1, 4713 bc.
  ! *  algorithm     : calculate the answer in one step for efficiency
  !
  ! :param(in)  integer(4) year: year                (1900 .. 2100)
  ! :param(in)  integer(4)  mon: month               (1 .. 12)
  ! :param(in)  integer(4)  day: day                 (1 .. 28,29,30,31)
  ! :param(in)  integer(4)   hr: universal time hour (0 .. 23)
  ! :param(in)  integer(4)  min: universal time min  (0 .. 59)
  ! :param(in)  real(8)     sec: universal time sec  (0.0 .. 59.999)
  ! :return     real(8)    jday: julian date         (days from 4713 bc)
  real(DP) function jday(year, mon, day, hr, minute, sec)
    implicit none
    integer(SP), intent(in) :: year, mon, day, hr, minute
    real(DP),    intent(in) :: sec

    jday = (367.0_DP * year &
       & - int(7.0_DP * (year + int((mon + 9.0_DP) / 12.0_DP)) * 0.25_DP) &
       & + int(275.0_DP * mon / 9.0_DP) &
       & + day + 1721013.5_DP &
       & + ((sec / 60.0_DP + minute) / 60.0_DP + hr) / 24.0_DP)
  end function jday

  ! Subroutine days2mdhms
  ! * this procedure converts the day of the year, days, to the equivalent month
  !   day, hour, minute and second.
  ! * algorithm: set up array for the number of days per month
  !              find leap year - use 1900 because 2000 is a leap year
  !              loop through a temp value while the value is < the days
  !              perform int conversions to the correct day and month
  !              convert remainder into h m s using type conversions

  ! :pamra(in)  integer(4)   year: year    (1900 .. 2100)
  ! :param(in)  real(8)      days: Fractional days into the year of 
  !                                the epoch moment
  ! :param(out) integer(4)  month: month   (1 .. 12)
  ! :param(out) integer(4)    day: day     (1 .. 28,29,30,31)
  ! :param(out) integer(4)   hour: hour    (0 .. 23)
  ! :param(out) integer(4) minute: minute  (0 .. 59)
  ! :param(out) real(8)    second: second  (0.0 .. 59.999)
  subroutine days2mdhms(year, days, month, day, hour, minute, second)
    implicit none
    integer(SP), intent(in)  :: year
    real(DP),    intent(in)  :: days
    integer(SP), intent(out) :: month, day, hour, minute
    real(DP),    intent(out) :: second
    integer(SP) :: lmonth(12), dayofyr, i, inttemp
    real(DP)    :: temp

    dayofyr = int(days)

    ! ----------------- find month and day of month ----------------
    if (mod(year, 4) == 0) then
      lmonth = (/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    else
      lmonth = (/31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
    end if

    i = 1
    inttemp = 0
    do while (dayofyr > inttemp + lmonth(i) .and. i < 13)
      inttemp = inttemp + lmonth(i)
      i = i + 1
    end do

    month = i
    day   = dayofyr - inttemp

    ! ----------------- find hours minutes and seconds -------------
    temp   = (days - dayofyr) * 24.0_DP
    hour   = int(temp)
    temp   = (temp - hour) * 60.0_DP
    minute = int(temp)
    second = (temp - minute) * 60.0_DP
  end subroutine days2mdhms

  ! 日の加減算
  !
  ! :param(in)  type(t_time)     gc: Gregorian Calendar
  ! :param(in)  real(8)           d: days
  ! :param(out) type(t_time) gc_dst: Gregorian Calendar (計算後)
  subroutine add_day(gc, d, gc_dst)
    implicit none
    type(t_time), intent(in)  :: gc
    real(DP),     intent(in)  :: d
    type(t_time), intent(out) :: gc_dst
    real(DP) :: jd

    call gc2jd(gc, jd)
    call jd2gc(jd + d, gc_dst)
  end subroutine add_day

  ! GC(Gregoria Calendar) -> JD(Julian Day)
  !
  ! :param(in)  type(t_time) utc
  ! :param(out) real(8)       jd
  subroutine gc2jd(utc, jd)
    implicit none
    type(t_time), intent(in)  :: utc
    real(DP),     intent(out) :: jd
    integer  :: ye, mo, da, ho, mi, se, ms
    real(DP) :: d, t

    ye = utc%year
    mo = utc%month
    da = utc%day
    ho = utc%hour
    mi = utc%minute
    se = utc%second
    ms = utc%msecond

    if (mo < 3) then
      ye= ye - 1
      mo= mo + 12
    end if
    d = int(365.25_DP * ye)      &
    & + int(ye / 400.0_DP)       &
    & - int(ye / 100.0_DP)       &
    & + int(30.59_DP * (mo - 2)) &
    & + da + 1721088.5_DP
    t = (ms / (3600.0_DP * 1.0e3_DP) &
    & + se / 3600.0_DP               &
    & + mi / 60.0_DP                 &
    & + ho) / 24.0_DP
    jd = d + t
  end subroutine gc2jd

  ! EOP （地球回転パラメータ）取得
  ! * 予め作成しておいた CSV データから取得する
  !
  ! :param(in)  type(t_time) utc: UTC
  ! :param(out) type(t_eop)  eop: EOP
  subroutine get_eop(utc, eop)
    implicit none
    type(t_time), intent(in)  :: utc
    type(t_eop),  intent(out) :: eop
    character(10)  :: utc_t
    character(128) :: buf
    integer(SP)    :: ios
    character(10) :: date
    character(1)  :: flag_pm, flag_dut, flag_nut
    real(DP)      :: mjd, pm_x, pm_x_e, pm_y, pm_y_e
    real(DP)      :: dut1, dut1_e, lod, lod_e
    real(DP)      :: nut_x, nut_x_e, nut_y, nut_y_e

    ! 対象の UTC 年月日
    write (utc_t, '(I4, "-", I0.2, "-", I0.2)') utc%year, utc%month, utc%day

    open (unit   = UID_E,       &
        & iostat = ios,         &
        & file   = F_EOP,       &
        & action = "read",      &
        & form   = "formatted", &
        & status = "old")
    if (ios /= 0) then
      print '("[ERROR:", I0 ,"] Failed to open file: ", A)', ios, F_EOP
      stop
    end if

    do
      eop = t_eop( &
        & "0000-00-00", 0.0_DP, &
        & "", 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        & "", 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
        & "", 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP  &
      & )
      read (UID_E, '(A)', iostat = ios) buf
      if (ios < 0) then
        exit
      else if (ios > 0) then
        print '("[ERROR:", I0 ,"] Failed to read file: ", A)', ios, F_EOP
      end if
      read (buf, *, iostat = ios) eop
      if (ios /= 0) exit
      if (utc_t == eop%date) exit
    end do

    close(UID_E)
  end subroutine get_eop

  ! DAT (= TAI - UTC)（うるう秒の総和）取得
  ! * 予め取得しておいたテキストファイルから取得する
  !
  ! :param(in)  type(t_time) utc: UTC
  ! :param(out) integer(4)   dat: DAT
  subroutine get_dat(utc, dat)
    implicit none
    type(t_time), intent(in)  :: utc
    integer(SP),  intent(out) :: dat
    character(10)  :: utc_t, utc_w
    character(128) :: buf
    integer(SP)    :: ios
    integer(SP)    :: d, m, y, dat_w, dat_p
    real(DP)       :: mjd

    ! 対象の UTC 年月日
    write (utc_t, '(I4, "-", I0.2, "-", I0.2)') utc%year, utc%month, utc%day

    open (unit   = UID_D,       &
        & iostat = ios,         &
        & file   = F_DAT,       &
        & action = "read",      &
        & form   = "formatted", &
        & status = "old")
    if (ios /= 0) then
      print '("[ERROR:", I0 ,"] Failed to open file: ", A)', ios, F_DAT
      stop
    end if

    do
      read (UID_D, '(A)', iostat = ios) buf
      if (ios < 0) then
        exit
      else if (ios > 0) then
        print '("[ERROR:", I0 ,"] Failed to read file: ", A)', ios, F_DAT
      end if
      if (buf(1:1) == "#") cycle
      read (buf, *) mjd, d, m, y, dat_w
      write (utc_w, '(I4, "-", I0.2, "-", I0.2)') y, m, d
      if (utc_t < utc_w) exit
      dat_p = dat_w
    end do
    dat = dat_p

    close(UID_D)
  end subroutine get_dat

  ! ベクトル外積の計算
  !
  ! :param(in) real(8) a(3): ベクトルA
  ! :param(in) real(8) b(3): ベクトルB
  ! :return    real(8) c(3): ベクトル
  function v_cross(a, b) result(c)
    implicit none
    real(DP), intent(in) :: a(3), b(3)
    real(DP) :: c(3)

    c = (/ &
      & a(2) * b(3) - a(3) * b(2), &
      & a(3) * b(1) - a(1) * b(3), &
      & a(1) * b(2) - a(2) * b(1)  &
    & /)
  end function v_cross

  ! 日付文字列の整形
  ! * type(t_time)型 -> YYYY-MM-DD HH:MM:SS.MMM
  !
  ! :param(in) type(t_time)  d
  ! :return    character(23) f
  function date_fmt(d) result(f)
    type(t_time), intent(in) :: d
    character(23) :: f

    write (f, FMT_DT_2) &
      & d%year, d%month, d%day, d%hour, d%minute, d%second, d%msecond
  end function date_fmt
end module ext


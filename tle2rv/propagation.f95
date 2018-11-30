!*******************************************************************************
! The sgp4 procedures for analytical propagation of a satellite.
!
! I have made the rather unorthodox decision to leave as much of this C++
! code alone as possible: if a line of code would run without change in
! Python, then I refused to re-indent it or remove its terminal semicolon,
! so that in the future it will be easier to keep updating this file as
! the original author's C++ continues to improve.  Thus, 5-space
! indentation (!) prevails in this file.
!
! I have even kept all of the C++ block comments (by turning them into
! Python string constants) to make this easier to navigate and maintain,
! as well as to make it more informative for people who encounter this
! code for the first time here in its Python form.
!
!   date          name            version
!   2018.11.23    mk-mode.com     1.00 新規作成
!
! Copyright(C) 2018 mk-mode.com All Rights Reserved.
!*******************************************************************************
!
module propagation
  use const, only : SP, DP, PI, PI2, DEG2RAD, MIN_D, t_cst
  use model, only : t_sat
  use ext,   only : t_time, jday
  implicit none
  private
  public :: propagate, sgp4init, sgp4

contains
  ! ============================
  ! Private subroutines/functions
  ! ============================

  ! Subroutine initl
  ! * this procedure initializes the spg4 propagator. all the initialization is
  !   consolidated here instead of having multiple loops inside other routines.
  !
  ! :param(in)    integer(4)       satn: satellite number
  ! :param(in)    type(t_cst) gravconst: constants
  ! :param(in)    real(8)          ecco: eccentricity  (0.0 - 1.0)
  ! :param(in)    real(8)         epoch: epoch time in days from jan 0, 1950. 0 hr
  ! :param(in)    real(8)         inclo: inclination of satellite
  ! :param(in)    character(1)  opsmode: mode of operation afspc or improved 'a', 'i'
  ! :param(inout) real(8)            no: mean motion of satellite
  ! :param(inout) character(1)   method: flag for deep space  ('d', 'n')
  ! :param(out)   real(8)          ainv: 1.0 / a
  ! :param(out)   real(8)            ao: semi major axis
  ! :param(out)   real(8)         con41:
  ! :param(out)   real(8)         con42: 1.0 - 5.0 cos(i)
  ! :param(out)   real(8)         cosio: cosine of inclination
  ! :param(out)   real(8)        cosio2: cosio squared
  ! :param(out)   real(8)         eccsq: eccentricity squared
  ! :param(out)   real(8)        omeosq: 1.0 - ecco * ecco
  ! :param(out)   real(8)          posq: semi-parameter squared
  ! :param(out)   real(8)            rp: radius of perigee
  ! :param(out)   real(8)        rteosq: square root of (1.0 - ecco*ecco)
  ! :param(out)   real(8)         sinio: sine of inclination
  ! :param(out)   real(8)          gsto: gst at time of observation  (rad)
  ! :param(out)   real(8)            no: mean motion of satellite
  subroutine initl( &
    & satn, gravconst, ecco, epoch, inclo, opsmode,             &
    & no, method, ainv, ao, con41, con42, cosio, cosio2, eccsq, &
    & omeosq, posq, rp, rteosq, sinio, gsto                     &
  & )
    implicit none
    integer(SP),  intent(in)    :: satn
    type(t_cst),  intent(in)    :: gravconst
    real(DP),     intent(in)    :: ecco, epoch, inclo
    character(1), intent(in)    :: opsmode
    real(DP),     intent(inout) :: no
    character(1), intent(inout) :: method
    real(DP),     intent(out)   :: ainv, ao, con41, con42, cosio, cosio2
    real(DP),     intent(out)   :: eccsq, omeosq, posq, rp, rteosq
    real(DP),     intent(out)   :: sinio, gsto
    integer(SP) :: ds70
    real(DP)    :: x2o3, ak, d1, del_, adel, po, ts70
    real(DP)    :: tfrac, c1, thgr70, fk5r, c1p2p

    ! sgp4fix use old way of finding gst

    ! ----------------------- earth constants ----------------------
    x2o3 = 2.0_DP / 3.0_DP

    ! ------------- calculate auxillary epoch quantities ----------
    eccsq  = ecco * ecco
    omeosq = 1.0_DP - eccsq
    rteosq = sqrt(omeosq)
    cosio  = cos(inclo)
    cosio2 = cosio * cosio

    ! ------------------ un-kozai the mean motion -----------------
    ak   = (gravconst%xke / no)**x2o3
    d1   = 0.75_DP * gravconst%j2 * (3.0_DP * cosio2 - 1.0_DP) &
       & / (rteosq * omeosq)
    del_ = d1 / (ak * ak)
    adel = ak * (1.0_DP - del_ * del_ - del_ &
       & * (1.0_DP / 3.0_DP + 134.0_DP * del_ * del_ / 81.0_DP))
    del_ = d1 / (adel * adel)
    no   = no / (1.0_DP + del_)

    ao     = (gravconst%xke / no) ** x2o3
    sinio  = sin(inclo)
    po     = ao * omeosq
    con42  = 1.0_DP - 5.0_DP * cosio2
    con41  = - con42 - cosio2 - cosio2
    ainv   = 1.0_DP / ao
    posq   = po * po
    rp     = ao * (1.0_DP - ecco)
    method = "n"

    ! sgp4fix modern approach to finding sidereal time
    if (opsmode == "a") then
      ! sgp4fix use old way of finding gst
      ! count integer number of days from 0 jan 1970
      ts70  = epoch - 7305.0_DP
      ds70  = int(ts70 + 1.0e-8_DP)
      tfrac = ts70 - real(ds70, DP)
      ! find greenwich location at epoch
      c1     = 1.72027916940703639e-2_DP
      thgr70 = 1.7321343856509374_DP
      fk5r   = 5.07551419432269442e-15_DP
      c1p2p  = c1 + PI2
      gsto   = mod( &
        & thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r, &
        & PI2 &
      & )
      if (gsto < 0.0_DP) gsto = gsto + PI2
    else
      gsto = gstime(epoch + 2433281.5_DP)
    end if
  end subroutine initl

  ! Subroutine dscom
  ! * this procedure provides deep space common items used by both the secular
  !   and periodics subroutines.  input is provided as shown. this routine
  !   used to be called dpper, but the functions inside weren't well organized.
  !
  ! :param(in)    real(8) epoch:
  ! :param(in)    real(8)  eccp: eccentricity
  ! :param(in)    real(8) argpp: argument of perigee
  ! :param(in)    real(8)    tc:
  ! :param(in)    real(8) inclp: inclination
  ! :param(in)    real(8) nodep: right ascension of ascending node
  ! :param(in)    real(8)    np: mean motion
  ! :param(inout) real(8) e3, ee2, peo, pgho, pho, pinco, plo,
  !                       se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3,
  !                       sl2, sl3, sl4, xgh2, xgh3, xgh4, xh2, xh3,
  !                       xi2, xi3, xl2, xl3, xl4, zmol, zmos
  ! :param(out)   real(8) snodm, cnodm, sinim, cosim, sinomm, cosomm, day,
  !                       em, emsq, gam, rtemsq, s1, s2, s3, s4, s5,
  !                       s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1,
  !                       sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23,
  !                       sz31, sz32, sz33, nm, z1, z2, z3,
  !                       z11, z12, z13, z21, z22, z23, z31, z32, z33
  subroutine dscom( &
    & epoch, eccp, argpp, tc, inclp, nodep, np,        &
    & e3, ee2, peo, pgho, pho, pinco, plo,             &
    & se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3,  &
    & sl2, sl3, sl4, xgh2, xgh3, xgh4, xh2, xh3,       &
    & xi2, xi3, xl2, xl3, xl4, zmol, zmos,             &
    & snodm, cnodm, sinim, cosim, sinomm, cosomm, day, &
    & em, emsq, gam, rtemsq, s1, s2, s3, s4, s5,       &
    & s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1,  &
    & sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23,    &
    & sz31, sz32, sz33, nm, z1, z2, z3,                &
    & z11, z12, z13, z21, z22, z23, z31, z32, z33      &
  & )
    implicit none
    real(DP), intent(in)    :: epoch, eccp, argpp, tc, inclp, nodep, np
    real(DP), intent(inout) :: e3, ee2, peo, pgho, pho, pinco, plo
    real(DP), intent(inout) :: se2, se3, sgh2, sgh3, sgh4, sh2, sh3, si2, si3
    real(DP), intent(inout) :: sl2, sl3, sl4, xgh2, xgh3, xgh4, xh2, xh3
    real(DP), intent(inout) :: xi2, xi3, xl2, xl3, xl4, zmol, zmos
    real(DP), intent(out)   :: snodm, cnodm, sinim, cosim, sinomm, cosomm, day
    real(DP), intent(out)   :: em, emsq, gam, rtemsq
    real(DP), intent(out)   :: s1, s2, s3, s4, s5, s6, s7
    real(DP), intent(out)   :: ss1, ss2, ss3, ss4, ss5, ss6, ss7
    real(DP), intent(out)   :: sz1, sz2, sz3, sz11, sz12, sz13
    real(DP), intent(out)   :: sz21, sz22, sz23, sz31, sz32, sz33
    real(DP), intent(out)   :: nm, z1, z2, z3
    real(DP), intent(out)   :: z11, z12, z13, z21, z22, z23, z31, z32, z33
    ! constants
    real(DP), parameter :: zes    =  0.01675_DP
    real(DP), parameter :: zel    =  0.05490_DP
    real(DP), parameter :: c1ss   =  2.9864797e-6_DP
    real(DP), parameter :: c1l    =  4.7968065e-7_DP
    real(DP), parameter :: zsinis =  0.39785416_DP
    real(DP), parameter :: zcosis =  0.91744867_DP
    real(DP), parameter :: zcosgs =  0.1945905_DP
    real(DP), parameter :: zsings = -0.98088458_DP
    ! variables
    integer(SP) :: lsflg
    real(DP)    :: betasq, ctem, stem, xnodce, zx, zy
    real(DP)    :: zcosgl, zcoshl, zcosil, zsingl, zsinhl, zsinil
    real(DP)    :: cc, xnoi, zcosg, zcosh, zcosi, zsing, zsinh, zsini
    real(DP)    :: a(10), x(8)

    ! --------------------- local variables ------------------------
    nm     = np
    em     = eccp
    snodm  = sin(nodep)
    cnodm  = cos(nodep)
    sinomm = sin(argpp)
    cosomm = cos(argpp)
    sinim  = sin(inclp)
    cosim  = cos(inclp)
    emsq   = em * em
    betasq = 1.0_DP - emsq
    rtemsq = sqrt(betasq)

    ! ----------------- initialize lunar solar terms ---------------
    peo    = 0.0_DP
    pinco  = 0.0_DP
    plo    = 0.0_DP
    pgho   = 0.0_DP
    pho    = 0.0_DP
    day    = epoch + 18261.5_DP + tc / 1440.0_DP
    xnodce = mod(4.5236020 - 9.2422029e-4 * day, PI2)
    stem   = sin(xnodce)
    ctem   = cos(xnodce)
    zcosil = 0.91375164_DP - 0.03568096_DP * ctem
    zsinil = sqrt(1.0_DP - zcosil * zcosil)
    zsinhl = 0.089683511_DP * stem / zsinil
    zcoshl = sqrt(1.0_DP - zsinhl * zsinhl)
    gam    = 5.8351514_DP + 0.0019443680_DP * day
    zx     = 0.39785416_DP * stem / zsinil
    zy     = zcoshl * ctem + 0.91744867_DP * zsinhl * stem
    zx     = atan2(zx, zy)
    zx     = gam + zx - xnodce
    zcosgl = cos(zx)
    zsingl = sin(zx)

    ! ------------------------- do solar terms ---------------------
    zcosg = zcosgs
    zsing = zsings
    zcosi = zcosis
    zsini = zsinis
    zcosh = cnodm
    zsinh = snodm
    cc    = c1ss
    xnoi  = 1.0_DP / nm

    do lsflg = 1, 2
      a( 1) =   zcosg * zcosh + zsing * zcosi * zsinh
      a( 3) =  -zsing * zcosh + zcosg * zcosi * zsinh
      a( 7) =  -zcosg * zsinh + zsing * zcosi * zcosh
      a( 8) =   zsing * zsini
      a( 9) =   zsing * zsinh + zcosg * zcosi * zcosh
      a(10) =   zcosg * zsini
      a( 2) =   cosim * a(7) + sinim * a( 8)
      a( 4) =   cosim * a(9) + sinim * a(10)
      a( 5) =  -sinim * a(7) + cosim * a( 8)
      a( 6) =  -sinim * a(9) + cosim * a(10)

      x( 1) =  a(1) * cosomm + a(2) * sinomm
      x( 2) =  a(3) * cosomm + a(4) * sinomm
      x( 3) = -a(1) * sinomm + a(2) * cosomm
      x( 4) = -a(3) * sinomm + a(4) * cosomm
      x( 5) =  a(5) * sinomm
      x( 6) =  a(6) * sinomm
      x( 7) =  a(5) * cosomm
      x( 8) =  a(6) * cosomm

      z31 = 12.0_DP * x(1) * x(1) - 3.0_DP * x(3) * x(3)
      z32 = 24.0_DP * x(1) * x(2) - 6.0_DP * x(3) * x(4)
      z33 = 12.0_DP * x(2) * x(2) - 3.0_DP * x(4) * x(4)
      z1  =  3.0_DP *  (a(1) * a(1) + a(2) * a(2)) + z31 * emsq
      z2  =  6.0_DP *  (a(1) * a(3) + a(2) * a(4)) + z32 * emsq
      z3  =  3.0_DP *  (a(3) * a(3) + a(4) * a(4)) + z33 * emsq
      z11 = -6.0_DP * a(1) * a(5) &
        & + emsq *  (-24.0 * x(1) * x(7) -6.0 * x(3) * x(5))
      z12 = -6.0_DP *  (a(1) * a(6) + a(3) * a(5)) &
        & + emsq * (-24.0_DP * (x(2) * x(7) + x(1) * x(8)) &
        & - 6.0_DP * (x(3) * x(6) + x(4) * x(5)))
      z13 = -6.0_DP * a(3) * a(6) &
        & + emsq * (-24.0_DP * x(2) * x(8) - 6.0_DP * x(4) * x(6))
      z21 =  6.0_DP * a(2) * a(5) &
        & + emsq * ( 24.0_DP * x(1) * x(5) - 6.0_DP * x(3) * x(7))
      z22 =  6.0_DP *  (a(4) * a(5) + a(2) * a(6)) &
        & + emsq * (24.0_DP * (x(2) * x(5) + x(1) * x(6)) &
        & - 6.0_DP * (x(4) * x(7) + x(3) * x(8)))
      z23 =  6.0_DP * a(4) * a(6) &
        & + emsq * (24.0_DP * x(2) * x(6) - 6.0_DP * x(4) * x(8))
      z1  = z1 + z1 + betasq * z31
      z2  = z2 + z2 + betasq * z32
      z3  = z3 + z3 + betasq * z33
      s3  = cc * xnoi
      s2  = -0.5_DP * s3 / rtemsq
      s4  = s3 * rtemsq
      s1  = -15.0_DP * em * s4
      s5  = x(1) * x(3) + x(2) * x(4)
      s6  = x(2) * x(3) + x(1) * x(4)
      s7  = x(2) * x(4) - x(1) * x(3)

      ! ----------------------- do lunar terms -------------------
      if (lsflg == 1) then
        ss1   = s1
        ss2   = s2
        ss3   = s3
        ss4   = s4
        ss5   = s5
        ss6   = s6
        ss7   = s7
        sz1   = z1
        sz2   = z2
        sz3   = z3
        sz11  = z11
        sz12  = z12
        sz13  = z13
        sz21  = z21
        sz22  = z22
        sz23  = z23
        sz31  = z31
        sz32  = z32
        sz33  = z33
        zcosg = zcosgl
        zsing = zsingl
        zcosi = zcosil
        zsini = zsinil
        zcosh = zcoshl * cnodm + zsinhl * snodm
        zsinh = snodm * zcoshl - cnodm * zsinhl
        cc    = c1l
      end if
    end do

    zmol = mod(4.7199672_DP + 0.22997150_DP  * day - gam, PI2)
    zmos = mod(6.2565837_DP + 0.017201977_DP * day, PI2)

    ! ------------------------ do solar terms ----------------------
    se2  =   2.0_DP * ss1 * ss6
    se3  =   2.0_DP * ss1 * ss7
    si2  =   2.0_DP * ss2 * sz12
    si3  =   2.0_DP * ss2 * (sz13 - sz11)
    sl2  =  -2.0_DP * ss3 * sz2
    sl3  =  -2.0_DP * ss3 * (sz3 - sz1)
    sl4  =  -2.0_DP * ss3 * (-21.0_DP - 9.0_DP * emsq) * zes
    sgh2 =   2.0_DP * ss4 * sz32
    sgh3 =   2.0_DP * ss4 * (sz33 - sz31)
    sgh4 = -18.0_DP * ss4 * zes
    sh2  =  -2.0_DP * ss2 * sz22
    sh3  =  -2.0_DP * ss2 * (sz23 - sz21)

    ! ------------------------ do lunar terms ----------------------
    ee2  =   2.0_DP * s1 * s6
    e3   =   2.0_DP * s1 * s7
    xi2  =   2.0_DP * s2 * z12
    xi3  =   2.0_DP * s2 * (z13 - z11)
    xl2  =  -2.0_DP * s3 * z2
    xl3  =  -2.0_DP * s3 * (z3 - z1)
    xl4  =  -2.0_DP * s3 * (-21.0 - 9.0 * emsq) * zel
    xgh2 =   2.0_DP * s4 * z32
    xgh3 =   2.0_DP * s4 * (z33 - z31)
    xgh4 = -18.0_DP * s4 * zel
    xh2  =  -2.0_DP * s2 * z22
    xh3  =  -2.0_DP * s2 * (z23 - z21)
  end subroutine dscom

  ! Subroutine dsinit
  ! * this procedure provides deep space contributions to mean motion dot due
  !   to geopotential resonance with half day and one day orbits.
  !
  ! :param(in)    type(t_cst) gravconst
  ! :param(in)    real(8)     cosim, argpo,
  !                           s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4,
  !                           ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33,
  !                           t, tc, gsto, mo, mdot, no, nodeo, nodedot, xpidot,
  !                           z1, z3, z11, z13, z21, z23, z31, z33,
  !                           ecco, eccsq
  ! :param(inout) integer(4)  irez
  ! :param(inout) real(8)     emsq, em, argpm, inclm, mm, nm, nodem, atime,
  !                           d2201, d2211, d3210, d3222, d4410, d4422,
  !                           d5220, d5232, d5421, d5433,
  !                           dedt, didt, dmdt, dnodt, domdt,
  !                           del1,  del2,  del3,  xfact, xlamo, xli, xni
  ! :param(out)   real(8)     dndt
  subroutine dsinit( &
    & gravconst,                                                   &
    & cosim, argpo, s1, s2, s3, s4, s5, sinim, ss1, ss2, ss3, ss4, &
    & ss5, sz1, sz3, sz11, sz13, sz21, sz23, sz31, sz33, t, tc,    &
    & gsto, mo, mdot, no, nodeo, nodedot, xpidot,                  &
    & z1, z3, z11, z13, z21, z23, z31, z33, ecco, eccsq,           &
    & emsq, em, argpm, inclm, mm, nm, nodem, irez,  atime,         &
    & d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232,      &
    & d5421, d5433, dedt, didt, dmdt, dnodt, domdt,                &
    & del1, del2, del3, xfact, xlamo, xli, xni,                    &
    & dndt                                                         &
  & )
    implicit none
    type(t_cst), intent(in)    :: gravconst
    real(DP),    intent(in)    :: cosim, argpo, s1, s2, s3, s4, s5
    real(DP),    intent(in)    :: sinim, ss1, ss2, ss3, ss4
    real(DP),    intent(in)    :: ss5, sz1, sz3, sz11, sz13
    real(DP),    intent(in)    :: sz21, sz23, sz31, sz33, t, tc
    real(DP),    intent(in)    :: gsto, mo, mdot, no, nodeo, nodedot
    real(DP),    intent(in)    :: xpidot, z1, z3, z11, z13, z21, z23
    real(DP),    intent(in)    :: z31, z33, ecco, eccsq
    integer(SP), intent(inout) :: irez
    real(DP),    intent(inout) :: emsq, em, argpm, inclm, mm, nm, nodem
    real(DP),    intent(inout) :: d2201, d2211, d3210, d3222, d4410
    real(DP),    intent(inout) :: d4422, d5220, d5232, d5421, d5433
    real(DP),    intent(inout) :: atime, dedt, didt, dmdt, dnodt, domdt
    real(DP),    intent(inout) :: del1, del2, del3, xfact, xlamo, xli, xni
    real(DP),    intent(out)   :: dndt
    real(DP) :: q22, q31, q33, root22, root32, root44, root52, root54
    real(DP) :: rptim, x2o3, znl, zns, ses, sghs, sgs, shs, sis, sls
    real(DP) :: sghl, shll, theta, aonv, cosisq, emo, emsqo, eoc, g201
    real(DP) :: g200, g211, g300, g310, g322, g410, g422, g520, g521
    real(DP) :: g532, g533, f220, f221, f311, f321, f322, f330, f441
    real(DP) :: f442, f522, f523, f542, f543
    real(DP) :: sini2, ainv2, xno2, temp, temp1

    q22    = 1.7891679e-6_DP
    q31    = 2.1460748e-6_DP
    q33    = 2.2123015e-7_DP
    root22 = 1.7891679e-6_DP
    root44 = 7.3636953e-9_DP
    root54 = 2.1765803e-9_DP
    rptim  = 4.37526908801129966e-3_DP  ! equates to 7.29211514668855e-5 rad/sec
    root32 = 3.7393792e-7_DP
    root52 = 1.1428639e-7_DP
    x2o3   = 2.0_DP / 3.0_DP
    znl    = 1.5835218e-4_DP
    zns    = 1.19459e-5_DP

    ! -------------------- deep space initialization ------------
    irez = 0
    if (0.0034906585_DP < nm .and. nm < 0.0052359877_DP) then
      irez = 1
    else if (8.26e-3_DP <= nm .and. nm <= 9.24e-3_DP &
      & .and. em >= 0.5) then
      irez = 2
    end if

    ! ------------------------ do solar terms -------------------
    ses  =  ss1 * zns * ss5
    sis  =  ss2 * zns * (sz11 + sz13)
    sls  = -zns * ss3 * (sz1 + sz3 - 14.0_DP - 6.0_DP * emsq)
    sghs =  ss4 * zns * (sz31 + sz33 - 6.0_DP)
    shs  = -zns * ss2 * (sz21 + sz23)
    ! sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2_DP .or. inclm > PI - 5.2359877e-2_DP) then
      shs = 0.0_DP
    else if (sinim /= 0.0_DP) then
      shs = shs / sinim
    end if
    sgs  = sghs - cosim * shs

    ! ------------------------- do lunar terms ------------------
    dedt = ses + s1 * znl * s5
    didt = sis + s2 * znl * (z11 + z13)
    dmdt = sls - znl * s3 * (z1 + z3 - 14.0_DP - 6.0_DP * emsq)
    sghl = s4 * znl * (z31 + z33 - 6.0_DP)
    shll = -znl * s2 * (z21 + z23)
    ! sgp4fix for 180 deg incl
    if (inclm < 5.2359877e-2_DP .or. inclm > PI - 5.2359877e-2_DP) then
      shll = 0.0_DP
    end if
    domdt = sgs + sghl
    dnodt = shs
    if (sinim /= 0.0_DP) then
      domdt = domdt - cosim / sinim * shll
      dnodt = dnodt + shll / sinim
    end if

    ! ----------- calculate deep space resonance effects --------
    dndt   = 0.0_DP
    theta  = mod(gsto + tc * rptim, PI2)
    em     = em + dedt * t
    inclm  = inclm + didt * t
    argpm  = argpm + domdt * t
    nodem  = nodem + dnodt * t
    mm     = mm + dmdt * t

    ! -------------- initialize the resonance terms -------------
    if (irez /= 0) aonv = (nm / gravconst%xke)**x2o3

    ! ---------- geopotential resonance for 12 hour orbits ------
    if (irez == 2) then
      cosisq = cosim * cosim
      emo    = em
      em     = ecco
      emsqo  = emsq
      emsq   = eccsq
      eoc    = em * emsq
      g201   = -0.306_DP - (em - 0.64_DP) * 0.440_DP

      if (em <= 0.65_DP) then
        g211 =    3.616_DP  -  13.2470_DP * em +   16.2900_DP * emsq
        g310 =  -19.302_DP  + 117.3900_DP * em -  228.4190_DP * emsq &
           & +  156.5910_DP * eoc
        g322 =  -18.9068_DP + 109.7927_DP * em -  214.6334_DP * emsq &
           & +  146.5816_DP * eoc
        g410 =  -41.122_DP  + 242.6940_DP * em -  471.0940_DP * emsq &
           & +  313.9530_DP * eoc
        g422 = -146.407_DP  + 841.8800_DP * em - 1629.014_DP  * emsq &
           & + 1083.4350_DP * eoc
        g520 = -532.114_DP  + 3017.977_DP * em - 5740.032_DP  * emsq &
           & + 3708.2760_DP * eoc
      else
        g211 =   -72.099_DP +   331.819_DP * em -   508.738_DP * emsq &
           & +   266.724_DP * eoc
        g310 =  -346.844_DP +  1582.851_DP * em -  2415.925_DP * emsq &
           & +  1246.113_DP * eoc
        g322 =  -342.585_DP +  1554.908_DP * em -  2366.899_DP * emsq &
           & +  1215.972_DP * eoc
        g410 = -1052.797_DP +  4758.686_DP * em -  7193.992_DP * emsq &
           & +  3651.957_DP * eoc
        g422 = -3581.690_DP + 16178.110_DP * em - 24462.770_DP * emsq &
           & + 12422.520_DP * eoc
        if (em > 0.715_DP) then
          g520 = -5149.66_DP + 29936.92_DP * em - 54087.36_DP * emsq &
             & + 31324.56_DP * eoc
        else
          g520 =  1464.74_DP -  4664.75_DP * em +  3763.64_DP * emsq
        end if
      end if

      if (em < 0.7_DP) then
        g533 = -919.22770_DP + 4988.6100_DP * em - 9064.7700_DP * emsq &
           & + 5542.21_DP  * eoc
        g521 = -822.71072_DP + 4568.6173_DP * em - 8491.4146_DP * emsq &
           & + 5337.524_DP * eoc
        g532 = -853.66600_DP + 4690.2500_DP * em - 8624.7700_DP * emsq &
           & + 5341.4_DP  * eoc
      else
        g533 = -37995.780_DP + 161616.52_DP * em - 229838.20_DP * emsq &
           & + 109377.94_DP * eoc
        g521 = -51752.104_DP + 218913.95_DP * em - 309468.16_DP * emsq &
           & + 146349.42_DP * eoc
        g532 = -40023.880_DP + 170470.89_DP * em - 242699.48_DP * emsq &
           & + 115605.82_DP * eoc
      end if

      sini2 = sinim * sinim
      f220 = 0.75_DP * (1.0_DP + 2.0_DP * cosim + cosisq)
      f221 = 1.5_DP * sini2
      f321 = 1.875_DP * sinim &
         & * (1.0_DP - 2.0_DP * cosim - 3.0_DP * cosisq)
      f322 = -1.875_DP * sinim &
         & * (1.0_DP + 2.0_DP * cosim - 3.0_DP * cosisq)
      f441 = 35.0_DP * sini2 * f220
      f442 = 39.3750_DP * sini2 * sini2
      f522 = 9.84375_DP * sinim &
         & * (sini2 * (1.0_DP - 2.0_DP * cosim - 5.0_DP * cosisq) &
         & + 0.33333333_DP * (-2.0_DP + 4.0_DP * cosim + 6.0_DP * cosisq))
      f523 = sinim * (4.92187512_DP * sini2 &
         & * (-2.0_DP - 4.0_DP * cosim + 10.0_DP * cosisq) &
         & + 6.56250012_DP * (1.0_DP + 2.0_DP * cosim - 3.0_DP * cosisq))
      f542 = 29.53125_DP * sinim * (2.0_DP - 8.0_DP * cosim &
         & + cosisq * (-12.0_DP + 8.0_DP * cosim + 10.0_DP * cosisq))
      f543 = 29.53125_DP * sinim * (-2.0_DP - 8.0_DP * cosim &
         & + cosisq * (12.0_DP + 8.0_DP * cosim - 10.0_DP * cosisq))
      xno2  = nm * nm
      ainv2 = aonv * aonv
      temp1 = 3.0_DP * xno2 * ainv2
      temp  = temp1 * root22
      d2201 = temp * f220 * g201
      d2211 = temp * f221 * g211
      temp1 = temp1 * aonv
      temp  = temp1 * root32
      d3210 = temp * f321 * g310
      d3222 = temp * f322 * g322
      temp1 = temp1 * aonv
      temp  = 2.0_DP * temp1 * root44
      d4410 = temp * f441 * g410
      d4422 = temp * f442 * g422
      temp1 = temp1 * aonv
      temp  = temp1 * root52
      d5220 = temp * f522 * g520
      d5232 = temp * f523 * g532
      temp  = 2.0_DP * temp1 * root54
      d5421 = temp * f542 * g521
      d5433 = temp * f543 * g533
      xlamo = mod(mo + nodeo + nodeo-theta - theta, PI2)
      xfact = mdot + dmdt + 2.0_DP * (nodedot + dnodt - rptim) - no
      em    = emo
      emsq  = emsqo

    ! ---------------- synchronous resonance terms --------------
    else if (irez == 1) then
      g200  = 1.0_DP + emsq * (-2.5_DP + 0.8125_DP * emsq)
      g310  = 1.0_DP + 2.0_DP * emsq
      g300  = 1.0_DP + emsq * (-6.0_DP + 6.60937_DP * emsq)
      f220  = 0.75_DP * (1.0_DP + cosim) * (1.0_DP + cosim)
      f311  = 0.9375_DP * sinim * sinim * (1.0_DP + 3.0_DP * cosim) &
          & - 0.75_DP * (1.0_DP + cosim)
      f330  = 1.0_DP + cosim
      f330  = 1.875_DP * f330 * f330 * f330
      del1  = 3.0_DP * nm * nm * aonv * aonv
      del2  = 2.0_DP * del1 * f220 * g200 * q22
      del3  = 3.0_DP * del1 * f330 * g300 * q33 * aonv
      del1  = del1 * f311 * g310 * q31 * aonv
      xlamo = mod(mo + nodeo + argpo - theta, PI2)
      xfact = mdot + xpidot - rptim + dmdt + domdt + dnodt - no
    end if

    ! ------------ for sgp4, initialize the integrator ----------
    xli   = xlamo
    xni   = no
    atime = 0.0_DP
    nm    = no + dndt
  end subroutine dsinit

  ! Subroutine dspace
  ! * this procedure provides deep space contributions to mean elements for
  !   perturbing third body.  these effects have been averaged over one
  !   revolution of the sun and moon.  for earth resonance effects, the
  !   effects have been averaged over no revolutions of the satellite.
  !   (mean motion)
  !
  ! :param(in)    integer(4) irez
  ! :param(in)    real(8)    d2201, d2211, d3210, d3222, d4410, d4422,
  !                          d5220, d5232, d5421, d5433, dedt, del1,
  !                          del2, del3, didt, dmdt, dnodt, domdt,
  !                          argpo, argpdot, t, gsto, xfact, xlamo, no, tc
  ! :param(inout) real(8)    atime, em, argpm, inclm, xli, mm, xni,
  !                          nodem, nm
  ! :param(out)   real(8)    dndt
  subroutine dspace( &
    & irez,                                                           &
    & d2201, d2211, d3210, d3222, d4410, d4422, d5220, d5232, d5421,  &
    & d5433, dedt, del1, del2, del3, didt, dmdt, dnodt, domdt, argpo, &
    & argpdot, t, gsto, xfact, xlamo, no, tc,                         &
    & atime, em, argpm, inclm, xli, mm, xni, nodem, nm,               &
    & dndt                                                            &
  & )
    implicit none
    integer(SP), intent(in)    :: irez
    real(DP),    intent(in)    :: d2201, d2211, d3210, d3222, d4410
    real(DP),    intent(in)    :: d4422, d5220, d5232, d5421, d5433
    real(DP),    intent(in)    :: dedt, del1, del2, del3, didt, dmdt
    real(DP),    intent(in)    :: dnodt, domdt, argpo, argpdot
    real(DP),    intent(in)    :: t, tc, gsto, xfact, xlamo, no
    real(DP),    intent(inout) :: atime, em, argpm, inclm, xli, mm
    real(DP),    intent(inout) :: xni, nodem, nm
    real(DP),    intent(out)   :: dndt
    integer(SP) :: iret, iretn
    real(DP)    :: fasx2, fasx4, fasx6, g22, g32, g44, g52, g54
    real(DP)    :: rptim, stepp, stepn, step2, theta, ft, delt
    real(DP)    :: x2li, xomi, x2omi, xl, xldot, xnddt, xndt

    fasx2 = 0.13130908_DP
    fasx4 = 2.8843198_DP
    fasx6 = 0.37448087_DP
    g22   = 5.7686396_DP
    g32   = 0.95240898_DP
    g44   = 1.8014998_DP
    g52   = 1.0508330_DP
    g54   = 4.4108898_DP
    rptim = 4.37526908801129966e-3  ! equates to 7.29211514668855e-5 rad/sec
    stepp =    720.0_DP
    stepn =   -720.0_DP
    step2 = 259200.0_DP

    ! ----------- calculate deep space resonance effects -----------
    dndt  = 0.0_DP
    theta = mod(gsto + tc * rptim, PI2)
    em    = em + dedt * t

    inclm = inclm + didt * t
    argpm = argpm + domdt * t
    nodem = nodem + dnodt * t
    mm    = mm + dmdt * t

    ! ------------------------- epoch restart ----------------------  */
    ! sgp4fix for propagator problems
    ! the following integration works for negative time steps and periods
    ! the specific changes are unknown because the original code was so convoluted
    !
    ! sgp4fix take out atime = 0.0 and fix for faster operation
    ft = 0.0_DP
    if (irez /= 0) then
      ! sgp4fix streamline check
      if (atime == 0.0_DP .or. t * atime <= 0.0_DP &
        & .or. abs(t) < abs(atime)) then
        atime = 0.0_DP
        xni   = no
        xli   = xlamo
      end if

      ! sgp4fix move check outside loop
      if (t > 0.0_DP) then
        delt = stepp
      else
        delt = stepn
      end if

      iretn = 381  ! added for do loop
      iret  =   0  ! added for loop
      do while (iretn == 381)
        ! ------------------- dot terms calculated -------------
        ! ----------- near - synchronous resonance terms -------
        if (irez /= 2) then
          xndt  = del1 * sin(xli - fasx2) &
              & + del2 * sin(2.0_DP * (xli - fasx4)) &
              & + del3 * sin(3.0_DP * (xli - fasx6))
          xldot = xni + xfact
          xnddt = del1 * cos(xli - fasx2) &
              & + 2.0_DP * del2 * cos(2.0_DP * (xli - fasx4)) &
              & + 3.0_DP * del3 * cos(3.0_DP * (xli - fasx6))
          xnddt = xnddt * xldot
        else
          ! --------- near - half-day resonance terms --------
          xomi  = argpo + argpdot * atime
          x2omi = xomi + xomi
          x2li  = xli + xli
          xndt  = (d2201 * sin(x2omi + xli  - g22) &
              & +  d2211 * sin(        xli  - g22) &
              & +  d3210 * sin( xomi + xli  - g32) &
              & +  d3222 * sin(-xomi + xli  - g32) &
              & +  d4410 * sin(x2omi + x2li - g44) &
              & +  d4422 * sin(        x2li - g44) &
              & +  d5220 * sin( xomi + xli  - g52) &
              & +  d5232 * sin(-xomi + xli  - g52) &
              & +  d5421 * sin( xomi + x2li - g54) &
              & +  d5433 * sin(-xomi + x2li - g54))
          xldot = xni + xfact
          xnddt = (d2201 * cos(x2omi + xli  - g22) &
              & +  d2211 * cos(        xli  - g22) &
              & +  d3210 * cos( xomi + xli  - g32) &
              & +  d3222 * cos(-xomi + xli  - g32) &
              & +  d5220 * cos( xomi + xli  - g52) &
              & +  d5232 * cos(-xomi + xli  - g52) &
              & +  2.0_DP &
              & * (d4410 * cos(x2omi + x2li - g44) &
              & +  d4422 * cos(       x2li  - g44) &
              & +  d5421 * cos( xomi + x2li - g54) &
              & +  d5433 * cos(-xomi + x2li - g54)))
          xnddt = xnddt * xldot
        end if

        ! ----------------------- integrator -------------------
        ! sgp4fix move end checks to end of routine
        if (abs(t - atime) >= stepp) then
          iret  = 0
          iretn = 381
        else
          ft    = t - atime
          iretn = 0
        end if

        if (iretn == 381) then
          xli   = xli + xldot * delt + xndt * step2
          xni   = xni + xndt * delt + xnddt * step2
          atime = atime + delt
        end if
      end do

      nm = xni + xndt  * ft + xnddt * ft * ft * 0.5_DP
      xl = xli + xldot * ft + xndt  * ft * ft * 0.5_DP
      if (irez /= 1) then
        mm   = xl - 2.0_DP * nodem + 2.0_DP * theta
        dndt = nm - no
      else
        mm   = xl - nodem - argpm + theta
        dndt = nm - no
      end if

      nm = no + dndt
    end if
  end subroutine dspace

  ! Subroutine dpper
  ! * this procedure provides deep space long period periodic contributions
  !   to the mean elements.  by design, these periodics are zero at epoch.
  !   this used to be dscom which included initialization, but it's really a
  !   recurring function.
  !
  ! :param(in)    type(t_sat)   satrec
  ! :param(in)    real(8)        inclo: inclination - needed for lyddane modification
  ! :param(in)    character(1)    init
  ! :param(in)    character(1) opsmode
  ! :param(inout) real(8)      ep, xincp, nodep, argpp, mp
  subroutine dpper( &
    & satrec, inclo, init, opsmode, ep, xincp, nodep, argpp, mp &
  & )
    implicit none
    type(t_sat),  intent(in)    :: satrec
    real(DP),     intent(in)    :: inclo
    character(1), intent(in)    :: init, opsmode
    real(DP),     intent(inout) :: ep, xincp, nodep, argpp, mp
    real(DP) :: e3, ee2, peo, pgho, pho, pinco, plo, se2, se3, sgh2
    real(DP) :: sgh3, sgh4, sh2, sh3, si2, si3, sl2, sl3, sl4, t
    real(DP) :: xgh2, xgh3, xgh4, xh2, xh3, xi2, xi3, xl2, xl3, xl4
    real(DP) :: xls, xnoh, zmol, zmos, zns, zes, znl, zel, zm, f2, f3
    real(DP) :: sel, ses, sghl, sghs, shll, shs, sil, sinzf, sis, sll
    real(DP) :: sls, pe, pgh, ph, pinc, pl, zf, cosip, cosop, sinip
    real(DP) :: sinop, inclp, alfdp, betdp, dalf, dbet

    ! Copy satellite attributes into local variables for convenience
    ! and symmetry in writing formulae.
    e3    = satrec%e3
    ee2   = satrec%ee2
    peo   = satrec%peo
    pgho  = satrec%pgho
    pho   = satrec%pho
    pinco = satrec%pinco
    plo   = satrec%plo
    se2   = satrec%se2
    se3   = satrec%se3
    sgh2  = satrec%sgh2
    sgh3  = satrec%sgh3
    sgh4  = satrec%sgh4
    sh2   = satrec%sh2
    sh3   = satrec%sh3
    si2   = satrec%si2
    si3   = satrec%si3
    sl2   = satrec%sl2
    sl3   = satrec%sl3
    sl4   = satrec%sl4
    t     = satrec%t
    xgh2  = satrec%xgh2
    xgh3  = satrec%xgh3
    xgh4  = satrec%xgh4
    xh2   = satrec%xh2
    xh3   = satrec%xh3
    xi2   = satrec%xi2
    xi3   = satrec%xi3
    xl2   = satrec%xl2
    xl3   = satrec%xl3
    xl4   = satrec%xl4
    zmol  = satrec%zmol
    zmos  = satrec%zmos

    ! ---------------------- constants -----------------------------
    zns = 1.19459e-5_DP
    zes = 0.01675_DP
    znl = 1.5835218e-4_DP
    zel = 0.05490_DP

    ! --------------- calculate time varying periodics -----------
    zm    = zmos + zns * t
    ! be sure that the initial call has time set to zero
    if (init == "y") zm = zmos
    zf    = zm + 2.0_DP * zes * sin(zm)
    sinzf = sin(zf)
    f2    =  0.5_DP * sinzf * sinzf - 0.25_DP
    f3    = -0.5_DP * sinzf * cos(zf)
    ses   = se2* f2 + se3 * f3
    sis   = si2 * f2 + si3 * f3
    sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf
    sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
    shs   = sh2 * f2 + sh3 * f3
    zm    = zmol + znl * t
    if (init == "y") zm = zmol
    zf    = zm + 2.0_DP * zel * sin(zm)
    sinzf = sin(zf)
    f2    =  0.5_DP * sinzf * sinzf - 0.25_DP
    f3    = -0.5_DP * sinzf * cos(zf)
    sel   = ee2 * f2 + e3 * f3
    sil   = xi2 * f2 + xi3 * f3
    sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf
    sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
    shll  = xh2 * f2 + xh3 * f3
    pe    = ses + sel
    pinc  = sis + sil
    pl    = sls + sll
    pgh   = sghs + sghl
    ph    = shs + shll

    if (init == "n") then
      pe    = pe - peo
      pinc  = pinc - pinco
      pl    = pl - plo
      pgh   = pgh - pgho
      ph    = ph - pho
      inclp = inclp + pinc
      ep    = ep + pe
      sinip = sin(inclp)
      cosip = cos(inclp)

      ! ----------------- apply periodics directly ------------ */
      ! sgp4fix for lyddane choice
      ! strn3 used original inclination - this is technically feasible
      ! gsfc used perturbed inclination - also technically feasible
      ! probably best to readjust the 0.2 limit value and limit discontinuity
      ! 0.2 rad = 11.45916 deg
      ! use next line for original strn3 approach and original inclination
      ! if (inclo >= 0.2)
      ! use next line for gsfc version and perturbed inclination
      if (inclp >= 0.2_DP) then
        ph    = ph / sinip
        pgh   = pgh - cosip * ph
        argpp = argpp + pgh
        nodep = nodep + ph
        mp    = mp + pl
      else
        ! ---- apply periodics with lyddane modification ----
        sinop  = sin(nodep)
        cosop  = cos(nodep)
        alfdp  = sinip * sinop
        betdp  = sinip * cosop
        dalf   =  ph * cosop + pinc * cosip * sinop
        dbet   = -ph * sinop + pinc * cosip * cosop
        alfdp  = alfdp + dalf
        betdp  = betdp + dbet
        nodep  = mod(nodep, PI2)
        ! sgp4fix for afspc written intrinsic functions
        ! nodep used without a trigonometric function ahead
        if (nodep < 0.0_DP .and. opsmode == "a") then
          nodep = nodep + PI2
        end if
        xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep
        xnoh  = nodep
        nodep = atan2(alfdp, betdp)
        ! sgp4fix for afspc written intrinsic functions
        ! nodep used without a trigonometric function ahead
        if (nodep < 0.0_DP .and. opsmode == "a") nodep = nodep + PI2
        if (abs(xnoh - nodep) > PI) then
          if (nodep < xnoh) then
            nodep = nodep + PI2
          else
            nodep = nodep - PI2
          end if
        end if
        mp = mp + pl
        argpp = xls - mp - cosip * nodep
      end if
    end if
  end subroutine dpper

  ! Function gstime
  ! * this function finds the greenwich sidereal time.
  !
  ! :param(in) real(8)  jdut1: julian date in ut1  (days from 4713 bc)
  ! :return    real(8) gstime: greenwich sidereal time  (0 to 2pi rad)
  real(8) function gstime(jdut1)
    implicit none
    real(DP), intent(in) :: jdut1
    real(DP) :: tut1

    tut1 = (jdut1 - 2451545.0_DP) / 36525.0_DP
    gstime = 67310.54841_DP &
         & + ((876600.0_DP * 3600.0_DP + 8640184.812866_DP) &
         & + (0.093104_DP &
         & - 6.2e-6_DP &
         & * tut1) * tut1) * tut1
    gstime = mod(gstime * DEG2RAD / 240.0_DP, PI2)
    if (gstime < 0.0_DP) gstime = gstime + PI2
  end function gstime

  ! Subroutine sgp4
  ! * this procedure is the sgp4 prediction model from space command. this is an
  !   updated and combined version of sgp4 and sdp4, which were originally
  !   published separately in spacetrack report #3. this version follows the
  !   methodology from the aiaa paper (2006) describing the history and
  !   development of the code.
  !
  ! :param(in)    type(t_cst) gravconst: constants
  ! :param(in)    real(8)        tsince: time eince epoch (minutes)
  ! :param(inout) type(t_sat)    satrec: initialised structure from sgp4init call.
  ! :param(out)   real(8)          r(3): position vector  (km)
  ! :param(out)   real(8)          v(3): velocity         (km/sec)
  subroutine sgp4(gravconst, tsince, satrec, r, v)
    implicit none
    type(t_cst), intent(in)    :: gravconst
    real(DP),    intent(in)    :: tsince
    type(t_sat), intent(inout) :: satrec
    real(DP),    intent(out)   :: r(3), v(3)
    integer(SP) :: ktr
    real(DP)    :: mrt, temp, temp4, tempa, tempe, templ, x2o3
    real(DP)    :: vkmpersec, xmdf, argpdf, argpm, am, em, mm, nm
    real(DP)    :: nodedf, nodem, delomg, delm, delmtemp, inclm
    real(DP)    :: t2, t3, t4, tc, dndt, emsq, xlm, argpp, ep, mp
    real(DP)    :: nodep, xincp, sinim, sinip, cosim, cosip
    real(DP)    :: axnl, aynl, xl, eo1, tem5, u, coseo1, sineo1
    real(DP)    :: ecose, esine, el2, pl, betal, rdotl, rl, rvdotl
    real(DP)    :: cosu, cos2u, sinu, sin2u, su, temp1, temp2, cosisq
    real(DP)    :: mvt, rvdot, xinc, xnode, ux, uy, uz, vx, vy, vz
    real(DP)    :: cnod, cosi, cossu, sini, sinsu, snod, xmx, xmy
    real(DP)    :: mr

    satrec%error = 0
    mrt = 0.0_DP

    ! ------------------ set mathematical constants --------------- */
    ! sgp4fix divisor for divide by zero check on inclination
    ! the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    ! 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    temp4 = 1.5e-12_DP
    x2o3  = 2.0_DP / 3.0_DP
    ! sgp4fix identify constants and allow alternate values
    vkmpersec = gravconst%radiusearthkm * gravconst%xke / 60.0_DP

    ! --------------------- clear sgp4 error flag -----------------
    satrec%t = tsince

    ! ------- update for secular gravity and atmospheric drag -----
    xmdf   = satrec%mo + satrec%mdot * satrec%t
    argpdf = satrec%argpo + satrec%argpdot * satrec%t
    nodedf = satrec%nodeo + satrec%nodedot * satrec%t
    argpm  = argpdf
    mm     = xmdf
    t2     = satrec%t * satrec%t
    nodem  = nodedf + satrec%nodecf * t2
    tempa  = 1.0_DP - satrec%cc1 * satrec%t
    tempe  = satrec%bstar * satrec%cc4 * satrec%t
    templ  = satrec%t2cof * t2

    if (satrec%isimp /= 1) then
      delomg = satrec%omgcof * satrec%t
      ! sgp4fix use mutliply for speed instead of pow
      delmtemp = 1.0_DP + satrec%eta * cos(xmdf)
      delm   = satrec%xmcof &
           & * (delmtemp * delmtemp * delmtemp - satrec%delmo)
      temp   = delomg + delm
      mm     = xmdf + temp
      argpm  = argpdf - temp
      t3     = t2 * satrec%t
      t4     = t3 * satrec%t
      tempa  = tempa - satrec%d2 * t2 - satrec%d3 * t3 &
                   & - satrec%d4 * t4
      tempe  = tempe + satrec%bstar * satrec%cc5 * (sin(mm) &
                   & - satrec%sinmao)
      templ  = templ + satrec%t3cof * t3 + t4 * (satrec%t4cof &
                   & + satrec%t * satrec%t5cof)
    end if

    nm    = satrec%no
    em    = satrec%ecco
    inclm = satrec%inclo
    if (satrec%method == "d") then
      tc = satrec%t
      call dspace( &
        ! (in)
        & satrec%irez, &
        & satrec%d2201, satrec%d2211, satrec%d3210, satrec%d3222,   &
        & satrec%d4410, satrec%d4422, satrec%d5220, satrec%d5232,   &
        & satrec%d5421, satrec%d5433, satrec%dedt,  satrec%del1,    &
        & satrec%del2,  satrec%del3,  satrec%didt,  satrec%dmdt,    &
        & satrec%dnodt, satrec%domdt, satrec%argpo, satrec%argpdot, &
        & satrec%t,     satrec%gsto,  satrec%xfact, satrec%xlamo,   &
        & satrec%no,    tc, &
        ! (inout)
        & satrec%atime, em, argpm, inclm, satrec%xli, &
        & mm, satrec%xni, nodem, nm, &
        ! (out)
        & dndt &
      & )
    end if

    if (nm <= 0.0_DP) then
      satrec%error = 2
      return
    end if

    ! mean motion less than 0.0
    am = (gravconst%xke / nm)**x2o3 * tempa * tempa
    nm = gravconst%xke / (am**1.5_DP)
    em = em - tempe

    ! fix tolerance for error recognition
    if (em >= 1.0_DP .or. em < -0.001_DP .or. am < 0.95_DP) then
      satrec%error = 1
      return
    end if

    ! sgp4fix fix tolerance to avoid a divide by zero
    if (em < 1.0e-6_DP) em = 1.0e-6_DP
    mm   = mm + satrec%no * templ
    xlm  = mm + argpm + nodem
    emsq = em * em
    temp = 1.0_DP - emsq

    nodem = mod(nodem, PI2)
    argpm = mod(argpm, PI2)
    xlm   = mod(xlm, PI2)
    mm    = mod(xlm - argpm - nodem, PI2)
    do while (mm < 0.0_DP)
      mm = mm + PI2
    end do

    ! ----------------- compute extra mean quantities -------------
    sinim = sin(inclm)
    cosim = cos(inclm)

    ! -------------------- add lunar-solar periodics --------------
    ep    = em
    xincp = inclm
    argpp = argpm
    nodep = nodem
    mp    = mm
    sinip = sinim
    cosip = cosim
    if (satrec%method == "d") then
      call dpper( &
        ! (in)
        & satrec, satrec%inclo, "n", satrec%opsmode, &
        ! (inout)
        & ep, xincp, nodep, argpp, mp &
      & )

      if (xincp < 0.0_DP) then
        xincp = -xincp
        nodep = nodep + pi
        argpp = argpp - pi
      end if

      if (ep < 0.0_DP .or. ep > 1.0_DP) then
        satrec%error = 3
        return
      end if
    end if

    ! -------------------- long period periodics ------------------
    if (satrec%method == "d") then
      sinip = sin(xincp)
      cosip = cos(xincp)
      satrec%aycof = -0.5_DP * gravconst%j3oj2 * sinip
      ! sgp4fix for divide by zero for xincp = 180 deg
      if (abs(cosip + 1.0_DP) > 1.5e-12_DP) then
        satrec%xlcof = -0.25_DP * gravconst%j3oj2 * sinip &
                   & * (3.0_DP + 5.0_DP * cosip) / (1.0_DP + cosip)
      else
        satrec%xlcof = -0.25_DP * gravconst%j3oj2 * sinip &
                   & * (3.0_DP + 5.0_DP * cosip) / temp4
      end if
    end if

    axnl = ep * cos(argpp)
    temp = 1.0_DP / (am * (1.0_DP - ep * ep))
    aynl = ep * sin(argpp) + temp * satrec%aycof
    xl   = mp + argpp + nodep + temp * satrec%xlcof * axnl

    ! --------------------- solve kepler's equation ---------------
    u    = mod(xl - nodep, PI2)
    eo1  = u
    tem5 = 9999.9_DP
    ktr  = 1
    ! sgp4fix for kepler iteration
    ! the following iteration needs better limits on corrections
    do while (abs(tem5) >= 1.0e-12_DP .and. ktr <= 10)
      sineo1 = sin(eo1)
      coseo1 = cos(eo1)
      tem5   = 1.0_DP - coseo1 * axnl - sineo1 * aynl
      tem5   = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5
      if (abs(tem5) >= 0.95_DP) then
        if (tem5 > 0.0_DP) then
          tem5 = 0.95_DP
        else
          tem5 = -0.95_DP
        end if
      end if
      eo1 = eo1 + tem5
      ktr = ktr + 1
    end do

    ! ------------- short period preliminary quantities -----------
    ecose = axnl * coseo1 + aynl * sineo1
    esine = axnl * sineo1 - aynl * coseo1
    el2   = axnl * axnl + aynl * aynl
    pl    = am * (1.0_DP - el2)
    if (pl < 0.0_DP) then
      satrec%error = 4
      return
    else
      rl     = am * (1.0_DP - ecose)
      rdotl  = sqrt(am) * esine/rl
      rvdotl = sqrt(pl) / rl
      betal  = sqrt(1.0_DP - el2)
      temp   = esine / (1.0_DP + betal)
      sinu   = am / rl * (sineo1 - aynl - axnl * temp)
      cosu   = am / rl * (coseo1 - axnl + aynl * temp)
      su     = atan2(sinu, cosu)
      sin2u  = (cosu + cosu) * sinu
      cos2u  = 1.0_DP - 2.0_DP * sinu * sinu
      temp   = 1.0_DP / pl
      temp1  = 0.5_DP * gravconst%j2 * temp
      temp2  = temp1 * temp

      ! -------------- update for short period periodics ------------
      if (satrec%method == "d") then
        cosisq = cosip * cosip
        satrec%con41  = 3.0_DP * cosisq - 1.0_DP
        satrec%x1mth2 = 1.0_DP  - cosisq
        satrec%x7thm1 = 7.0_DP * cosisq - 1.0_DP
      end if

      mrt   = rl * (1.0_DP - 1.5_DP * temp2 * betal * satrec%con41) &
          & + 0.5_DP * temp1 * satrec%x1mth2 * cos2u
      su    = su - 0.25_DP * temp2 * satrec%x7thm1 * sin2u
      xnode = nodep + 1.5_DP * temp2 * cosip * sin2u
      xinc  = xincp + 1.5_DP * temp2 * cosip * sinip * cos2u
      mvt   = rdotl - nm * temp1 * satrec%x1mth2 * sin2u / gravconst%xke
      rvdot = rvdotl + nm * temp1 * (satrec%x1mth2 * cos2u &
          & + 1.5_DP * satrec%con41) / gravconst%xke

      ! --------------------- orientation vectors -------------------
      sinsu =  sin(su)
      cossu =  cos(su)
      snod  =  sin(xnode)
      cnod  =  cos(xnode)
      sini  =  sin(xinc)
      cosi  =  cos(xinc)
      xmx   = -snod * cosi
      xmy   =  cnod * cosi
      ux    =  xmx * sinsu + cnod * cossu
      uy    =  xmy * sinsu + snod * cossu
      uz    =  sini * sinsu
      vx    =  xmx * cossu - cnod * sinsu
      vy    =  xmy * cossu - snod * sinsu
      vz    =  sini * cossu

      ! --------- position and velocity (in km and km/sec) ----------
      mr = mrt * gravconst%radiusearthkm
      r = (/mr * ux, mr * uy, mr * uz/)
      v = (/ &
        & (mvt * ux + rvdot * vx) * vkmpersec, &
        & (mvt * uy + rvdot * vy) * vkmpersec, &
        & (mvt * uz + rvdot * vz) * vkmpersec  &
      & /)
    end if

    ! sgp4fix for decaying satellites
    if (mrt < 1.0) then
      satrec%error = 6
      return
    end if
  end subroutine sgp4

  ! ============================
  ! Public subroutines/functions
  ! ============================

  ! Return a position and velocity vector for a given date and time.
  !
  ! :param(in)    type(t_cst)  gravconst: constants
  ! :param(inout) type(t_sat)     satrec: common values for subsequent calls
  ! :param(in)    type(t_time)       ut1: UT1
  ! :param(out)   real(8)           r(3): positions
  ! :param(out)   real(8)           v(3): velocities
  subroutine propagate(gravconst, satrec, ut1, r, v)
    implicit none
    type(t_cst),  intent(in)    :: gravconst
    type(t_sat),  intent(inout) :: satrec
    type(t_time), intent(in)    :: ut1
    real(DP),     intent(out)   :: r(3), v(3)
    real(DP) :: j, m

    j = jday(ut1%year, ut1%month, ut1%day, &
      & ut1%hour, ut1%minute, &
      & real(ut1%second, DP) + ut1%msecond * 1.0e-3)
    m = (j - satrec%jdsatepoch) * MIN_D
    call sgp4(gravconst, m, satrec, r, v)
  end subroutine propagate

  ! Subroutine sgp4init
  ! * this procedure initializes variables for sgp4.
  !
  ! :param(in)    type(t_cst) gravconst: constants
  ! :param(in)    character(1)  opsmode: mode of operation afspc or improved 'a', 'i'
  ! :param(in)    integer(4)       satn: satellite number
  ! :param(in)    real(8)         epoch: epoch time in days from jan 0, 1950. 0 hr
  ! :param(in)    real(8)        xbstar: sgp4 type drag coefficient              kg/m2er
  ! :param(in)    real(8)         xecco: eccentricity
  ! :param(in)    real(8)        xargpo: argument of perigee (output if ds)
  ! :param(in)    real(8)        xinclo: inclination
  ! :param(in)    real(8)           xmo: mean anomaly (output if ds)
  ! :param(in)    real(8)           xno: mean motion
  ! :param(in)    real(8)        xnodeo: right ascension of ascending node
  ! :param(inout) type(t_sat)    satrec: common values for subsequent calls
  subroutine sgp4init( &
    & gravconst, opsmode, satn, epoch,                        &
    & xbstar, xecco, xargpo, xinclo, xmo, xno, xnodeo, satrec &
  & )
    implicit none
    type(t_cst),  intent(in)    :: gravconst
    character(1), intent(in)    :: opsmode
    integer(SP),  intent(in)    :: satn
    real(DP),     intent(in)    :: epoch, xbstar, xecco, xargpo
    real(DP),     intent(in)    :: xinclo, xmo, xno, xnodeo
    type(t_sat),  intent(inout) :: satrec
    integer(SP) :: error
    real(DP)    :: ss, qzms2ttemp, qzms2t, x2o3, ainv, ao, con42
    real(DP)    :: cosio, cosio2, cc1sq, eccsq, omeosq, posq, rp
    real(DP)    :: rteosq, sinio, perige, qzms24, sfour, qzms24temp
    real(DP)    :: pinvsq, coef, coef1, eeta, etasq, psisq, tsi
    real(DP)    :: cc2, cc3, cosio4, temp, temp1, temp2, temp3, temp4
    real(DP)    :: xhdot1, xpidot, delmotemp, tc, inclm, argpm, nodem
    real(DP)    :: mm, dndt, snodm, cnodm, sinim, cosim, sinomm, cosomm
    real(DP)    :: day, em, emsq, gam, rtemsq
    real(DP)    :: s1, s2, s3, s4, s5, s6, s7, ss1, ss2, ss3, ss4, ss5
    real(DP)    :: ss6, ss7, sz1, sz2, sz3, sz11, sz12, sz13, sz21, sz22
    real(DP)    :: sz23, sz31, sz32, sz33, nm, z1, z2, z3, z11, z12
    real(DP)    :: z13, z21, z22, z23, z31, z32, z33
    real(DP)    :: r(3), v(3)

    ! ------------------------ initialization --------------------- */
    ! sgp4fix divisor for divide by zero check on inclination
    ! the old check used 1.0 + cos(pi-1.0e-9), but then compared it to
    ! 1.5 e-12, so the threshold was changed to 1.5e-12 for consistency
    temp4 = 1.5e-12_DP

    ! sgp4fix - note the following variables are also passed directly via satrec%
    ! it is possible to streamline the sgp4init call by deleting the "x"
    ! variables, but the user would need to set the satrec%* values first. we
    ! include the additional assignments in case twoline2rv is not used.
    satrec%bstar = xbstar
    satrec%ecco  = xecco
    satrec%argpo = xargpo
    satrec%inclo = xinclo
    satrec%mo    = xmo
    satrec%no    = xno
    satrec%nodeo = xnodeo

    ! sgp4fix add opsmode
    satrec%opsmode = opsmode

    ! ------------------------ earth constants -----------------------
    ! sgp4fix identify constants and allow alternate values
    ss = 78.0_DP / gravconst%radiusearthkm + 1.0_DP
    ! sgp4fix use multiply for speed instead of pow
    qzms2ttemp = (120.0_DP - 78.0_DP) / gravconst%radiusearthkm
    qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp
    x2o3   = 2.0_DP / 3.0_DP

    satrec%init = "y"
    satrec%t    = 0.0_DP

    call initl( &
      & satn, gravconst, satrec%ecco, epoch, satrec%inclo,  &
      & satrec%opsmode, satrec%no, satrec%method,           &
      & ainv, ao, satrec%con41, con42, cosio, cosio2,       &
      & eccsq, omeosq, posq, rp, rteosq, sinio, satrec%gsto &
    & )

    ! sgp4fix remove this check as it is unnecessary
    ! the mrt check in sgp4 handles decaying satellite cases even if the starting
    ! condition is below the surface of te earth
    if (rp < 1.0_DP) then
      satrec%error = 5
      return
    end if

    if (omeosq >= 0.0_DP .or. satrec%no >= 0.0_DP) then
      satrec%isimp = 0
      if (rp < 220.0_DP / gravconst%radiusearthkm + 1.0_DP) then
        satrec%isimp = 1
      end if
      sfour  = ss
      qzms24 = qzms2t
      perige = (rp - 1.0_DP) * gravconst%radiusearthkm

      ! for perigees below 156 km, s and qoms2t are altered -
      if (perige < 156.0_DP) then
        sfour = perige - 78.0
        if (perige < 98.0_DP) sfour = 20.0_DP
        ! sgp4fix use multiply for speed instead of pow
        qzms24temp = (120.0_DP - sfour) / gravconst%radiusearthkm
        qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp
        sfour  = sfour / gravconst%radiusearthkm + 1.0_DP
      end if

      pinvsq = 1.0_DP / posq
      tsi    = 1.0_DP / (ao - sfour)
      satrec%eta = ao * satrec%ecco * tsi
      etasq  = satrec%eta * satrec%eta
      eeta   = satrec%ecco * satrec%eta
      psisq  = abs(1.0_DP - etasq)
      coef   = qzms24 * tsi**4.0_DP
      coef1  = coef / psisq**3.5_DP
      cc2    = coef1 * satrec%no &
           & * (ao * (1.0_DP + 1.5_DP * etasq + eeta * (4.0_DP + etasq)) &
           & + 0.375_DP * gravconst%j2 * tsi / psisq * satrec%con41 &
           & * (8.0_DP + 3.0_DP * etasq * (8.0_DP + etasq)))
      satrec%cc1 = satrec%bstar * cc2
      cc3    = 0.0_DP
      if (satrec%ecco > 1.0e-4_DP) then
        cc3 = -2.0_DP * coef * tsi * gravconst%j3oj2 * satrec%no * sinio &
          & / satrec%ecco
      end if
      satrec%x1mth2 = 1.0_DP - cosio2
      satrec%cc4 = 2.0_DP * satrec%no * coef1 * ao * omeosq &
               & * (satrec%eta * (2.0_DP + 0.5_DP * etasq) &
               & + satrec%ecco * (0.5_DP + 2.0_DP * etasq) &
               & - gravconst%j2 * tsi / (ao * psisq) &
               & * (-3.0_DP * satrec%con41 * (1.0_DP - 2.0_DP * eeta &
               & + etasq * (1.5_DP - 0.5_DP * eeta)) &
               & + 0.75_DP * satrec%x1mth2 * (2.0_DP * etasq &
               & - eeta * (1.0_DP + etasq)) * cos(2.0_DP * satrec%argpo)))
      satrec%cc5 = 2.0_DP * coef1 * ao * omeosq &
               & * (1.0_DP + 2.75_DP * (etasq + eeta) + eeta * etasq)
      cosio4 = cosio2 * cosio2
      temp1  = 1.5_DP * gravconst%j2 * pinvsq * satrec%no
      temp2  = 0.5_DP * temp1 * gravconst%j2 * pinvsq
      temp3  = -0.46875_DP * gravconst%j4 * pinvsq * pinvsq * satrec%no
      satrec%mdot    = satrec%no + 0.5_DP * temp1 * rteosq * satrec%con41 &
                   & + 0.0625_DP * temp2 * rteosq &
                   & * (13.0_DP - 78.0_DP * cosio2 + 137.0_DP * cosio4)
      satrec%argpdot = (-0.5_DP * temp1 * con42 + 0.0625_DP * temp2 &
                   & * (7.0_DP - 114.0_DP * cosio2 + 395.0_DP * cosio4) &
                   & + temp3 * (3.0_DP - 36.0_DP * cosio2 + 49.0_DP * cosio4))
      xhdot1 = -temp1 * cosio
      satrec%nodedot = xhdot1 + (0.5_DP * temp2 * (4.0_DP - 19.0_DP * cosio2) &
                   & + 2.0_DP * temp3 * (3.0_DP - 7.0_DP * cosio2)) * cosio
      xpidot = satrec%argpdot + satrec%nodedot
      satrec%omgcof = satrec%bstar * cc3 * cos(satrec%argpo)
      satrec%xmcof  = 0.0_DP
      if (satrec%ecco > 1.0e-4_DP) then
        satrec%xmcof = -x2o3 * coef * satrec%bstar / eeta
      end if
      satrec%nodecf = 3.5_DP * omeosq * xhdot1 * satrec%cc1
      satrec%t2cof  = 1.5_DP * satrec%cc1
      ! sgp4fix for divide by zero with xinco = 180 deg
      if (abs(cosio + 1.0_DP) > 1.5e-12_DP) then
        satrec%xlcof = -0.25_DP * gravconst%j3oj2 * sinio &
                   & * (3.0_DP + 5.0_DP * cosio) / (1.0_DP + cosio)
      else
        satrec%xlcof = -0.25_DP * gravconst%j3oj2 * sinio &
                   & * (3.0_DP + 5.0_DP * cosio) / temp4
      end if
      satrec%aycof  = -0.5_DP * gravconst%j3oj2 * sinio
      ! sgp4fix use multiply for speed instead of pow
      delmotemp     = 1.0_DP + satrec%eta * cos(satrec%mo)
      satrec%delmo  = delmotemp * delmotemp * delmotemp
      satrec%sinmao = sin(satrec%mo)
      satrec%x7thm1 = 7.0_DP * cosio2 - 1.0_DP

      ! --------------- deep space initialization -------------
      if (PI2 / satrec%no >= 225.0_DP) then
        satrec%method = "d"
        satrec%isimp  = 1
        tc    =  0.0_DP
        inclm = satrec%inclo

        call dscom( &
          ! (in)
          & epoch, satrec%ecco, satrec%argpo, tc, satrec%inclo, &
          & satrec%nodeo, satrec%no,                            &
          ! (inout)
          & satrec%e3, satrec%ee2, satrec%peo, satrec%pgho,     &
          & satrec%pho, satrec%pinco, satrec%plo,               &
          & satrec%se2, satrec%se3,                             &
          & satrec%sgh2, satrec%sgh3, satrec%sgh4,              &
          & satrec%sh2, satrec%sh3, satrec%si2, satrec%si3,     &
          & satrec%sl2, satrec%sl3, satrec%sl4,                 &
          & satrec%xgh2, satrec%xgh3, satrec%xgh4,              &
          & satrec%xh2, satrec%xh3, satrec%xi2, satrec%xi3,     &
          & satrec%xl2, satrec%xl3, satrec%xl4,                 &
          & satrec%zmol, satrec%zmos,                           &
          ! (out)
          & snodm, cnodm, sinim, cosim, sinomm, cosomm, day,    &
          & em, emsq, gam, rtemsq, s1, s2, s3, s4, s5,          &
          & s6, s7, ss1, ss2, ss3, ss4, ss5, ss6, ss7, sz1,     &
          & sz2, sz3, sz11, sz12, sz13, sz21, sz22, sz23,       &
          & sz31, sz32, sz33, nm, z1, z2, z3,                   &
          & z11, z12, z13, z21, z22, z23, z31, z32, z33         &
        & )

        argpm = 0.0_DP
        nodem = 0.0_DP
        mm    = 0.0_DP

        call dsinit( &
          ! (in)
          & gravconst,                                              &
          & cosim, emsq, satrec%argpo, s1, s2, s3, s4, s5, sinim,   &
          & ss1, ss2, ss3, ss4, ss5, sz1, sz3, sz11, sz13, sz21,    &
          & sz23, sz31, sz33, satrec%t, tc, satrec%gsto, satrec%mo, &
          & satrec%mdot, satrec%no, satrec%nodeo, satrec%nodedot,   &
          & xpidot, z1, z3, z11, z13, z21, z23, z31, z33,           &
          & satrec%ecco, eccsq, &
          ! (inout)
          & em, argpm, inclm, mm, nm, nodem, satrec%irez,           &
          & satrec%atime, satrec%d2201, satrec%d2211, satrec%d3210, &
          & satrec%d3222, satrec%d4410, satrec%d4422, satrec%d5220, &
          & satrec%d5232, satrec%d5421, satrec%d5433, satrec%dedt,  &
          & satrec%didt, satrec%dmdt, satrec%dnodt, satrec%domdt,   &
          & satrec%del1, satrec%del2, satrec%del3, satrec%xfact,    &
          & satrec%xlamo, satrec%xli, satrec%xni,                   &
          ! (out)
          & dndt                                                    &
        & )
      end if

      ! ----------- set variables if not deep space -----------
      if (satrec%isimp /= 1) then
        cc1sq        = satrec%cc1 * satrec%cc1
        satrec%d2    = 4.0_DP * ao * tsi * cc1sq
        temp         = satrec%d2 * tsi * satrec%cc1 / 3.0_DP
        satrec%d3    = (17.0_DP * ao + sfour) * temp
        satrec%d4    = 0.5_DP * temp * ao * tsi &
                   & * (221.0_DP * ao + 31.0_DP * sfour) * satrec%cc1
        satrec%t3cof = satrec%d2 + 2.0_DP * cc1sq
        satrec%t4cof = 0.25_DP * (3.0_DP * satrec%d3 + satrec%cc1 &
                   & * (12.0_DP * satrec%d2 + 10.0_DP * cc1sq))
        satrec%t5cof = 0.2_DP * (3.0_DP * satrec%d4 &
                   & + 12.0_DP * satrec%cc1 * satrec%d3 &
                   & + 6.0_DP * satrec%d2 * satrec%d2 &
                   & + 15.0_DP * cc1sq * (2.0_DP * satrec%d2 + cc1sq))
      end if
    end if

    call sgp4(gravconst, 0.0_DP, satrec, r, v)

    satrec%init = "n"
  end subroutine sgp4init
end module propagation



program thermo

    use mod_thermo
    use mod_fuel
    use mod_fluid

    implicit none

    integer ( kind = 4 ) :: option
    real    ( kind = 8 ) :: l
    real    ( kind = 8 ) :: pow
    real    ( kind = 8 ) :: peak
    real    ( kind = 8 ) :: a
    integer ( kind = 4 ), parameter :: lstep = 100
    integer ( kind = 4 ) :: nrods
    real    ( kind = 8 ) :: qp(lstep)

    real    ( kind = 8 ) :: mdot
    real    ( kind = 8 ) :: tin
    real    ( kind = 8 ) :: pin
    integer ( kind = 4 ) :: nassemb
    integer ( kind = 4 ) :: nshape
    real    ( kind = 8 ) :: h(lstep)



    option = 1
    l = 144
    pow = 3300.0D+0
    peak = 1.9D+0
    peak = 2.7D+0
    a = 1.0D+0
    qp(1:lstep) = 0.0D+0
    nrods = 7*7*764

    mdot = 1.025E+08 * 0.90D+0
    tin = 537.0D+0
    pin = 2250.0D+0
    h(1:lstep) = 0.0D+0

    nassemb = 764
    nshape = 7*7


    call fuelshape ( option, l, pow, peak, a, lstep, nrods, qp )

    call enth_rise ( mdot, qp, tin, pin, nrods, l, lstep, nassemb, nshape, h )





end program thermo

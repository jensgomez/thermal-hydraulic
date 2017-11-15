
module mod_thermo

    implicit none

    contains


    subroutine p_sat ( t, psat )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the saturation pressure as a function
!   of temperature (valid for 89.965 deg C < T < 373.253 deg C).
!   Correlations taken from Appendix I
!
!   t,    input,  real ( kind = 8 ), the temperature
!                                    units: deg F
!
!   psat, output, real ( kind = 8 ), the saturation pressure for the given temp.
!                                    units: psia
!

        implicit none

        real ( kind = 8 ), intent ( in )    :: t
        real ( kind = 8 ), intent ( inout ) :: psat


!
!       Variables used in this subroutine
        real ( kind = 8 ) :: tc
        real ( kind = 8 ) :: psat_mpa


!       Convert temperature from deg F to deg C
        tc = (5.0D+0/9.0D+0)*(t - 32.0D+0)

        if ( tc < 89.965D+0 .or. tc > 373.253) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( tc >= 89.965D+0 .and. tc < 139.781D+0 ) then


            psat_mpa = ((tc+57.0D+0)/236.2315D+0)**5.602972D+0


        else if ( tc >= 139.781D+0 .and. tc < 203.662D+0 ) then

            psat_mpa = ((tc+28.0D+0)/207.9248D+0)**4.778504D+0


        else if ( tc >= 203.662D+0 .and. tc < 299.407D+0 ) then

            psat_mpa = ((tc+5.0D+0)/185.0779D+0)**4.304376D+0


        else if ( tc >= 299.407D+0 .and. tc < 355.636D+0 ) then

            psat_mpa = ((tc+5.0D+0)/185.0779D+0)**4.304376D+0


        else if ( tc >= 299.407D+0 .and. tc < 355.636D+0 ) then

            psat_mpa = ((tc+16.0D+0)/195.1819D+0)**4.460843D+0


        else if ( tc >= 355.563D+0 .and. tc <= 373.253D+0 ) then

            psat_mpa = ((tc+50.0D+0)/227.2963D+0)**4.960785D+0


        end if


!       Convert from MPa to psia
        psat = psat_mpa * 145.038D+0


    end subroutine p_sat


    subroutine t_sat ( p, tsat )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the saturation temperature as a function
!   of pressure (valid for 0.070 MPa < P < 21.85 MPa).
!   Correlations taken from Appendix II
!
!   p,    input,  real ( kind = 8 ), the pressure
!                                    units: psia
!
!   tsat, output, real ( kind = 8 ), the saturation temp for the given pres.
!                                    units: deg F
!
!


        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( inout ) :: tsat

!       Local variables
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat_si

!       Convert pressures from psia to MPa
        p_si = p * 0.00689476D+0


!       Correlations from Appendix II
        if ( p_si < 0.070D+0 .and. p_si > 21.85D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si >= 0.070D+0 .and. p_si < 0.359D+0 ) then

            tsat_si = 236.2315D+0 * p_si**(0.1784767D+0) - 57.0D+0


        else if ( p_si >= 0.359D+0 .and. p_si <= 1.676D+0 ) then

            tsat_si = 207.9248D+0 * p_si**(0.2092705D+0) - 28.0D+0


        else if ( p_si > 1.676D+0 .and. p_si <= 8.511D+0 ) then

            tsat_si = 185.0779D+0 * p_si**(0.2323217D+0) - 5.0D+0


        else if ( p_si > 8.511D+0 .and. p_si < 17.69D+0 ) then

            tsat_si = 195.1819D+0 * p_si**(0.2241729D+0) - 16.0D+0


        else if ( p_si >= 17.69D+0 .and. p_si <= 21.85D+0 ) then

            tsat_si = 227.2963D+0 * p_si**(0.201581D+0) - 50.0D+0

        end if


!       Convert the saturation temperature into deg F
        tsat = tsat_si*(9.0D+0/5.0D+0) + 32.0D+0

    end subroutine t_sat


    subroutine liq_dens ( psat, p, t, rho )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid density as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!   psat,  input, real ( kind = 8 ), the saturation pressure
!                                    units: psia
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated density
!                                    units: lbm/ft**3
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: psat
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: rho


!
!       Variables used locally in subroutine
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: rho_si
        real ( kind = 8 ) :: rho_t


!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for rho_t
        if ( psat_si < 0.075D+0 .or. psat_si > 21.500D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( psat_si >= 0.075D+0 .and. psat_si <= 1.00D+0 ) then

            rho_t = ( 1.2746977E-04 * psat_si**(0.4644339D+0) + 0.001D+0 )**-1.0D+0


        else if ( psat_si > 1.00D+0 .and. psat_si <= 3.88D+0 ) then

            rho_t = ( 1.0476071E-04 * psat_si**(0.5651090D+0) + 0.001022D+0 )**-1.0D+0


        else if ( psat_si > 3.880D+0 .and. psat_si <= 8.840D+0 ) then

            rho_t = ( 3.2836717E-05 * psat_si + 1.12174735E-03 )**-1.0D+0


        else if ( psat_si > 8.840D+0 .and. psat_si <= 14.463D+0 ) then

            rho_t = ( 3.3551046E-04 * exp(5.8403566E-02*psat_si) + 0.00085D+0 )**-1.0D+0


        else if ( psat_si > 14.463D+0 .and. psat_si < 18.052D+0 ) then

            rho_t = ( 3.1014626E-08 * psat_si**(3.284754D+0) + 0.001430D+0 )**-1.0D+0


        else if ( psat_si >= 18.052D+0 .and. psat_si < 20.204D+0 ) then

            rho_t = ( 1.5490787E-11 * psat_si**(5.7205D+0) + 0.001605D+0 )**-1.0D+0


        else if ( psat_si >= 20.204D+0 .and. psat_si <= 21.500D+0 ) then

            rho_t = ( 4.1035988E-24 * psat_si**(15.03329D+0) + 0.00189D+0 )**-1.0D+0


        end if


!       Density equation from Section 4.1
        rho_si = rho_t + ((170.0D+0/(375.0D+0 - t_si)) - 0.20D+0)*( p_si - psat_si )

!       Convert from kg/m**3 to lb/ft**3
        rho = rho_si * 0.062428D+0



    end subroutine liq_dens


    subroutine liq_enth ( psat, p, t, h )

! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid enthalpy as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!   psat,  input, real ( kind = 8 ), the saturation pressure
!                                    units: psia
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated enthalpy
!                                    units: BTU/lbm
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: psat
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: h


!       Variables used locally in this subroutine
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: hs



!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for hs
        if ( psat_si < 0.075D+0 .and. psat_si > 21.70D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '

        else if ( psat_si > 0.075D+0 .and. psat_si < 0.942D+0 ) then

            hs = 912.1779D+0 * psat_si**(0.2061637D+0) - 150.0D+0


        else if ( psat_si >= 0.942D+0 .and. psat_si < 4.020D+0 ) then

            hs = 638.0621D+0 * psat_si**(0.2963197D+0) + 125.0D+0


        else if ( psat_si >= 4.020D+0 .and. psat_si < 9.964D+0 ) then

            hs = 373.7665D+0 * psat_si**(0.4235532D+0) + 415.0D+0


        else if ( psat_si >= 9.964D+0 .and. psat_si < 16.673D+0 ) then

            hs = 75.38673D+0 * psat_si**(0.8282384D+0) + 900.0D+0


        else if ( psat_si >= 16.673D+0 .and. psat_si < 20.396D+0 ) then

            hs = 0.1150827D+0 * psat_si**(2.711412D+0) + 1440.0D+0


        else if ( psat_si >= 20.396D+0 .and. psat_si <= 21.700D+0 ) then

            hs = 9.1417257E-14 * psat_si**(11.47287D+0) + 1752.0D+0

        end if


!       Enthalpy equation from Section 4.2
        h = hs + (1.4D+0 - ( 169.0D+0 / (369.0D0 - t_si)))*(p_si - psat_si)


!       Convert from kJ/kg to BTU/lbm
        h = h * (0.4299D+0)


    end subroutine liq_enth




    subroutine liq_cp ( psat, p, t, cp )


! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the liquid specific heat as a function
!   of temperature and pressure
!
!   This subroutine is valid for: 91.79 deg C < T < 373.253 deg C &
!                                  saturation < P < 21.5 MPa
!
!
!   psat,  input, real ( kind = 8 ), the saturation pressure
!                                    units: psia
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!    rho, output, real ( kind = 8 ), the calculated specific heat
!                                    units: BTU/(lbm-degF)
!
        implicit none

        real ( kind = 8 ), intent ( in ) :: psat
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: cp


!       Variables used locally in this subroutine
        real ( kind = 8 ) :: psat_si
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: cp_si
        real ( kind = 8 ) :: cpf


!       Convert pressures from psia to MPa
        psat_si = psat * 0.00689476D+0
        p_si = p * 0.00689476D+0

!       Convert temperature from deg F to deg C
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


!       Correlations for cpf
        if ( psat_si < 0.030D+0 .and. psat_si > 21.50D+0 ) then

            write ( *, * ) ' THE SATURATION PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( psat_si > 0.030D+0 .and. psat_si <= 0.671D+0 ) then

            cpf = 0.247763D+0 * psat_si**(0.5704026D+0) + 4.150D+0


        else if ( psat_si > 0.671D+0 .and. psat_si <= 2.606D+0 ) then

            cpf = 0.1795305D+0 * psat_si**(0.8957323D+0) + 4.223D+0


        else if ( psat_si > 2.606D+0 .and. psat_si <= 6.489D+0 ) then

            cpf = 0.09359843D+0 * psat_si**(1.239114D+0) + 4.340D+0


        else if ( psat_si > 6.489D+0 .and. psat_si <= 11.009D+0 ) then

            cpf = 0.01068888D+0 * psat_si**(2.11376D+0) + 4.740D+0


        else if ( psat_si > 11.009D+0 .and. psat_si <= 14.946D+0 ) then

            cpf = 1.333058E-04 * psat_si**(3.707294D+0) + 5.480D+0


        else if ( psat_si > 14.946D+0 .and. psat_si <= 18.079D+0 ) then

            cpf = 6.635658E-03 * ( psat_si - 10.0D+0 )**3.223323D+0 + 7.350D+0


        else if ( psat_si > 18.079D+0 .and. psat_si <= 21.50D+0 ) then

            cpf = 4.6844786E-06 * exp( 0.7396875D+0 * psat_si) + 10.020D+0


        end if



!       Calculate specific heat per equation in Section 4.4
        cp_si = cpf + ( p_si - psat_si )*( 0.0018D+0 - 76.D+0/( 364.0D+0 - t_si)**1.8D+0)



!       Convert from kJ/(kg-k) to BTU/(lbm-degF)
        cp = cp_si * 0.2388D+0


    end subroutine liq_cp



    subroutine vap_vol ( p, tsat, t, v, rho )
! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the vapor specific volume as a function
!   of temperature and pressure
!
!   This subroutine is valid for: saturation deg C < T < 450.0 deg C &
!                                  0.085 MPa < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!   tsat,  input, real ( kind = 8 ), the saturation temperature
!                                    units: deg F
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!
!      v, output, real ( kind = 8 ), the calculated specific volume
!                                    units: ft**3/kg
!
!    rho, output, real ( kind = 8 ), the calculated density
!                                    units: kg/m**3
!
!


        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: tsat
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: v
        real ( kind = 8 ), intent ( inout ) :: rho


!       Variables used locally in subroutine
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: v_si
        real ( kind = 8 ) :: vg


!       Convert pressure from psia to MPa
        p_si = p * 0.00689476D+0

!       Convert temperatures from deg F to deg C
        tsat_si = (5.0D+0/9.0D+0)*(tsat - 32.0D+0)
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)


        if ( p_si < 0.085D+0 .and. p_si > 21.50D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si > 0.085D+0 .and. p_si < 1.112D+0 ) then

            vg = ( 5.126076D+0 * p_si**(0.9475862D+0) + 0.012D+0 )**-1.0D+0


        else if ( p_si >= 1.112D+0 .and. p_si < 3.932D+0 ) then

            vg = ( 4.630832D+0 * p_si**(1.038819D+0) + 0.520D+0 )**-1.0D+0


        else if ( p_si >= 3.932D+0 .and. p_si < 8.996D+0 ) then

            vg = ( 2.868721D+0 * p_si**(1.252148D+0) + 3.800D+0 )**-1.0D+0


        else if ( p_si >= 8.996D+0 .and. p_si < 14.628D+0 ) then

            vg = ( 0.5497653D+0 * p_si**(1.831182D+0) + 18.111 )**-1.0D+0


        else if ( p_si >= 14.628D+0 .and. p_si <= 18.210D+0 ) then

            vg = ( 8.5781572E-03 * p_si**(3.176484D+0) + 50.0D+0 )**-1.0D+0


        else if ( p_si > 18.210D+0 .and. p_si <= 20.253D+0 ) then

            vg = ( 3.5587113E-06 * p_si**(5.660939D+0) + 88.0D+0 )**-1.0D+0


        else if ( p_si > 20.253D+0 .and. p_si <= 21.50D+0 ) then

            vg = ( 3.558734E-16 * p_si**(13.03774D+0) + 138.0D+0 )**-1.0D+0

        end if


!       Formula to calculate specific volume per Section 5.1
        v_si = vg + (t_si - tsat_si) * ( 0.000466D+0/p_si - &
        ((0.12D+0/(t_si+100.0D+0) - 0.000106) * p_si**(0.1D+0))/ &
        sqrt ( 1.96E-08 * (t_si + 8.0D+0)**4.0D+0 - p_si**2.0D+0 ))


!       Convert from m**3/kg to ft**3/lbm
        v = v_si * 16.0185D+0

!       Calculate vapor density in lbm/ft**3
        rho = 1.0D+0/v



    end subroutine vap_vol


    subroutine vap_enth ( p, tsat, t, h )
! -------------------------------------------------------------------------- 80
!
!   This subroutine calculates the vapor enthalpy as a function
!   of temperature and pressure
!
!   This subroutine is valid for: saturation deg C < T < 450.0 deg C &
!                                  0.085 MPa < P < 21.5 MPa
!
!
!      p,  input, real ( kind = 8 ), the pressure
!                                    units: psia
!
!   tsat,  input, real ( kind = 8 ), the saturation temperature
!                                    units: deg F
!
!
!      t,  input, real ( kind = 8 ), the temperature
!                                    units: deg F
!
!
!      h, output, real ( kind = 8 ), the calculated enthalpy
!                                    units: BTU/lbm
!
!

        implicit none

        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: tsat
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( inout ) :: h


!       Variables used locally in subroutine
        real ( kind = 8 ) :: p_si
        real ( kind = 8 ) :: tsat_si
        real ( kind = 8 ) :: t_si
        real ( kind = 8 ) :: h_si
        real ( kind = 8 ) :: hg


!       Convert pressure from psia to MPa
        p_si = p * 0.00689476D+0

!       Convert temperatures from deg F to deg C
        tsat_si = (5.0D+0/9.0D+0)*(tsat - 32.0D+0)
        t_si = (5.0D+0/9.0D+0)*(t - 32.0D+0)



!       Correlations for hg from Section 5.2
        if ( p_si < 0.075D+0 .and. p_si > 21.55D+0 ) then

            write ( *, * ) ' THE PRESSURE USED IS OUTSIDE THE DOMAIN '
            write ( *, * ) ' THE SUBROUTINE IS EXITING '


        else if ( p_si > 0.075D+0 .and. p_si <= 0.348D+0 ) then

            hg = -4.0381938E-06 * ( 3.0D+0 - p_si )**15.72364D+0 &
            + 2750.0D+0


        else if ( p_si > 0.348D+0 .and. p_si <= 1.248D+0 ) then

            hg = -0.5767304 * exp( -1.66153D+0 * ( p_si - 3.2D+0 )) &
            + 2800.0D+0


        else if ( p_si > 1.248D+0 .and. p_si < 2.955D+0 ) then

            hg = -7.835986D+0 * ( p_si - 3.001D+0 )**2.0D+0 &
            - 2.934312D+0 * ( p_si - 3.001D+0 ) + 2803.71D+0


        else if ( p_si >= 2.955D+0 .and. p_si <= 6.522D+0 ) then

            hg = -1.347244D+0 * ( p_si - 2.999D+0 )**2.0D+0 &
            - 2.326913D+0 * ( p_si - 2.999D+0 ) + 2803.35D+0


        else if ( p_si > 6.522D+0 .and. p_si < 16.497D+0 ) then

            hg = -0.9219176D+0 * ( p_si - 9.00D+0 )**2.0D+0 &
            - 16.388835D+0 * ( p_si - 9.00D+0) + 2742.03D+0


        else if ( p_si >= 16.497D+0 .and. p_si < 20.193D+0 ) then

            hg = -3.532177D+0 * ( p_si - 8.00D+0 )**2.0D+0 &
            + 29.81305D+0 * ( p_si - 8.00D+0 ) + 2565.0D+0


        else if ( p_si >= 20.193D+0 .and. p_si <= 21.550D+0 ) then

            hg = -22.92521D+0 * ( p_si - 18.0D+0 )**2.0D+0 &
            + 44.23671D+0 * ( p_si - 18.0D+0 ) + 2415.01D+0


        end if


!       Formula to calculate enthalpy from Section 5.2
        h_si = hg + ( t_si - tsat_si ) * (( 0.28D+0 * exp( -0.008D+0 * ( t_si - 162.0D+0 ))) &
        + ( - 100.0D+0/t_si - 2.225D+0) &
        + ( 4.5D+0 * p_si / sqrt( 7.4529E-06 * t_si**3 - p_si**2 )))



!       Convert enthalpy from kJ/kg to BTU/lb
        h = h_si * 0.4299D+0


       end subroutine vap_enth



end module mod_thermo



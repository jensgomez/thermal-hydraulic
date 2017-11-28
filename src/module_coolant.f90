module mod_coolant

    use mod_thermo

    implicit none

    contains


    subroutine coolant ( mdot, qp, tin, pin, nrods, l, lstep, nassemb, nshape )

! ************************************************************************** 80
!
! Discussion:
!   This subroutine contains all the input parameters needed to calculate
!   the thermodynamic and hydraulic properties of the coolant
!
! Arguments:
!
!    mdot, input,     real ( kind = 8 ), this is the mass flow rate of the
!                                        coolant entering the core
!                                        units: lb/hr
!
!
!      qp, input,     real ( kind = 8 ), this is the linear heat generation rate
!                                        for the channel
!                                        units: BTU/(hr-ft)
!
!
!     tin, input,     real ( kind = 8 ), this is the inlet temperature of the fluid
!                                        units: deg F
!
!
!     pin, input,     real ( kind = 8 ), this is the inlet pressure of the fluid
!                                        units: psi
!
!
!   nrods, input,  integer ( kind = 4 ), this is the number of fuel rods in the core
!                                        units: dimensionless
!
!
!       l, input,     real ( kind = 8 ), this is the active fuel length
!                                        units: inches
!
!
!   lstep, input,  integer ( kind = 4 ), this is the number of axial elevations in the channel
!                                        units: dimensionless
!
!
! nassemb, input,  integer ( kind = 4 ), the number of fuel assemblies
!                                        units: dimensionless
!
!
!  nshape, input,  integer ( kind = 4 ), the fuel assembly shape
!                                        units: dimensionless

        implicit none


!       Input arguments
        real    ( kind = 8 ), intent ( in ) :: mdot
        real    ( kind = 8 ), intent ( in ) :: qp(lstep)
        real    ( kind = 8 ), intent ( in ) :: tin
        real    ( kind = 8 ), intent ( in ) :: pin
        integer ( kind = 4 ), intent ( in ) :: nrods
        real    ( kind = 8 ), intent ( in ) :: l
        integer ( kind = 4 ), intent ( in ) :: lstep
        integer ( kind = 4 ), intent ( in ) :: nassemb
        integer ( kind = 4 ), intent ( in ) :: nshape


!       Local variables
        real    ( kind = 8 ) :: dzft
        real    ( kind = 8 ) :: qdot
        real    ( kind = 8 ) :: rho
        real    ( kind = 8 ) :: psat
        real    ( kind = 8 ) :: h(lstep)
        real    ( kind = 8 ) :: hin
        real    ( kind = 8 ) :: dh
        real    ( kind = 8 ) :: x
        integer ( kind = 4 ) :: i
        real    ( kind = 8 ) :: s
        real    ( kind = 8 ) :: area
        integer ( kind = 4 ) :: opt
        real    ( kind = 8 ) :: void
        real    ( kind = 8 ) :: mdotchannel

        h(1:lstep) = 0.0D+0



!       Calculating coolant density at entrance
        call p_sat ( tin, psat )
        call liq_dens ( psat, pin, tin, rho )

!       Calculating volumetric flow rate
        qdot = mdot / rho

!       Calculating coolant enthalpy at entrance
        call liq_enth ( psat, pin, tin, hin )

!       Set the enthalpy at the entrance equal to the first value of the enthalpy array
        h(1) = hin

!       Calculate the length differential
        dzft = (l/12.0D+0)/lstep


        !       Open file stream
        open ( unit = 10, file = 'coolant_information.out', status = 'replace')
        write ( 10, '(a)' ) ' COOLANT SUBROUTINE CALLED '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' )         '*** INPUTS LISTED BELOW *** '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a, g14.6, a)' )  ' ACTIVE FUEL LENGTH                      = ', l ,      ' INCHES'
        write ( 10, '(a, 2x, i4, a)' ) ' THE NUMBER OF FUEL ASSEMBLIES           = ', nassemb, ' '
        write ( 10, '(a, 2x, i2, a2, i2, a)' ) ' THE FUEL SHAPE FOR ONE ASSEMBLY         = ', int(sqrt(dble(nshape))), ' X ', &
         int(sqrt(dble(nshape))), ' '
        write ( 10, '(a, g14.6, a)' ) ' THE MASS FLOW RATE                      = ', mdot  ,   ' LB/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE VOLUMETR FLOW RATE                  = ', qdot  ,   ' FT**3/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET TEMPERATURE                   = ', tin   ,   ' DEG F '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET PRESSURE                      = ', pin   ,   ' PSI '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT ENTHALY AT THE INLET        = ', hin  ,   ' BTU/LBM '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT DENSITY AT THE INLET        = ', rho   ,   ' LBM/FT**3 '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' ) ' '


        write ( 10, '(a)') ' THE FOLLOWING INFORMATION IS CALCULATED AS THE FLUID TRAVELS UP THE FUEL CHANNEL '
        write ( 10, '(a)' ) ' '


!       Calculate various coolant properties in this loop
        do i = 2, lstep

            call enthalpy ( mdot, qp(i), nassemb, nshape, dzft, dh )
            h(i) = h(i-1) + dh

            call quality ( tin, pin, h(i), x )


!           Following area value from Todreas and Kazimi Appendix K
!           Assembly flow area for BWR-5 = 9.718E-03 m**2 =
            area = 0.1046D+0

!           Calculate the mass flow rate into the channel in lb/hr
            mdotchannel = mdot/(nassemb)
            write(*,*) mdotchannel

!           Opt = 1 for CISE slip correlation
            opt = 1
            call slip ( x, tin, pin, mdotchannel , area, opt, s, void )


            write ( 10, '(a, g14.6, a)' ) ' THE AXIAL ELEVATION             = ', dzft*i , ' FT '
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED ENTHALPY         = ', h(i)   , ' BTU/LBM'
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED QUALITY          = ', x*100.0D+0 , ' %'
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED SLIP FACTOR      = ', s
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED VOID FRACTION    = ', void
            write ( 10, '(a)' ) ' '



        end do




    end subroutine coolant













    subroutine enthalpy ( mdot, qp, nassemb, nshape, dzft, dh )

! ************************************************************************** 80
!
! Discussion:
!   This subroutine is created to calculate the enthalpy rise of the coolant
!   within the fuel channel
!
! Arguments:
!
!    mdot, input,     real ( kind = 8 ), this is the mass flow rate of the
!                                        coolant entering the core
!                                        units: lb/hr
!
!
!      qp, input,     real ( kind = 8 ), this is the linear heat generation rate
!                                        for the channel
!                                        units: BTU/(hr-ft)
!
!
! nassemb, input,  integer ( kind = 4 ), the number of fuel assemblies
!                                        units: dimensionless
!
!
!  nshape, input,  integer ( kind = 4 ), the fuel assembly shape
!                                        units: dimensionless
!
!
!    dzft, input,     real ( kind = 4 ), the axial length step
!                                        units: ft
!
!
!      dh, output,    real ( kind = 8 ), the calculated enthalpy rise
!                                        units: BTU/(lbm)
!

        implicit none


!       Input arguments
        real    ( kind = 8 ), intent ( in ) :: mdot
        real    ( kind = 8 ), intent ( in ) :: qp
        integer ( kind = 4 ), intent ( in ) :: nassemb
        integer ( kind = 4 ), intent ( in ) :: nshape
        real    ( kind = 8 ), intent ( in ) :: dzft
        real    ( kind = 8 ), intent ( inout ) :: dh

!       Local variables
        real    ( kind = 8 ) :: mchannel
        real    ( kind = 8 ) :: q

!       Calculate the channel mass flow rate
        mchannel = mdot/nassemb

!       Calculate the power from the linear heat generation rate
        q = qp*dzft*nshape

!       Calculate the enthalpy rise
        dh = q/mchannel


    end subroutine enthalpy




    subroutine quality ( t, p, h, x )

! ************************************************************************** 80
!
! Discussion:
!   This subroutine is created to calculate the quality of the coolant
!   within the fuel channel
!
! Arguments:
!
!       t, input,     real ( kind = 8 ), this is the temperature of fluid
!                                        units: deg F
!
!
!       p, input,     real ( kind = 8 ), this is the pressure of the coolant
!                                        units: psia
!
!
!       h, input,     real ( kind = 8 ), this is the enthalpy of the coolant
!                                        units: BTU/lbm
!
!
!       x, input,     real ( kind = 8 ), this is the calculate quality
!                                        units: dimensionless
!


        implicit none

!       Input arguments
        real ( kind = 8 ), intent ( in ) :: t
        real ( kind = 8 ), intent ( in ) :: p
        real ( kind = 8 ), intent ( in ) :: h
        real ( kind = 8 ), intent ( inout ) :: x

!       Local variables
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: tsat
        real ( kind = 8 ) :: h_l
        real ( kind = 8 ) :: h_v
        character ( len = 25 ) :: state

!       Call saturated liquid and saturated vapor subroutines
        call p_sat ( t, psat )
        call liq_enth ( psat, p, t, h_l )

        call t_sat ( p, tsat )
        call vap_enth ( p, tsat, t, h_v )


!       Calculate the quality of the coolant
!       The logicals are used so that NaN errors are avoided
        if ( isnan(h_v) .eqv. .true. ) then

            x = -1.0D+0

        else if ( isnan(h_l) .eqv. .true. ) then

            x = 1.5D+0

        else if ( isnan(h_l) .and. isnan(h_v) .eqv. .false. ) then

            x = ( h - h_l ) / ( h_v - h_l )

        end if



!       Determine state of coolant
        if ( x < 0.0D+0 ) then

            x = 0.0D+0
            state = ' COMPRESSED LIQUID '


        else if ( x == 0.0D+0 ) then

            state = ' SATURATED LIQUID '


        else if ( x > 0.0D+0 .and. x < 1.0D+0 ) then

            state = ' LIQUID - VAPOR MIX '


        else if ( x == 1.0D+0 ) then

            state = ' SATURATED VAPOR '



        else if ( x > 1.0D+0 ) then

            x = 1.0D+0
            state = ' SUPERHEATED VAPOR '

        end if


    end subroutine quality



    subroutine slip ( x, t, p, mdot, area, opt, s, void )
! ************************************************************************** 80
!
! Discussion:
!   This subroutine is created to calculate the slip factor for use in the
!   calculation of the void fraction
!
! Arguments:
!
!       x, input,     real ( kind = 8 ), this is the quality of the coolant
!                                        units: dimensionless
!
!
!       t, input,     real ( kind = 8 ), this is the temperature of the coolant
!                                        units: deg F
!
!
!       p, input,     real ( kind = 8 ), this is the pressure of the coolant
!                                        units: psia
!
!
!    mdot, input,     real ( kind = 8 ), the mass flow rate
!                                        units: lb/hr
!
!
!    area, input,     real ( kind = 8 ), the flow area
!                                        units: ft**2
!
!
!     opt, input,   integer ( kind = 4), user option to select which slip correlation
!                                        to use
!
!
!       s, output,     real ( kind = 8 ), this is the calculated slip factor
!                                         units: dimensionless
!
!    void, output,     real ( kind = 8 ), the void fraction
!                                         units: dimensionless
!
!

        implicit none


!       Input arguments
        real    ( kind = 8 ) :: x
        real    ( kind = 8 ) :: t
        real    ( kind = 8 ) :: p
        real    ( kind = 8 ) :: mdot
        real    ( kind = 8 ) :: area
        integer ( kind = 4 ) :: opt
        real    ( kind = 8 ) :: s
        real    ( kind = 8 ) :: void


!       Local variables
        real ( kind = 8 ) :: plocal
        real ( kind = 8 ) :: rho_l
        real ( kind = 8 ) :: rho_v
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: tsat
        real ( kind = 8 ) :: d
        real ( kind = 8 ) :: v
        real ( kind = 8 ) :: sigma
        real ( kind = 8 ) :: rey
        real ( kind = 8 ) :: web
        real ( kind = 8 ) :: g
        real ( kind = 8 ) :: mu
        real ( kind = 8 ) :: e1
        real ( kind = 8 ) :: e2
        real ( kind = 8 ) :: beta
        real ( kind = 8 ) :: y
        real ( kind = 8 ), parameter :: pi = 4.0D+0*atan(1.0D+0)

!       Calculate the liquid and vapor density
        call p_sat ( t, psat )
        call liq_dens ( psat, p, t, rho_l )

        call t_sat ( p, tsat )
        call vap_vol ( p, tsat, t, v, rho_v )

!       Calculate the liquid surface tension
        call surtension ( t , sigma )

!       Calculate flow diameter
        d = sqrt(4.0D+0*area/pi)

!       Calculate the viscosity
        mu = viscosity ( t, p )


!       If the opt parameter = 1, implement CISE correlation
        if ( opt == 1 ) then

!           Calculate mass flux (lb/sec**ft**2)
            g = mdot/(area*3600.0D+0)

!           Calculate the Reynolds number
            rey = g*d/mu

!           Calculate the Webber number
            web = d*g**2/(rho_l*sigma)

!           CISE correlation parameters below
            e1 = 1.578D+0*rey**(-0.19D+0)*(rho_l/rho_v)**(0.22D+0)

            e2 = 0.0273*web*rey**(-0.51D+0)* &
            (rho_l/rho_v)**(-0.08D+0)

            beta = rho_l*x/(rho_l*x + rho_v*(1.0D+0 - x))

            y = beta/(1.0D+0 - beta)

            s = 1.0D+0 + e1*sqrt(y/(1.0D+0+y*e2) - y*e2)

!       Implement Chisholm correlation
        else if ( opt == 2 ) then

            s = sqrt ( 1.0D+0 - x*( 1.0D+0 - rho_l/rho_v ))

        end if

        void = 1.0D+0/(1.0D+0 + (s*((1.0D+0-x)/x)*(rho_v/rho_l)))


    end subroutine slip




end module mod_coolant


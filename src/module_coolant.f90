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


            write ( 10, '(a, g14.6, a)' ) ' THE AXIAL ELEVATION             = ', dzft*i , ' FT '
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED ENTHALPY         = ', h(i)   , ' BTU/LBM'
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED QUALITY          = ', x*100.0D+0 , ' %'
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




    subroutine quality ( t, p, h , x )

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



end module mod_coolant


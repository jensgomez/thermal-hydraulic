module mod_fluid

    use mod_thermo

    implicit none

    contains


    subroutine enth_rise ( mdot, qp, tin, pin, nrods, l, lstep, nassemb, nshape, h )

! ************************************************************************** 80
!
! Discussion:
!   This subroutine is created to calculate the enthalpy rise within the
!   fuel channel
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
!
!
!       h, output,    real ( kind = 8 ), the calculated enthalpy rise
!                                        units: BTU/(lbm)
!

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
        real    ( kind = 8 ), intent ( inout ) :: h(lstep)

!       Local variables
        real    ( kind = 8 ) :: psat
        real    ( kind = 8 ) :: rho
        real    ( kind = 8 ) :: qdot
        real    ( kind = 8 ) :: mchannel
        real    ( kind = 8 ) :: dzft
        real    ( kind = 8 ) :: q(lstep)
        real    ( kind = 8 ) :: cp(lstep)
        real    ( kind = 8 ) :: t(lstep)

        integer ( kind = 4 ) :: i

!       Assign the enthalpy to zeros
        h(1:lstep) = 0.0D+0

!       Calling the saturated pressure subroutine
        call p_sat ( tin, psat )

!       Calculate the entrance coolant density
        call liq_dens ( psat, pin, tin, rho )

!       Calculate the inlet temperature enthalpy
        call liq_enth ( psat, pin, tin, h(1) )

!       Calculate the inlet specific heat
        call liq_cp ( psat, pin, tin, cp(1) )


!       The entrance volumetric flow rate
        qdot = mdot / rho

!       Calculate the channel flow
        mchannel = mdot/nassemb

!       The length differential
        dzft = (l/12.0D+0)/lstep


        t(1) = tin

        do i = 2, lstep

            call p_sat ( t(i-1), psat)
            call liq_cp ( psat, pin, t(i), cp(i) )

            q(i) = qp(i)*dzft*nshape
            h(i) = h(i-1) + q(i)/mchannel

            call liq_cp ( psat, pin, t(i-1), cp(1) )

            t(i) = q(i)/(mchannel*cp(i)) + t(i-1)


        end do





!       Open file stream
        open ( unit = 10, file = "fluid_information.out" )
        write ( 10, '(a)' ) ' FLUID  SUBROUTINE CALLED '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' )         '*** INPUTS LISTED BELOW *** '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a, g14.6, a)' ) ' ACTIVE FUEL LENGTH                      = ', l ,      ' INCHES'
        write ( 10, '(a, 2x, i4, a)' ) ' THE NUMBER OF FUEL ASSEMBLIES           = ', nassemb, ' '
        write ( 10, '(a, 2x, i2, a2, i2, a)' ) ' THE FUEL SHAPE FOR ONE ASSEMBLY         = ', int(sqrt(dble(nshape))), ' X ', &
         int(sqrt(dble(nshape))), ' '
        write ( 10, '(a, g14.6, a)' ) ' THE MASS FLOW RATE                      = ', mdot  ,   ' LB/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE VOLUMETR FLOW RATE                  = ', qdot  ,   ' FT**3/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET TEMPERATURE                   = ', tin   ,   ' DEG F '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET PRESSURE                      = ', pin   ,   ' PSI '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT ENTHALY AT THE INLET        = ', h(1)  ,   ' BTU/LBM '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT DENSITY AT THE INLET        = ', rho   ,   ' LBM/FT**3 '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT SPECIFIC HEAT AT THE INLET  = ', cp(1) ,   ' BTU/(LBM-DEG F) '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' ) ' '


        write ( 10, '(a)') ' THE FOLLOWING INFORMATION IS CALCULATED AS THE FLUID TRAVELS UP THE FUEL CHANNEL '
        write ( 10, '(a)' ) ' '
        do i = 1, lstep

            write ( 10, '(a, g14.6, a)' ) ' THE AXIAL ELEVATION             = ', dzft*i , ' FT '
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED ENTHALPY         = ', h(i)   , ' BTU/LBM'
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATED SPECIFIC HEAT    = ', cp(i)  , ' BTU/LBM'
            write ( 10, '(a, g14.6, a)' ) ' THE CALCULATE TEMPERATURE       = ', t(i) ,   ' DEG F '
            write ( 10, '(a)' ) ' '
            write ( 10, '(a)' ) ' '


        end do





        end subroutine enth_rise


end module mod_fluid


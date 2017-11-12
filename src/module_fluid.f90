module mod_fluid

    use mod_thermo

    implicit none

    contains


    subroutine enth_rise ( mdot, qp, tin, pin, nrods, l, lstep, h )

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
!       h, output,    real ( kind = 8 ), the calculated enthalpy rise
!                                        units: BTU/(lbm)
!

        implicit none


!       Input arguments
        real ( kind = 8 ), intent ( in ) :: mdot
        real ( kind = 8 ), intent ( in ) :: qp(lstep)
        real ( kind = 8 ), intent ( in ) :: tin
        real ( kind = 8 ), intent ( in ) :: pin
        integer ( kind = 4 ), intent ( in ) :: nrods
        real ( kind = 8 ), intent ( in ) :: l
        integer ( kind = 4 ), intent ( in ) :: lstep
        real ( kind = 8 ), intent ( inout ) :: h(lstep)

!       Local variables
        real ( kind = 8 ) :: psat
        real ( kind = 8 ) :: rho
        real ( kind = 8 ) :: qdot

!       Assign the enthalpy to zeros
        h(1:lstep) = 0.0D+0

!       Calling the saturated pressure subroutine
        call p_sat ( tin, psat )

!       Calculate the entrance coolant density
        call liq_dens ( psat, pin, tin, rho )

!       Calculate the inlet temperature enthalpy
        call liq_enth ( psat, pin, tin, h(1) )

!       The entrance volumetric flow rate
        qdot = mdot / rho

!       Open file stream
        open ( unit = 10, file = "fluid_information.out" )
        write ( 10, '(a)' ) ' FLUID  SUBROUTINE CALLED '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' )         '*** INPUTS LISTED BELOW *** '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a, g14.6, a)' ) ' ACTIVE FUEL LENGTH                = ', l ,      ' INCHES'
        write ( 10, '(a, g14.6, a)' ) ' THE MASS FLOW RATE                = ', mdot ,   ' LB/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE VOLUMETR FLOW RATE            = ', qdot ,   ' FT**3/HR '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET TEMPERATURE             = ', tin  ,   ' DEG F '
        write ( 10, '(a, g14.6, a)' ) ' THE INLET PRESSURE                = ', pin  ,   ' PSI '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT ENTHALY AT THE INLET  = ', h(1) ,   ' BTU/LBM '
        write ( 10, '(a, g14.6, a)' ) ' THE COOLANT DENSITY AT THE INLET  = ', rho  ,   ' LBM/FT**3 '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' ) ' '










        end subroutine enth_rise


































end module mod_fluid

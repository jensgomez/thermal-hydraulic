module mod_fuel

    implicit none

    contains

    subroutine fuelshape ( option, l, pow, peak, a, lstep, nrods, qp )

! ************************************************************************** 80
!
! Discussion:
!   This subroutine is created to calculate the linear heat generation rate
!   for a fuel channel with user inputs
!
!
! Arguments:
!
!   option,  input, integer ( kind = 4 ), this determines what fuel  shape
!                                         will be used within the subroutine.
!                                         units: N/A
!
!
!        l,  input,    real ( kind = 8 ), this is the length of the fuel
!                                         unit: inches
!
!
!      pow,  input,    real ( kind = 8 ), this is the rate power of the entire active fuel
!                                         units: MWt
!
!
!     peak,  input,    real ( kind = 8 ), this is the peaking factor constant for the fuel
!                                         units: dimensionless
!
!
!       a,   input,    real ( kind = 8 ), this is a constant that can move the fuel shape axially
!                                         units: dimensionless
!
!
!   lstep,   input, integer ( kind = 4 ), this is the number of axially steps that the user
!                                         wants to model
!                                         units: dimensionless
!
!
!   nrods,   input, integer ( kind = 4 ), this is the number of fuel rods in the core
!                                         units: dimensionless
!
!
!     qp,   output,    real ( kind  = 8), this is the linear heat generation rate as a function of
!                                         axial elevation
!                                         units: BTU/(hr-ft)
!
!

        implicit none

!       Subroutine arguments
        integer ( kind = 4 ), intent ( in ) :: option
        real (kind = 8 ), intent ( in ) :: l
        real (kind = 8 ), intent ( in ) :: pow
        real (kind = 8 ), intent ( in ) :: peak
        real (kind = 8 ), intent ( in ) :: a
        integer ( kind = 4 ), intent ( in ) :: lstep
        integer ( kind = 4 ), intent ( in ) :: nrods
        real (kind = 8 ), intent ( inout ) :: qp(lstep)

!       Locally used variables
        integer ( kind = 4 ) :: i
        real ( kind = 8 ), parameter :: pi = acos( -1.0 )
        real ( kind = 8 ), parameter :: mwt2btuhr = 3.412E+06
        real ( kind = 8 ) :: dz
        real ( kind = 8 ) :: z(lstep)
        real ( kind = 8 ) :: q2qp
        real ( kind = 8 ) :: qpmean



!       Open file stream
        open(unit = 10, file = "Fuel_information.out")
        write ( 10, '(a)' ) 'FUEL SUBROUTINE CALLED '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a)' )         '*** INPUTS LISTED BELOW *** '
        write ( 10, '(a)' ) ' '
        write ( 10, '(a, g14.6, a)' ) ' ACTIVE FUEL LENGTH        = ', l,     ' INCHES'
        write ( 10, '(a, g14.6, a)' ) ' THE FUEL RATED POWER      = ', pow,   ' MWT '
        write ( 10, '(a, g14.6, a)' ) ' THE FUEL PEAKING FACTOR   = ', peak,  ' '
        write ( 10, '(a)' ) ' '


        dz = l/lstep

!       Create axial elevation array
        do i = 1,lstep
            z(i) = dz*i
        end do


!       Q to q' factor
        q2qp = pow/(nrods*l/12.0D+0)


        if ( option == 1 ) then

            write ( 10, '(a)' ) ' '
            write ( 10, '(a)' ) ' SINE FUEL SHAPE SELECTED '
            write ( 10, '(a)' ) ' '
            write ( 10, '(a)' ) ' '

            do i = 1, lstep

                qp(i) = exp(-a*z(i)/l) * sin(pi*z(i)/l) * mwt2btuhr * peak * q2qp

            end do


        else if ( option == 2 ) then

            write ( 10, '(a)' ) ' '
            write ( 10, '(a)' ) ' FLAT FUEL SHAPE SELECTED '
            write ( 10, '(a)' ) ' '
            write ( 10, '(a)' ) ' '

            do i = 1, lstep

                qp(i) = mwt2btuhr * peak * q2qp

            end do

        end if


        do i = 1, lstep

            write ( 10 , '(a, i5)')        ' THE LSTEP VALUE                 = ', i
            write ( 10 , '(a, g14.6, a)')  ' THE DZ VALUE                    = ', dz*i  , '   INCHES '
            write ( 10 , '(a, g14.6, a)')  ' THE LINEAR HEAT GENERATION RATE = ', qp(i) , '   BTU/(HR-FT) '
            write ( 10 , '(a, g14.6, a)')  ' THE LINEAR HEAT GENERATION RATE = ', qp(i)*9.615E-4 , '   KW/M '
            write ( 10 , '(a)' ) ' '
            write ( 10 , '(a)' ) ' '

        end do


        qpmean = sum(qp)/lstep

        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a, g14.6, a)')  ' THE AVERAGE LINEAR HEAT GENERATION RATE  = ', qpmean , '   BTU/(HR-FT)'
        write ( 10 , '(a, g14.6, a)')  ' THE AVERAGE LINEAR HEAT GENERATION RATE  = ', qpmean*9.615E-04 , '   KW/M '
        write ( 10 , '(a, g14.6, a)')  ' THE AVERAGE LINEAR HEAT GENERATION RATE  = ', qpmean*2.931E-04 , '   KW/FT '
        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a, g14.6, a)')  ' THE MAXIUMUM LINEAR HEAT GENERATION RATE = ', maxval(qp) , '   BTU/(HR-FT)'
        write ( 10 , '(a, g14.6, a)')  ' THE MAXIUMUM LINEAR HEAT GENERATION RATE = ', maxval(qp)*9.615E-04 , '   KW/M '
        write ( 10 , '(a, g14.6, a)')  ' THE MAXIUMUM LINEAR HEAT GENERATION RATE = ', maxval(qp)*2.931E-04 , '   KW/FT '
        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a)' ) ' '
        write ( 10 , '(a, f6.2, a, i3, a)') ' THE MAXIUMUM LINEAR HEAT GENERATION &
                                                 &RATE OCCURS AT = ',  maxloc(qp)*dz, ' INCHES ( LSTEP = ', maxloc(qp), ' )'

    end subroutine fuelshape


end module mod_fuel



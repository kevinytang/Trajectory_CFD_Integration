!> Gravity models
!> Author: Kevin Tang
module mod_gravity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: EarthPointMass

contains

    ! Model 1: Earth point-mass gravity
    function EarthPointMass(alt) result(g)
        implicit none

        ! Input altitude
        real(dp), intent(in) :: alt ! [m]

        ! Constants
        real(dp), parameter :: mu = 3.986004418e14_dp ! gravitational parameter of Earth [m^3/s^2]
        real(dp), parameter :: REar = 6.378137e6_dp ! radius of earth [m]

        ! Output gravitational acceleration
        real(dp) :: g ! [m/s^2]

        ! Point mass model
        g = mu / (REar + alt)**2

    end function EarthPointMass

end module mod_gravity

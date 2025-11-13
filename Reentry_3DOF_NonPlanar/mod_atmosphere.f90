!> Atmospheric models for density
!> Author: Kevin Tang
module mod_atmosphere
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: EarthAtmNASA

contains

    ! Model 1: Earth Atmosphere by NASA [https://www.grc.nasa.gov/www/k-12/airplane/atmosmet.html]
    subroutine EarthAtmNASA(alt, T, P, rho, a)
        implicit none

        ! Input altitude
        real(dp), intent(in) :: alt ! [m]

        ! Outputs
        real(dp), intent(out) :: T ! temperature [C]
        real(dp), intent(out) :: P ! pressure [KPa]
        real(dp), intent(out) :: rho ! density [kg/m^3]
        real(dp), intent(out) :: a ! speed of sound

        ! Constants
        real(dp), parameter :: T_ref = 273.15_dp ! reference temperature [K]
        real(dp), parameter :: S = 110.4_dp ! Sutherland temperature for air [K]
        real(dp), parameter :: gamma = 1.4_dp ! Heat capacity ratio assume CPG for now
        real(dp), parameter :: R = 287_dp ! Specific gas constant of air [J/kg/K]

        ! Upper stratosphere
        if (alt > 25000.0_dp) then
            T = -131.21_dp + 0.00299_dp * alt
            P = 2.488_dp * ((T + 273.1_dp) / 216.6_dp)**(-11.338_dp)

            ! Speed of sound [...] for now, assume CPG
            a = sqrt(gamma * R * (T + T_ref))

        ! Lower stratosphere
        else if (alt >= 11000.0_dp .and. alt <= 25000.0_dp) then
            T = -56.46_dp
            P = 22.65_dp * exp(1.73_dp - 0.000157_dp * alt)

            ! Speed of sound
            a = sqrt(gamma * R * (T + T_ref))

        ! Troposphere
        else
            T = 15.04_dp - 0.00649_dp * alt
            P = 101.29_dp * ((T + 273.1_dp) / 288.08_dp)**(5.256_dp)

            ! Speed of sound
            a = sqrt(gamma * R * (T + T_ref))

        end if

        ! Calculate air density at given altitude
        rho = P / (0.2869_dp * (T + 273.1_dp)) ! [kg/m^3]

    end subroutine EarthAtmNASA

end module mod_atmosphere

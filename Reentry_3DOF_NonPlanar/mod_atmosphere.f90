!> Atmospheric models for density
!> Author: Kevin Tang
module mod_atmosphere
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: EarthAtmNASA
    public :: EarthAtmNRLMSISE

    ! Interface to the FORTRAN 77 subroutine for the NRLMSISE-00 model
    interface
        subroutine GTD7(IYD, SEC ,ALT ,GLAT ,GLONG ,STL ,F107A ,F107 ,AP ,MASS ,D ,T)

            ! Inputs
            integer, intent(in) :: IYD, MASS
            real, intent(in) :: SEC, ALT, GLAT, GLONG, STL, F107A, F107, AP

            ! Outputs
            real, intent(out) :: D(9), T(2)
        end subroutine GTD7

        subroutine METERS(METER)
        ! Convert outputs to Kg & Meters if METER true

            logical, intent(in) :: METER
        end subroutine METERS
    end interface

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

    ! NRLMSISE Atmospheric Model
    subroutine EarthAtmNRLMSISE(iyd, sec_val, alt_m, glat, glong, f107a_val, f107_val, ap_val, rho_out, &
        T_out, P_out, a_out)
        implicit none
        ! Input arguments:
        !
        !   IYD    - Year and day in YYDDD format (day of year 1–365 or 366).
        !            (Year is ignored in the current model.)
        !   SEC    - Universal Time (seconds).
        !   ALT    - Altitude [km].
        !   GLAT   - Geodetic latitude [deg].
        !   GLONG  - Geodetic longitude [deg].
        !   STL    - Local apparent solar time [hours]. See notes below.
        !   F107A  - 81-day average of F10.7 solar flux (centered on day DDD).
        !   F107   - Daily F10.7 solar flux for the previous day.
        !   AP     - Magnetic index:
        !            If SW(9) ≠ -1: scalar daily Ap value.
        !            If SW(9) = -1: array containing:
        !              (1) Daily Ap
        !              (2) 3-hour Ap index for current time
        !              (3) 3-hour Ap index for 3 hours before current time
        !              (4) 3-hour Ap index for 6 hours before current time
        !              (5) 3-hour Ap index for 9 hours before current time
        !              (6) Average of eight 3-hour Ap indices from 12–33 hours prior
        !              (7) Average of eight 3-hour Ap indices from 36–57 hours prior
        !   MASS   - Mass number selector:
        !              0  = temperature only
        !             48  = all species
        !             17  = anomalous oxygen only
        !            Otherwise, only the selected species density is computed.
        !
        ! Notes on input arguments:
        !
        !   UT, local time, and longitude are used independently in the model and
        !   are not equally important in all situations. For physically consistent
        !   calculations, these variables should satisfy:
        !
        !       STL = SEC / 3600 + GLONG / 15
        !
        !   Equation-of-time corrections for apparent local time may be applied if
        !   available, but are of minor importance.
        !
        !   The F107 and F107A values used by the model correspond to the 10.7 cm
        !   solar radio flux at the actual Earth–Sun distance, not at 1 AU.
        !   Data sources are available at:
        !     ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
        !
        !   Below 80 km altitude, the effects of F107, F107A, and AP are small and
        !   poorly established. Recommended values in this regime are:
        !     F107 = 150, F107A = 150, AP = 4.
        !
        ! =========================
        ! Output arguments:
        !
        !   D(1) - He number density [cm^-3]
        !   D(2) - O  number density [cm^-3]
        !   D(3) - N2 number density [cm^-3]
        !   D(4) - O2 number density [cm^-3]
        !   D(5) - Ar number density [cm^-3]
        !   D(6) - Total mass density [g/cm^3]
        !   D(7) - H  number density [cm^-3]
        !   D(8) - N  number density [cm^-3]
        !   D(9) - Anomalous oxygen number density [cm^-3]
        !
        !   T(1) - Exospheric temperature [K]
        !   T(2) - Temperature at altitude ALT [K]
        !
        ! =========================
        ! Notes on output arguments:
        !
        !   To obtain output in SI units (m^-3 and kg/m^3), call:
        !
        !       call meters(.true.)
        !
        ! =========================

        ! Inputs
        integer, intent(in) :: iyd
        real(dp), intent(in) :: sec_val, alt_m, glat, glong, f107a_val, f107_val, ap_val

        ! Outputs
        real(dp), intent(out) :: rho_out, T_out, P_out, a_out

        ! Variables
        integer :: mass
        real :: stl_f, sec_f, alt_f, glat_f, glong_f, f107a_f, f107_f, ap_f, D(9), T_msis(2)

        ! Constants
        real(dp), parameter :: gamma = 1.4_dp      ! Heat capacity ratio
        real(dp), parameter :: R_air = 287.0_dp   ! Specific gas constant [J/kg/K]

        ! Convert inputs
        sec_f    = real(sec_val)
        alt_f    = real(alt_m / 1000.0_dp) ! Convert to km
        glat_f   = real(glat)
        glong_f  = real(glong)
        f107a_f  = real(f107a_val)
        f107_f   = real(f107_val)
        ap_f     = real(ap_val)
        stl_f = sec_f/3600.0 + glong_f/15.0      ! Local solar time

        ! Mass = 48 means compute all species
        mass = 48

        ! Call NRLMSISE-00
        call GTD7(iyd, sec_f, alt_f, glat_f, glong_f, stl_f, f107a_f, f107_f, ap_f, mass, D, T_msis)

        ! Extract outputs
        rho_out = real(D(6), dp) * 1000.0_dp ! g/cm³, convert to kg/m³
        T_out = real(T_msis(2), dp)
        P_out = rho_out * R_air * T_out / 1000.0_dp ! Convert to KPa
        a_out = sqrt(gamma * R_air * T_out)

    end subroutine EarthAtmNRLMSISE

end module mod_atmosphere

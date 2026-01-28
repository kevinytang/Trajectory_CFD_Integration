!> Non-Planar 3-DOF Reentry Trajectory Propagator
!> Author: Kevin Tang
program ReentryNonPlanar3DOF
    !==================== Declaration of Variables ====================
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_atmosphere, only: EarthAtmNASA, EarthAtmNRLMSISE
    use mod_gravity, only: EarthPointMass
    implicit none

    ! State variables
    real(dp) :: state(6), state_new(6), t

    ! Integration control
    integer, parameter :: steps = 100
    real(dp), parameter :: t0 = 0.0_dp, tend = 5.0_dp ! Propagate 5 seconds at each call
    real(dp), parameter :: dt = (tend - t0) / steps
    integer :: i

    ! Vehicle parameters
    real(dp), parameter :: m    = 2000.0_dp ! [kg]
    real(dp), parameter :: Aref = 12.566_dp ! [m^2]
    real(dp), parameter :: Lref = 4_dp ! [m] [...] for now. Use higher fidelity model later

    ! Thrust parameters (set to zero for ballistic reentry)
    real(dp), parameter :: T_thrust = 0.0_dp
    real(dp), parameter :: epsilon  = 0.0_dp

    ! Bank angle control
    real(dp), parameter :: sigma = 0.0_dp

    ! Planet constants
    real(dp), parameter :: REar = 6.378137e6_dp
    real(dp), parameter :: omega_earth = 7.2921159e-5_dp ! Earth rotation rate [rad/s]

    ! Aerodynamics variables
    real(dp) :: Temp, P, rho, a ! temperature [C] pressure [KPa] density [kg/m^3] viscosity [kg/m/s] speed of sound [m/s]

    ! Helper constants
    real(dp), parameter :: deg2rad = 3.14159265359_dp / 180.0_dp
    real(dp), parameter :: rad2deg = 180.0_dp / 3.14159265359_dp
    real(dp) :: g, L_D, beta_param ! L/D and ballistic coef

    ! Declare vars that match the namelist
    real(dp) :: step, time, V, gamma, psi, alt, lon, lat, Cl, Cd
    real(dp) :: Mach, temperature, Alpha
    namelist /current_states/ step, time, V, gamma, psi, alt, lon, lat, Cl, Cd
    namelist /cfd_variables/ Mach, temperature, Alpha

    ! Atmosphere model selection and NRLMSISE-00 parameters
    integer :: atm_model ! 1 = NASA, 2 = NRLMSISE-00
    integer :: iyd ! Year and day as YYDDD (e.g., 25172 = day 172 of 2025)
    real(dp) :: sec_ut ! UT seconds of day
    real(dp) :: f107a ! 81-day average F10.7 solar flux
    real(dp) :: f107 ! Daily F10.7 solar flux (previous day)
    real(dp) :: ap ! Geomagnetic Ap index
    namelist /atmosphere_input/ atm_model, iyd, sec_ut, f107a, f107, ap

    ! Initial conditions (read-only, preserved for reference)
    real(dp) :: v0, gamma0, psi0, alt0, lon0, lat0
    namelist /initial_conditions/ v0, gamma0, psi0, alt0, lon0, lat0

    ! Control settings (read-only, preserved for reference)
    real(dp) :: t_step, t_end_ctrl, tol
    namelist /control_settings/ t_step, t_end_ctrl, tol

    ! Current UT time tracker (updated each timestep)
    real(dp) :: current_sec_ut

    ! Read namelist once at startup
    integer :: iu, ios
    open(newunit=iu, file='../config.nml', status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'FATAL: cannot open ../config.nml'
    read(iu, nml=initial_conditions, iostat=ios)
    read(iu, nml=control_settings, iostat=ios)
    read(iu, nml=current_states, iostat=ios)
    read(iu, nml=atmosphere_input, iostat=ios)
    close(iu)
    if (ios /= 0) stop 'FATAL: error reading namelist /current_states/'

    !=============================== Main Body ===============================

    ! Initialize your state from the namelist vars
    state(1) = V
    state(2) = gamma * deg2rad
    state(3) = psi   * deg2rad
    state(4) = alt
    state(5) = lon   * deg2rad
    state(6) = lat   * deg2rad
    t = t0

    ! Initialize current UT time
    current_sec_ut = sec_ut

    ! Main integration loop
    i = 0
    do while (t < tend .and. state(4) > 0.0_dp)

        ! Re-read Cl/Cd each step
        if (Cd > 0.0_dp) then
            L_D = Cl / Cd
            beta_param = m / (Cd * Aref)
        else
            L_D = 0.0_dp ! Initial state
            beta_param = 1.0e10_dp  ! Very large = no drag (freefall)
        end if

        ! Get atmospheric properties based on the selected model
        call get_atmosphere(state(4), state(6)*rad2deg, state(5)*rad2deg, &
                           Temp, P, rho, a)
        ! Get gravity
        g   = EarthPointMass(state(4))

        ! Propagate one step (pass L/D & beta so CFD can drop in later)
        call rk4_nonplanar(t, state, state_new, dt, rho, g, m, beta_param, L_D, &
                           T_thrust, epsilon, sigma, REar, omega_earth)

        ! Update state
        state = state_new
        t = t + dt
        i = i + 1

        ! Update UT time for next iteration (sec_ut advances with simulation time)
        current_sec_ut = current_sec_ut + dt

        ! Handle day rollover (86400 seconds in a day)
        if (current_sec_ut >= 86400.0_dp) then
            current_sec_ut = current_sec_ut - 86400.0_dp
            iyd = iyd + 1  ! Increment day (simplified - doesn't handle year rollover)
        end if

    end do

    ! Update the state variables
    V  = state(1)
    gamma = state(2) * rad2deg
    psi= state(3) * rad2deg
    alt = state(4)
    lon = state(5) * rad2deg
    lat  = state(6) * rad2deg

    ! Update the cfd variables
    Mach = state(1) / a
    temperature = Temp + 273.15_dp
    Alpha = state(2) * rad2deg ! [rad]

    ! Update step and time
    step = step + 1.0_dp
    time = time + tend

    ! Update sec_ut to reflect propagation
    sec_ut = current_sec_ut

    open(newunit=iu, file='../config.nml', status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'FATAL: cannot open for writing'
    write(iu, nml=initial_conditions)
    write(iu, nml=control_settings)
    write(iu, nml=current_states)
    write(iu, nml=cfd_variables)
    write(iu, nml=atmosphere_input)
    close(iu)

contains

    !================================= Subroutines =================================
    subroutine get_atmosphere(alt_m, lat_deg, lon_deg, T_out, P_out, rho_out, a_out)
        implicit none
        real(dp), intent(in) :: alt_m, lat_deg, lon_deg
        real(dp), intent(out) :: T_out, P_out, rho_out, a_out

        real(dp):: T_K

        if (atm_model == 1) then
            ! NASA model
            call EarthAtmNASA(alt_m, T_out, P_out, rho_out, a_out)

        else
            ! NRLMSISE-00 model
            call EarthAtmNRLMSISE(iyd, current_sec_ut, alt_m, lat_deg, lon_deg, &
                                  f107a, f107, ap, rho_out, T_K, P_out, a_out)
            T_out = T_K - 273.15_dp ! Convert K to C

        end if

    end subroutine get_atmosphere

    subroutine rk4_nonplanar(t, state_in, state_out, dt, rho, g, m, beta, L_D, &
                             T_thrust, epsilon, sigma, REar, omega)
        implicit none
        real(dp), intent(in) :: t, dt, rho, g, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega
        real(dp), intent(in) :: state_in(6)
        real(dp), intent(out):: state_out(6)
        real(dp) :: k1(6), k2(6), k3(6), k4(6), state_t(6)
        real(dp) :: rho_t, g_t, h_t, dummy1, dummy2, dummy3

        call fun_nonplanar(t, state_in, k1, rho, g, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega)
        k1 = dt * k1

        state_t = state_in + 0.5_dp * k1
        h_t = state_t(4)
        call get_atmosphere(h_t, state_t(6)*rad2deg, state_t(5)*rad2deg, &
                    dummy1, dummy2, rho_t, dummy3)
        g_t   = EarthPointMass(h_t)
        call fun_nonplanar(t + 0.5_dp*dt, state_t, k2, rho_t, g_t, m, beta, L_D, T_thrust, epsilon, sigma, &
        REar, omega)
        k2 = dt * k2

        state_t = state_in + 0.5_dp * k2
        h_t = state_t(4)
        call get_atmosphere(h_t, state_t(6)*rad2deg, state_t(5)*rad2deg, &
                    dummy1, dummy2, rho_t, dummy3)
        g_t   = EarthPointMass(h_t)
        call fun_nonplanar(t + 0.5_dp*dt, state_t, k3, rho_t, g_t, m, beta, L_D, T_thrust, epsilon, sigma, &
        REar, omega)
        k3 = dt * k3

        state_t = state_in + k3
        h_t = state_t(4)
        call get_atmosphere(h_t, state_t(6)*rad2deg, state_t(5)*rad2deg, &
                    dummy1, dummy2, rho_t, dummy3)
        g_t   = EarthPointMass(h_t)
        call fun_nonplanar(t + dt, state_t, k4, rho_t, g_t, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega)
        k4 = dt * k4

        state_out = state_in + (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4) / 6.0_dp
    end subroutine rk4_nonplanar

    subroutine fun_nonplanar(t, state, dstate, rho, g, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega)
        implicit none
        real(dp), intent(in) :: t, rho, g, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega
        real(dp), intent(in) :: state(6)
        real(dp), intent(out):: dstate(6)
        real(dp) :: V, gamma, psi, h, theta, phi
        real(dp) :: r

        V     = state(1)
        gamma = state(2)
        psi   = state(3)
        h     = state(4)
        theta = state(5)
        phi   = state(6)

        r = REar + h

        ! --- Equations of Motion (non-planar, with aero & simple rotation terms) ---

        ! dV/dt
        dstate(1) = (T_thrust / m) * cos(epsilon) - rho * V**2.0_dp / (2.0_dp * beta) - g * sin(gamma) &
                    + omega**2.0_dp * r * cos(phi) * (sin(gamma) * cos(phi) - cos(gamma) * sin(phi) * sin(psi))

        ! dgamma/dt
        dstate(2) = T_thrust * sin(epsilon) / (V * m) * cos(sigma) + (V * cos(gamma)) / r &
                    + rho * V * L_D * cos(sigma) / (2.0_dp * beta) - (g * cos(gamma)) / V &
                    + 2.0_dp * omega * cos(phi) * cos(psi) &
                    + omega**2.0_dp * r * cos(phi) / V * ( cos(gamma) * cos(phi) + sin(gamma) * sin(psi) * sin(phi) )

        ! dpsi/dt
        dstate(3) = T_thrust * sin(epsilon) * sin(sigma) / (V * m * cos(gamma)) &
                    + rho * V * L_D * sin(sigma) / (2.0_dp * beta * cos(gamma)) &
                    + (V * cos(gamma) * cos(psi) * tan(phi))/r &
                    + 2.0_dp * omega * (tan(gamma) * cos(phi) * sin(psi) - sin(phi) ) &
                    - omega**2.0_dp * r * sin(phi) * cos(phi) * cos(psi)/(V * cos(gamma))

        ! dh/dt
        dstate(4) = V * sin(gamma)

        ! dtheta/dt
        dstate(5) = V * cos(gamma) * cos(psi) / (r * cos(phi))

        ! dphi/dt
        dstate(6) = V * cos(gamma) * sin(psi) / r

    end subroutine fun_nonplanar

end program ReentryNonPlanar3DOF

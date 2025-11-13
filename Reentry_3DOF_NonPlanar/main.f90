!> Non-Planar 3-DOF Reentry Trajectory Propagator
!> Author: Kevin Tang
program ReentryNonPlanar3DOF
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use mod_atmosphere, only: EarthAtmNASA
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
    namelist /trajectory_input/ step, time, V, gamma, psi, alt, lon, lat, Cl, Cd
    namelist /cfd_input/ Mach, temperature, Alpha

    ! Read namelist once at startup
    integer :: iu, ios
    open(newunit=iu, file='../config.nml', status='old', action='read', iostat=ios)
    if (ios /= 0) stop 'FATAL: cannot open ../config.nml'
    read(iu, nml=trajectory_input, iostat=ios)
    close(iu)
    if (ios /= 0) stop 'FATAL: error reading namelist /trajectory_input/'

    ! Initialize your state from the namelist vars
    state(1) = V
    state(2) = gamma * deg2rad
    state(3) = psi   * deg2rad
    state(4) = alt
    state(5) = lon   * deg2rad
    state(6) = lat   * deg2rad
    t = t0

    ! Main integration loop
    i = 0
    do while (t < tend .and. state(4) > 0.0_dp)

        ! Re-read Cl/Cd each step
        L_D = Cl / Cd ! lift-to-drag ratio
        beta_param = m / (Cd * Aref) ! ballistic coefficient m/(Cd*A)

        ! Environment
        call EarthAtmNASA(state(4), Temp, P, rho, a)
        g   = EarthPointMass(state(4))

        ! Propagate one step (pass L/D & beta so CFD can drop in later)
        call rk4_nonplanar(t, state, state_new, dt, rho, g, m, beta_param, L_D, &
                           T_thrust, epsilon, sigma, REar, omega_earth)

        ! Update state
        state = state_new
        t = t + dt
        i = i + 1

    end do

    ! Update the state variables
    V       = state(1)
    gamma = state(2) * rad2deg
    psi= state(3) * rad2deg
    alt     = state(4)
    lon = state(5) * rad2deg
    lat  = state(6) * rad2deg

    ! Update the cfd variables
    Mach = state(1) / a
    temperature = Temp + 273.15_dp
    Alpha = state(2) * rad2deg ! [rad]

    ! Update step and time
    step = step + 1.0_dp
    time = time + tend

    open(newunit=iu, file='../config.nml', status='replace', action='write', iostat=ios)
    if (ios /= 0) stop 'FATAL: cannot open for writing'
    write(iu, nml=trajectory_input)
    write(iu, nml=cfd_input)

close(iu)

contains

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
        call EarthAtmNASA(h_t, dummy1, dummy2, rho_t, dummy3)
        g_t   = EarthPointMass(h_t)
        call fun_nonplanar(t + 0.5_dp*dt, state_t, k2, rho_t, g_t, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega)
        k2 = dt * k2

        state_t = state_in + 0.5_dp * k2
        h_t = state_t(4)
        call EarthAtmNASA(h_t, dummy1, dummy2, rho_t, dummy3)
        g_t   = EarthPointMass(h_t)
        call fun_nonplanar(t + 0.5_dp*dt, state_t, k3, rho_t, g_t, m, beta, L_D, T_thrust, epsilon, sigma, REar, omega)
        k3 = dt * k3

        state_t = state_in + k3
        h_t = state_t(4)
        call EarthAtmNASA(h_t, dummy1, dummy2, rho_t, dummy3)
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

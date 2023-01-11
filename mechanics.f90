!-----------------------------------------------------------------------
!Module: mechanics
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! 
!! This module takes care of the physics of the current problem. 
!! This is where we will calculate r or the distance between planets and the 
!! primary. We will also calculate the total energy due to the system which 
!! should be more or less conserved. 
!! More importantly this is where we will define our 8 coupled second order 
!! ordinary differential equations that describe the positions of the planets.
!!----------------------------------------------------------------------
!! Included subroutines: calculate_r, calculate_energy
!!
!!----------------------------------------------------------------------
!! Included functions: planets_ode
!!
!-----------------------------------------------------------------------
module mechanics
use types
implicit none
private
public :: planets_ode, calculate_energy

contains


!-----------------------------------------------------------------------
!! function: planets_ode
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This function contains the 8 coupled ODEs that describe the position 
!! of the planets. It returns an f array after the ODEs have been solved 
!! for a specific moment in time. We need to be careful here to make sure 
!! we define this function the same way our interface is defined!
!!----------------------------------------------------------------------
!! Input: r     Array containing the initial or current positions and velocities 
!!              of the planets around the primary. 
!!        t     Array containing time in seconds of the system 
!!        work  Array containing masses of the planets and the primary
!!
!!----------------------------------------------------------------------
!! Output: f    Array containing calculated velocities and accelerations 
!!              of the system.
!!
!-----------------------------------------------------------------------
function planets_ode(r, t, work) result(f)
    implicit none
    real(dp), intent(in) :: r(:), t, work(:)
    real(dp), allocatable :: f(:) 
    real(dp) :: x1, x2, y1, y2, vx1, vy1, vx2, vy2, primary_m, planet_m1, planet_m2 
    real(dp) :: r1, r2, r12 
    allocate(f(1:size(r))) 
    !Setup initial positions and velocities
    x1 = r(1) 
    y1 = r(2) 
    x2 = r(3) 
    y2 = r(4) 
    vx1 = r(5) 
    vy1 = r(6) 
    vx2 = r(7) 
    vy2 = r(8) 
    
    primary_m = work(1) 
    planet_m1 = work(2) 
    planet_m2 = work(3) 
    !Need to calculate radius for the ODE before calculating anything
    call calculate_r(x1, y1, x2, y2, r1, r2, r12) 
    !f(1-4) are the velocities and f(5-8) are the accelerations that are calculated 
    !through newton's laws of gravity.
    !r(1-4) are the initial positions and r(5-8) are the initial velocities 
    !So our f array should be initialized as follows: 
    f(1) = r(5) 
    f(2) = r(6) 
    f(3) = r(7) 
    f(4) = r(8) 
    
    f(5) = (-primary_m * x1/(r1**3)) + (-planet_m2 * (x1 - x2)/(r12**3)) 
    f(6) = (-primary_m * y1/(r1**3)) + (-planet_m2 * (y1 - y2)/(r12**3)) 
    f(7) = (-primary_m * x2/(r2**3)) + (-planet_m1 * (x2 - x1)/(r12**3)) 
    f(8) = (-primary_m * y2/(r2**3)) + (-planet_m1 * (y2 - y1)/(r12**3))
    ! This is the function that will be sent to 
    ! solve_runge_kutta_4 as an argument.

    
end function planets_ode

!-----------------------------------------------------------------------
!! Subroutine: calculate_energy
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This subroutine calculates the energy of the system(which should be conserved) 
!! by calculating the kinetic terms and potential terms and adding them together.
!!----------------------------------------------------------------------
!! Input: solution      Array containing the solution to the 8 coupled 
!!                      ODEs due to Newton's Law of Gravity. 
!!        work          Array containing masses of the planets and the primary. 
!!        n_steps       Integer value containing number of points between initial 
!!                      final time.
!!
!!----------------------------------------------------------------------
!! Output: energy       Array containing energy of the system 
!!
!-----------------------------------------------------------------------
subroutine calculate_energy(solution, energy, n_steps, work_array)
    implicit none 
    real(dp), intent(in) :: solution(:, :), work_array(:) 
    real(dp), allocatable, intent(out) :: energy(:) 
    real(dp) :: kinetic_m1, kinetic_m2, fg_m1, fg_m2, fg_m12 
    real(dp) :: r1, r2, r12, x1, y1, x2, y2, primary_m, m1, m2, vx1, vy1, vx2, vy2
    integer, intent(in) :: n_steps 
    integer :: i 
    
    primary_m = work_array(1) 
    m1 = work_array(2) 
    m2 = work_array(3)
    
    allocate(energy(1:n_steps)) 
    do i = 1, n_steps 
        !Solution array is setup such that (1-4, :) are the positions 
        !and (5-8, : are the velocties. 
        x1 = solution(1, i) 
        y1 = solution(2, i) 
        x2 = solution(3, i) 
        y2 = solution(4, i) 
        vx1 = solution(5, i) 
        vy1 = solution(6, i) 
        vx2 = solution(7, i) 
        vy2 = solution(8, i) 
        !First we should calculate the r term. 
        call calculate_r(x1, y1, x2, y2, r1, r2, r12) 
        !To improve readability lets calculate each kinetic and potential 
        !term and then add them. 
        kinetic_m1 = 0.5_dp * m1 * (vx1**2 + vy1**2) 
        kinetic_m2 = 0.5_dp * m2 * (vx2**2 + vy2**2) 
        !Potential should be negative here. 
        fg_m1 = -m1 * primary_m/r1 
        fg_m2 = -m2 * primary_m/r2 
        fg_m12 = -m1 * m2/r12
        energy(i) = kinetic_m1 + kinetic_m2 + fg_m1 + fg_m2 + fg_m12
    end do 
    
end subroutine calculate_energy



!-----------------------------------------------------------------------
!! Subroutine: calculate_r
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!!
!! This module is more straightforward. We're defining this just to make the 
!! program easier to read. It takes in the positions of the planets to calculate the 
!! distances between them and the primaries. 
!!----------------------------------------------------------------------
!! Input: x1        x coordinate of planet 1 
!!        y1        y coordinate of planet 1 
!!        x2        x coordinate of planet 2 
!!        y2        y coordinate of planet 2 
!!
!!----------------------------------------------------------------------
!! Output: r1       distance from planet 1 to the primary 
!!         r2       distance from planet 2 to the primary 
!!         r12      distance from planet 1 to planet 2
!!
!-----------------------------------------------------------------------
subroutine calculate_r(x1, y1, x2, y2, r1, r2, r12)   
implicit none 
real(dp), intent(out) :: r1, r2, r12 
real(dp), intent(in) :: x1, y1, x2, y2 
!This is more or less straightforward. We use the Pythagorean Theorem to calculate 
!distance since it is a vector quantity.
r1 = sqrt(x1**2 + y1**2) 
r2 = sqrt(x2**2 + y2**2) 
r12 = sqrt((x1 - x2)**2 + (y1 - y2)**2)
end subroutine calculate_r
end module mechanics

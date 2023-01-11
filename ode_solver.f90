!-----------------------------------------------------------------------
!Module: ode_solver
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! 
!! This is the "math module." This module will solve the system of ODEs 
!! for this problem. In fact it is generalized with this interface so that 
!! it will solve ANY system of ODEs with fourth order Runge Kutta. To do 
!! that we defined the interface below.
!!----------------------------------------------------------------------
!! Included subroutines: solve_runge_kutta_4
!!
!!----------------------------------------------------------------------
!! Included interface: func 
!!
!-----------------------------------------------------------------------
module ode_solver
use types
implicit none
private

public :: solve_runge_kutta_4
!We want to send solve_runge_kutta_4 a function as an argument. 
!In order to do this we need to declare this interface which is defined exactly 
!like planets ODE to generalize this module for any system of ODEs.
interface
    function func(r, t, work) result(f) 
    use types, only : dp
    implicit none 
    real(dp), intent(in) :: r(:), work(:), t 
    real(dp), allocatable :: f(:)
    end function func
end interface

contains

!-----------------------------------------------------------------------
!! Subroutine: solve_runge_kutta_4
!-----------------------------------------------------------------------
!! By: Nathan Crawford 
!! 
!! This subroutine solves a system of ODEs using the fourth order Runge 
!! Kutta method. We define the k1, k2, k3 and k4 as arrays to avoid 
!! any unnecessary complexity. 
!!----------------------------------------------------------------------
!! Input: f                     This is a function or the system of ODEs as an argument 
!!        initial_condition     Array containing the initial conditions 
!!        final_time            maximum time 
!!        n_steps               number of points between initial and final time 
!!        work_array            array containing masses of planets 
!!        t                     array containing points in time 
!!
!!----------------------------------------------------------------------
!! Output: solution     two dimensional array containing positions and 
!!                      velocities of the 2 planets at every time t.
!!
!-----------------------------------------------------------------------
subroutine solve_runge_kutta_4(f, initial_condition, final_time, n_steps, work_array, t, solution)
    implicit none
    procedure(func) :: f  
    real(dp) :: initial_time, h, t_sol 
    real(dp), allocatable :: r_sol(:), k1(:), k2(:), k3(:), k4(:), dr(:) 
    real(dp), intent(in) :: final_time, initial_condition(:), work_array(:) 
    real(dp), allocatable, intent(out) :: t(:), solution(:, :) 
    integer :: i, n_variables, n_steps  
    
    initial_time = 0.0_dp 
    h = (final_time - initial_time)/n_steps 
    !n_variables are the amount of dependent variables 
    n_variables = size(initial_condition)
    
    allocate(t(n_steps)) 
    allocate(solution(n_variables, n_steps)) 
    allocate(r_sol(n_variables)) 
    !making k1, k2, k3 and k4 removes any unnecessary complexity
    allocate(k1(n_variables), k2(n_variables), k3(n_variables), k4(n_variables))
    !initialize all arrays
    t_sol = initial_time
    r_sol = initial_condition
    t = initial_time 
    
    do i = 1, n_steps 
        
        k1 = h*f(r_sol, t(i), work_array) 
        k2 = h*f(r_sol + 0.5_dp * k1, t(i) + 0.5_dp * h, work_array) 
        k3 = h*f(r_sol + 0.5_dp * k2, t(i) + 0.5_dp * h, work_array) 
        k4 = h*f(r_sol + k3, t(i) + h, work_array) 
        dr = (k1 + 2 * k2 + 2 * k3 + k4)/6 
        r_sol = r_sol + dr 
        t_sol = t_sol + h
        t(i) = t_sol
        solution(:, i) = r_sol 
    end do 
    
end subroutine solve_runge_kutta_4

    
end module ode_solver

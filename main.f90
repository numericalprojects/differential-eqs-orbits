! Program: planets
! By: Nathan Crawford 
! This program solves the equation of motion for 
! 2 planets and a primary in two dimensions(x and y). 
! The equations of motion involve Newton's laws of gravity. 
! The program has a read_write module which takes input for 
! the initial and final time, number of points in said time interval, 
! masses and initial positions and velocities of the planets and the primary, 
! and the name of the output file. 
! It should be noted this program also uses namelists so if you want to run this 
! program with different outputs then you can change the namelist file. 
! The program also has a solve runge kutta module which solves any system of 
! ODEs with fourth order Runge Kutta method. 
! The program has a mechanics module which takes care of the physics of this 
! specific problem. 
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
program planets

use types
use read_write, only : read_input, write_results
use ode_solver, only : solve_runge_kutta_4
use mechanics, only : calculate_energy, planets_ode

implicit none

real(dp) :: work_array(1:3), initial_condition(1:8)
real(dp) :: final_time
integer :: n_steps
real(dp), allocatable :: time(:), solution(:,:), energy(:)
character(len=1024) :: output_file

 call read_input(work_array, initial_condition, final_time, n_steps, output_file)

 call solve_runge_kutta_4(planets_ode, initial_condition, final_time, n_steps, work_array, time, solution)
 call calculate_energy(solution, energy, n_steps, work_array)

 call write_results(output_file, energy, time, solution, n_steps)


end program planets

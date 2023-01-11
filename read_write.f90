!-----------------------------------------------------------------------
!Module: read_write
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! 
!! This module reads input from a namelist or the default values and 
!! then writes positions of planets as well as energy as a function of time
!!----------------------------------------------------------------------
!! Included subroutines: read_input, write_results
!!
!!----------------------------------------------------------------------
!! Included functions:
!!
!-----------------------------------------------------------------------
module read_write
use types
implicit none

private
public :: read_input, write_results

contains

!-----------------------------------------------------------------------
!! Subroutine: read_input
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! 
!! This subroutine reads input from a namelist or uses default values
!!----------------------------------------------------------------------
!! Output: work_array   Array containing masses of planets and primary 
!!         initial_condition    Array containing initial conditions 
!!         final_time   maximum time 
!!         n_steps      number of points between initial and final time 
!!         output_file  file for results to be written into
!!
!-----------------------------------------------------------------------
subroutine read_input(work_array, initial_condition, final_time, n_steps, output_file)
    implicit none
    real(dp), intent(out) :: work_array(1:3)
    real(dp), intent(out) :: initial_condition(1:8)
    real(dp), intent(out) :: final_time
    integer, intent(out) :: n_steps
    character(len=*), intent(out) :: output_file
    real(dp) :: primary_mass, planet_mass_1, planet_mass_2
    real(dp) :: initial_pos_1(1:2), initial_pos_2(1:2)
    real(dp) :: initial_vel_1(1:2), initial_vel_2(1:2) 
    logical :: file_exists 
    integer :: n_arguments, file_unit, ierror
    character(len=200) :: namelist_file


    namelist /masses/ primary_mass, planet_mass_1, planet_mass_2
    namelist /initial_conditions/ initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2
    namelist /solution_parameters/ final_time, n_steps
    namelist /output/ output_file


    ! Set default values
    
    
    
     print *, "Greetings user, this program calculates the positions of 2 planets orbitting a primary" 
     print *, "according to Newton's Law of Gravity."
     primary_mass  = 1.0_dp
     planet_mass_1 = 1.0_dp
     planet_mass_2 = 0.0_dp
     initial_pos_1 = [1.0_dp, 0.0_dp]
     initial_pos_2 = [-1.0_dp, 1.0_dp]
     initial_vel_1 = [0.0_dp, 1.0_dp]
     initial_vel_2 = [0.0_dp, -1.0_dp]
     final_time = 30.0_dp
     n_steps = 30000
     output_file = 'planet_motion.dat'

  
    
    n_arguments = command_argument_count()
 

    if (n_arguments == 1) then
        call get_command_argument(1, namelist_file)
        inquire(file = trim(namelist_file), exist = file_exists)
        if (file_exists) then
            open(newunit=file_unit, file=namelist_file)
            read(file_unit, nml=masses, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading masses namelist"
                stop
            endif
            read(file_unit, nml=initial_conditions, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading initial_conditions namelist"
                stop
            endif
            read(file_unit, nml=solution_parameters, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading solution_parameters namelist"
                stop
            endif
            read(file_unit, nml=output, iostat=ierror)
            if(ierror /= 0) then
                print*, "Error reading output namelist"
                stop
            endif
        else
            print*, namelist_file, 'not found'
            stop
        endif
    elseif (n_arguments /= 0) then
        print*, 'Incorrect number of arguments. Program takes either 0 or 1 argument only'
        print*, 'See details in README.md'
        stop
    endif
    
    work_array = [primary_mass, planet_mass_1, planet_mass_2]
    initial_condition = [initial_pos_1, initial_pos_2, initial_vel_1, initial_vel_2]

end subroutine read_input

!-----------------------------------------------------------------------
!! Subroutine: write_results
!-----------------------------------------------------------------------
!! By: Nathan Crawford
!! 
!! This subroutine writes the time, positions of the planets and energy 
!! of the system into a specified output file.
!!----------------------------------------------------------------------
!! Input: energy    Array containing energy of the system
!!        t         Array containing each time t 
!!        solution  2 dimensional array containing positions and velocities 
!!                  of the planets 
!!        n_steps   number of points between initial and final time
!!
!-----------------------------------------------------------------------
subroutine write_results(output_file, energy, t, solution, n_steps)
    implicit none 
    real(dp), intent(in) :: energy(:), t(:), solution(:, :) 
    integer, intent(in) :: n_steps 
    integer :: i, unit
    character(len=*) :: output_file
    print * 
    print *, "Writing results into ", output_file 
    
    open(newunit = unit, file = output_file) 
    write(unit, *) 'time ', 'x1 ', 'y1 ', 'x2 ', 'y2 ', 'Energy' 
    do i = 1, n_steps 
        write(unit, *) t(i), solution(1, i), solution(2, i), solution(3, i), solution(4, i), energy(i) 
    end do 
    close(unit) 
    print *
    print *, "Results written into ", output_file 
   
end subroutine write_results


end module read_write

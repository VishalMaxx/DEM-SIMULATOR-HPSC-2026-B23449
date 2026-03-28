program particle_simulator
    use physics_data
    implicit none
    integer :: step, i
    real(8) :: t, final_t = 3.0d0
    real(8), dimension(4) :: gammas = (/10.0d0, 50.0d0, 100.0d0, 250.0d0/)

    print *, "Starting Damping Effect Study (Section 19.1)..."
    N = 1 
    dt = 0.001d0

    open(unit=30, file="damping_results.txt", status="replace")
    write(30, '(A10, A10, A15, A15)') "# Gamma", "Time", "Height_Z", "Velocity_Z"

    do i = 1, 4
        gamma_n = gammas(i)
        call allocate_particles()
        
        ! Initial State: High drop to see multiple bounces
        z(1) = 10.0d0
        vz(1) = 0.0d0
        
        do step = 1, int(final_t / dt)
            t = step * dt
            call apply_gravity()
            call compute_wall_contacts()
            call update_particles()
            
            ! Record every 10 steps to keep the file size reasonable
            if (mod(step, 10) == 0) then
                write(30, '(F10.1, F10.4, F15.6, F15.6)') gamma_n, t, z(1), vz(1)
            end if
        end do
        
        print '(A,F6.1,A)', "Finished simulation for Gamma = ", gamma_n, "..."
    end do

    close(30)
    print *, "Damping data saved to damping_results.txt"
end program particle_simulator

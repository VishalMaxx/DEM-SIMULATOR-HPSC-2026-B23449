program particle_simulator
    use physics_data
    use omp_lib
    implicit none
    integer :: i, step, total_steps = 200
    real(8) :: start_t, end_t, rand_val
    integer :: threads

    ! Request the current number of threads from the environment
    threads = omp_get_max_threads()
    N = 5000 
    
    print *, "Performance Test: N =", N, "| Threads =", threads
    call allocate_particles()

    ! Random initialization for a dense system
    do i = 1, N
        call random_number(rand_val); x(i) = rand_val * (Lx - 1.0) + 0.5
        call random_number(rand_val); y(i) = rand_val * (Ly - 1.0) + 0.5
        call random_number(rand_val); z(i) = rand_val * (Lz - 1.0) + 0.5
    end do

    start_t = omp_get_wtime()

    do step = 1, total_steps
        call apply_gravity()
        call compute_wall_contacts()
        call compute_particle_contacts()
        call update_particles()
    end do

    end_t = omp_get_wtime()
    print '(A,F10.4,A)', "Execution Time: ", end_t - start_t, " seconds"
    
    deallocate(x, y, z, vx, vy, vz, fx, fy, fz, mass, radius)
end program particle_simulator

module physics_data
    implicit none
    save
    ! --- Neighbour Search Parameters ---
    real(8) :: cell_size
    integer :: n_cells_x, n_cells_y, n_cells_z
    integer, allocatable :: head(:)      ! Stores the first particle in each cell
    integer, allocatable :: next_p(:)    ! Links to the next particle in the same cell
    ! --- Simulation Parameters ---
    integer :: N = 100                     
    real(8) :: dt = 0.001                  
    real(8) :: g = -9.81                   
    
    ! --- Physical Constants for Contact ---
    real(8) :: kn = 10000.0                
    real(8) :: gamma_n = 50.0              

    ! --- Box Dimensions ---
    real(8) :: Lx = 10.0, Ly = 10.0, Lz = 10.0

    ! --- Particle Arrays ---
    real(8), allocatable :: x(:), y(:), z(:)       
    real(8), allocatable :: vx(:), vy(:), vz(:)    
    real(8), allocatable :: fx(:), fy(:), fz(:)    
    real(8), allocatable :: mass(:), radius(:)     

! ==========================================================
contains  ! <--- THIS IS THE SECTION YOU ARE LOOKING FOR!
! ==========================================================

    ! Function to calculate Total Kinetic Energy [cite: 128]
    function compute_kinetic_energy() result(ke)
        real(8) :: ke
        integer :: i
        ke = 0.0
        do i = 1, N
            ke = ke + 0.5d0 * mass(i) * (vx(i)**2 + vy(i)**2 + vz(i)**2)
        end do
    end function compute_kinetic_energy

    subroutine allocate_particles()
        allocate(x(N), y(N), z(N), vx(N), vy(N), vz(N))
        allocate(fx(N), fy(N), fz(N), mass(N), radius(N))
        
        ! Grid setup: Cell size must be at least the max particle diameter
        cell_size = 1.1d0 * (2.0d0 * 0.5d0) ! slightly larger than diameter
        n_cells_x = ceiling(Lx / cell_size)
        n_cells_y = ceiling(Ly / cell_size)
        n_cells_z = ceiling(Lz / cell_size)
        
        allocate(head(n_cells_x * n_cells_y * n_cells_z))
        allocate(next_p(N))
        
        vx = 0.0; vy = 0.0; vz = 0.0
        fx = 0.0; fy = 0.0; fz = 0.0
        mass = 1.0; radius = 0.5
    end subroutine allocate_particles

    subroutine apply_gravity()
        integer :: i
        do i = 1, N
            fx(i) = 0.0; fy(i) = 0.0; fz(i) = mass(i) * g
        end do
    end subroutine apply_gravity

    subroutine compute_wall_contacts()
        integer :: i
        real(8) :: delta, force_n
        do i = 1, N
            delta = radius(i) - z(i) 
            if (delta > 0.0) then
                force_n = kn * delta - gamma_n * vz(i) 
                if (force_n < 0.0) force_n = 0.0 
                fz(i) = fz(i) + force_n
            end if
        end do
    end subroutine compute_wall_contacts

    subroutine compute_particle_contacts()
        use omp_lib
        integer :: i, j
        real(8) :: dx, dy, dz, dist, overlap, nx, ny, nz, v_rel_n, force_n
        !$omp parallel do private(j, dx, dy, dz, dist, overlap, nx, ny, nz, v_rel_n, force_n)
        do i = 1, N - 1
            do j = i + 1, N
                dx = x(j) - x(i); dy = y(j) - y(i); dz = z(j) - z(i)
                dist = sqrt(dx**2 + dy**2 + dz**2)
                overlap = (radius(i) + radius(j)) - dist
                if (overlap > 0.0) then
                    nx = dx / dist; ny = dy / dist; nz = dz / dist
                    v_rel_n = (vx(j) - vx(i)) * nx + (vy(j) - vy(i)) * ny + (vz(j) - vz(i)) * nz
                    force_n = kn * overlap - gamma_n * v_rel_n
                    if (force_n < 0.0) force_n = 0.0 
                    !$omp atomic
                    fx(i) = fx(i) - force_n * nx
                    !$omp atomic
                    fy(i) = fy(i) - force_n * ny
                    !$omp atomic
                    fz(i) = fz(i) - force_n * nz
                    !$omp atomic
                    fx(j) = fx(j) + force_n * nx
                    !$omp atomic
                    fy(j) = fy(j) + force_n * ny
                    !$omp atomic
                    fz(j) = fz(j) + force_n * nz
                end if
            end do
        end do
        !$omp end parallel do
    end subroutine compute_particle_contacts

    subroutine update_particles()
        integer :: i
        do i = 1, N
            vx(i) = vx(i) + (fx(i) / mass(i)) * dt
            vy(i) = vy(i) + (fy(i) / mass(i)) * dt
            vz(i) = vz(i) + (fz(i) / mass(i)) * dt
            x(i) = x(i) + vx(i) * dt
            y(i) = y(i) + vy(i) * dt
            z(i) = z(i) + vz(i) * dt
        end do
    end subroutine update_particles
    
    subroutine build_neighbor_list()
        integer :: i, cx, cy, cz, cell_idx
        
        head = 0  ! Reset the grid
        next_p = 0
        
        do i = 1, N
            ! Find which cell the particle is in
            cx = int(x(i) / cell_size) + 1
            cy = int(y(i) / cell_size) + 1
            cz = int(z(i) / cell_size) + 1
            
            ! Stay within bounds
            cx = max(1, min(cx, n_cells_x))
            cy = max(1, min(cy, n_cells_y))
            cz = max(1, min(cz, n_cells_z))
            
            ! 1D index for the 3D grid
            cell_idx = cx + (cy-1)*n_cells_x + (cz-1)*n_cells_x*n_cells_y
            
            ! Link the particle into the list
            next_p(i) = head(cell_idx)
            head(cell_idx) = i
        end do
    end subroutine build_neighbor_list
end module physics_data

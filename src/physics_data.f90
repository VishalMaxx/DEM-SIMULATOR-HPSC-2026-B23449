module physics_data
    use omp_lib
    implicit none
    save

    ! --- Parameters ---
    integer :: N = 5000                     
    real(8) :: dt = 0.001                  
    real(8) :: g = -9.81                   
    real(8) :: kn = 10000.0                
    real(8) :: gamma_n = 50.0              
    real(8) :: Lx = 10.0, Ly = 10.0, Lz = 10.0

    ! --- Particle Arrays ---
    real(8), allocatable :: x(:), y(:), z(:)       
    real(8), allocatable :: vx(:), vy(:), vz(:)    
    real(8), allocatable :: fx(:), fy(:), fz(:)    
    real(8), allocatable :: mass(:), radius(:)     

    ! --- Neighbour Search Arrays (Section 18) ---
    real(8) :: cell_size
    integer :: n_cells_x, n_cells_y, n_cells_z
    integer, allocatable :: head(:)      
    integer, allocatable :: next_p(:)    

contains

    subroutine allocate_particles()
        if (allocated(x)) then
            deallocate(x, y, z, vx, vy, vz, fx, fy, fz, mass, radius, head, next_p)
        end if
        allocate(x(N), y(N), z(N), vx(N), vy(N), vz(N))
        allocate(fx(N), fy(N), fz(N), mass(N), radius(N))
        
        cell_size = 1.1d0 * (2.0d0 * 0.5d0) 
        n_cells_x = ceiling(Lx / cell_size)
        n_cells_y = ceiling(Ly / cell_size)
        n_cells_z = ceiling(Lz / cell_size)
        
        allocate(head(n_cells_x * n_cells_y * n_cells_z))
        allocate(next_p(N))
        
        vx = 0.0; vy = 0.0; vz = 0.0; mass = 1.0; radius = 0.5
        fx = 0.0; fy = 0.0; fz = 0.0
    end subroutine allocate_particles

    subroutine apply_gravity()
        integer :: i
        do i = 1, N
            fx(i) = 0.0
            fy(i) = 0.0
            fz(i) = mass(i) * g
        end do
    end subroutine apply_gravity

    subroutine build_neighbor_list()
        integer :: i, cx, cy, cz, cell_idx
        head = 0  
        next_p = 0
        do i = 1, N
            cx = max(1, min(int(x(i) / cell_size) + 1, n_cells_x))
            cy = max(1, min(int(y(i) / cell_size) + 1, n_cells_y))
            cz = max(1, min(int(z(i) / cell_size) + 1, n_cells_z))
            cell_idx = cx + (cy-1)*n_cells_x + (cz-1)*n_cells_x*n_cells_y
            next_p(i) = head(cell_idx)
            head(cell_idx) = i
        end do
    end subroutine build_neighbor_list

    subroutine compute_particle_contacts()
        integer :: i, j, cx, cy, cz, nx_c, ny_c, nz_c, neigh_idx, dx_c, dy_c, dz_c
        real(8) :: dx, dy, dz, dist, overlap, nx, ny, nz, v_rel_n, force_n

        !$omp parallel do private(j, cx, cy, cz, nx_c, ny_c, nz_c, neigh_idx, &
        !$omp& dx_c, dy_c, dz_c, dx, dy, dz, dist, overlap, nx, ny, nz, v_rel_n, force_n)
        do i = 1, N
            cx = int(x(i) / cell_size) + 1
            cy = int(y(i) / cell_size) + 1
            cz = int(z(i) / cell_size) + 1

            do dx_c = -1, 1
                do dy_c = -1, 1
                    do dz_c = -1, 1
                        nx_c = cx + dx_c; ny_c = cy + dy_c; nz_c = cz + dz_c
                        if (nx_c < 1 .or. nx_c > n_cells_x) cycle
                        if (ny_c < 1 .or. ny_c > n_cells_y) cycle
                        if (nz_c < 1 .or. nz_c > n_cells_z) cycle

                        neigh_idx = nx_c + (ny_c-1)*n_cells_x + (nz_c-1)*n_cells_x*n_cells_y
                        j = head(neigh_idx)
                        do while (j > 0)
                            if (i < j) then
                                dx = x(j)-x(i); dy = y(j)-y(i); dz = z(j)-z(i)
                                dist = sqrt(dx**2 + dy**2 + dz**2)
                                overlap = (radius(i) + radius(j)) - dist
                                if (overlap > 0.0) then
                                    nx = dx/dist; ny = dy/dist; nz = dz/dist
                                    v_rel_n = (vx(j)-vx(i))*nx + (vy(j)-vy(i))*ny + (vz(j)-vz(i))*nz
                                    force_n = max(0.0d0, kn * overlap - gamma_n * v_rel_n)
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
                            end if
                            j = next_p(j)
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine compute_particle_contacts

    subroutine compute_wall_contacts()
        integer :: i
        real(8) :: delta, force_n
        do i = 1, N
            delta = radius(i) - z(i) 
            if (delta > 0.0) then
                force_n = max(0.0d0, kn * delta - gamma_n * vz(i))
                fz(i) = fz(i) + force_n
            end if
        end do
    end subroutine compute_wall_contacts

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

end module physics_data

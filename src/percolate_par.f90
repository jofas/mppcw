program percolate
  !
  ! Program percolate.
  !
  ! This program searches clusters in a L x L matrix that
  ! percolate. A cell is either full (0) or empty
  ! (1 ... |empty cells|). A cell has 4 neighbor cells and
  ! builds a cluster with them (and therefore with their
  ! neighbors' neighbors, making it a recursive relation-
  ! ship). A cluster percolates, if such a cluster starts
  ! from the left most column and finds its way to the
  ! right most column.
  !

  use io
  use uni
  use cart_comm

  use mpi_f08

  implicit none

  type(CLIResults) :: cli
  type(CartComm) :: comm
  integer :: L, M, N
  integer, dimension(:, :), allocatable :: map, chunk

  real :: start_t, end_t

  call init()

  start_t = mpi_wtime()

  call scatter(comm, map, chunk)

  call cluster()

  call gather(comm, map, chunk)

  end_t = mpi_wtime()

  if(comm%rank == 0) &
    print *, comm%w_size, ",", L, ",", end_t - start_t

  call finalize()

contains

  subroutine init()
    !
    ! Initialize all components (random number generator,
    ! MPI session, map (only on root process), chunk,
    ! 2d cartesion communicator and variables describing
    ! the format of the map/the chunk (L, M, N).
    !
    cli = read_from_cli()

    call rinit(cli%seed)

    L = cli%matrix_dimension

    call mpi_init()

    call init_cart_comm(comm, L, mpi_comm_world)

    M = comm%M
    N = comm%N

    if (comm%rank == 0) then
      call init_map()
    end if

    call init_chunk()
  end


  subroutine finalize()
    !
    ! Write .pgm output file and close the MPI session.
    !
    if (comm%rank == 0) then
      call pgm_write(cli%pgm_file_path, map, &
        cli%print_n_clusters)
    end if

    call mpi_finalize()
  end


  subroutine init_map()
    !
    ! Initialize the map.
    !
    ! Only called by the root process.
    !
    integer :: i, j, free_cell_count

    allocate(map(L, L))
    map(:, :) = 0

    free_cell_count = 0
    do i = 1, L
      do j = 1, L
        if (random_uniform() > cli%density_of_filled_cells) then
          free_cell_count = free_cell_count + 1
          map(i, j) = free_cell_count
        end if
      end do
    end do
  end


  subroutine init_chunk()
    allocate(chunk(0:M + 1, 0:N + 1))
    chunk(:, :) = 0
  end


  subroutine cluster()
    !
    ! Performs the clustering operation. First the halos of
    ! the chunks are swapped with its neighbors. Then a
    ! cluster step is performed (every cell of the inner
    ! chunk is set to its biggest neighbor).
    !
    ! Afterwards the sum of every inner cell of the chunk
    ! is reduced over every chunk and broadcasted to every
    ! process (mpi_allreduce). If the gathered sum did not
    ! change this cluster step, clustering is finished.
    !
    integer :: i, sum_, sum_old

    i = 1
    sum_old = 0
    do
      sum_ = 0

      call halo_swap(comm, chunk)

      call cluster_step()

      call mpi_allreduce(sum(chunk(1:M, 1:N)), sum_, 1, &
        mpi_integer, mpi_sum, comm%comm)

      if (sum_ == sum_old) exit

      sum_old = sum_

      i = i + 1
    end do
  end


  subroutine cluster_step()
    !
    ! Sets every inner cell of chunk to its biggest neighbor.
    !
    integer :: i, j

    forall(j=1:N, i=1:M, chunk(i, j) /= 0)
      chunk(i, j) = max( chunk(i, j + 1) &
                       , chunk(i, j - 1) &
                       , chunk(i + 1, j) &
                       , chunk(i - 1, j) &
                       , chunk(i, j) )
    end forall
  end
end

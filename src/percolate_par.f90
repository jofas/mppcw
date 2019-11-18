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

  cli = read_from_cli()

  call rinit(cli%seed)

  L = cli%matrix_dimension

  call mpi_init()

  call init_cart_comm(comm, L, mpi_comm_world)

  M = comm%M
  N = comm%N

  if (comm%rank == 0) then
    print *, "percolate: params are L = ", L, " rho = ", &
      cli%density_of_filled_cells, " seed = ", cli%seed

    call init_map()
  end if

  call init_chunk()

  call scatter(comm, map, chunk)

  call cluster()

  call gather(comm, map, chunk)

  if (comm%rank == 0) then
    call percolates()

    call pgm_write(cli%pgm_file_path, map, &
      cli%print_n_clusters)
  end if

  call mpi_finalize()

contains

  subroutine init_map()
    !
    ! Initialize the map. Only called by the root process.
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

    print *, "percolate: rho = ", cli%density_of_filled_cells, &
      " actual density = ", float(L ** 2 - free_cell_count) &
      / float(L ** 2)
  end


  subroutine init_chunk()
    allocate(chunk(0:M + 1, 0:N + 1))
    chunk(:, :) = 0
  end


  subroutine cluster()
    !
    ! Performs the clustering operation. First the halos of the chunks
    ! are swapped with its neighbors. Then a cluster step is performed
    ! (every cell of the inner map is set to its biggest neighbor).
    !
    ! Afterwards the sum of every inner cell of the chunk is reduced
    ! over every chunk and broadcasted to every process (mpi_allreduce).
    ! If the gathered sum did not change this cluster step, clustering
    ! is finished.
    !
    integer :: i, sum_, sum_old

    i = 1
    sum_old = 0
    do
      call halo_swap(comm, chunk)

      call cluster_step()

      call mpi_allreduce(sum(chunk(1:M, 1:N)), sum_, 1, &
        mpi_integer, mpi_sum, comm%comm)

      if (sum_ == sum_old) exit

      sum_old = sum_

      call print_map_average(i, sum_)

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


  subroutine print_map_average(i, sum_)
    integer, intent(in) :: i, sum_

    if (mod(i, int(L * cli%print_iter_factor)) == 0 &
      .and. comm%rank == 0) &
    then
      print *, "percolate: average cell value of map on step ", &
        i, " is ", float(sum_) / float(L ** 2)
    end if
  end


  subroutine percolates()
    !
    ! Test if the map percolates (a cluster goes from the
    ! left most column to the right most column of map).
    !
    integer :: i, j
    logical :: does_percolate

    do i = 1, L
      if (map(i, 1) > 0) then
        do j = 1, L
          if (map(i, 1) == map(j, L)) then
            does_percolate = .true.
            exit
          end if
        end do
      end if
    end do

    if (does_percolate) then
      print *, "percolate: cluster DOES percolate"
    else
      print *, "percolate: cluster does NOT percolate"
    end if
  end
end

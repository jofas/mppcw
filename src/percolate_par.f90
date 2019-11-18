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
  use cart
  use scatter_info

  use mpi_f08

  implicit none

  type(CLIResults) :: cli

  type(ScatterInfo) :: scatter_info
  type(CartComm) :: cart_comm
  type(Neighbors) :: neighbors_

  integer :: rank, w_size
  integer :: L

  integer, dimension(:, :), allocatable :: map, chunk

  type(MPI_Datatype) :: chunk_type


  cli = read_from_cli()

  call rinit(cli%seed)

  L = cli%matrix_dimension

  call mpi_init()

  call mpi_comm_rank(mpi_comm_world, rank)
  call mpi_comm_size(mpi_comm_world, w_size)

  if (rank == 0) then
    print *, "percolate: params are L = ", L, " rho = ", &
      cli%density_of_filled_cells, " seed = ", cli%seed
  end if

  call init( &
    map, chunk, chunk_type, scatter_info, cart_comm, &
    rank, w_size, L, cli%density_of_filled_cells &
  )

  call scatter( &
    scatter_info, map, chunk, chunk_type, rank, w_size, &
    cart_comm &
  )

  call cluster( &
    chunk, scatter_info%M, scatter_info%N, L, neighbors_, &
    int(L * cli%print_iter_factor), cart_comm &
  )

  call gather( &
    scatter_info, map, chunk, chunk_type, rank, w_size, &
    cart_comm &
  )

  if (rank == 0) then
    call percolates(map, L)

    call pgm_write(cli%pgm_file_path, map, &
      cli%print_n_clusters)
  end if

  call mpi_finalize()

contains

  subroutine init(map, chunk, chunk_type, scatter_info, cart_, rank, w_size, L, density)
    integer, dimension(:, :), allocatable, intent(out) :: map
    integer, dimension(:, :), allocatable, intent(out) :: chunk
    type(MPI_Datatype), intent(out) :: chunk_type
    type(ScatterInfo), intent(out) :: scatter_info
    type(CartComm), intent(out) :: cart_
    integer, intent(in) :: rank, w_size, L
    real, intent(in) :: density

    if (rank == 0) then
      call init_map(map, L, cli%density_of_filled_cells)
    end if

    call init_cart(cart_, rank, w_size, mpi_comm_world)
    call init_scatter_info(scatter_info, rank, w_size, L, cart_)

    call init_neighbors(neighbors_, scatter_info%M, scatter_info%N, cart_)
    call init_chunk(chunk, scatter_info%M, scatter_info%N)
    call init_chunk_type(chunk_type, scatter_info%M, scatter_info%N)
  end


  subroutine init_map(map, L, d)
    integer, dimension(:, :), allocatable, intent(out) :: map
    integer, intent(in) :: L
    real, intent(in) :: d

    integer :: i, j, free_cell_count

    allocate(map(L, L))
    map(:, :) = 0

    free_cell_count = 0
    do i = 1, L
      do j = 1, L
        if (random_uniform() > d) then
          free_cell_count = free_cell_count + 1
          map(i, j) = free_cell_count
        end if
      end do
    end do

    print *, "percolate: rho = ", d, " actual density = ", &
      float(L ** 2 - free_cell_count) / float(L ** 2)
  end


  subroutine init_chunk(chunk, M, N)
    integer, dimension(:, :), allocatable, intent(out) :: chunk
    integer, intent(in) :: M, N

    allocate(chunk(0:M + 1, 0:N + 1))
    chunk(:, :) = 0
  end


  subroutine init_chunk_type(chunk_type, M, N)
    type(MPI_Datatype), intent(out) ::  chunk_type
    integer, intent(in) :: M, N

    call mpi_type_vector(N, M, M+2, mpi_integer, chunk_type)
    call mpi_type_commit(chunk_type)
  end


  subroutine cluster(chunk, M, N, L, neighbors_, print_modulus, cart_)
    integer, dimension(0:M + 1, 0:N + 1), intent(inout) :: chunk
    integer, intent(in) :: M, N, L, print_modulus
    type(Neighbors), intent(in) :: neighbors_
    type(CartComm), intent(in) :: cart_

    integer :: i, sum_, sum_old

    i = 1
    sum_old = 0
    do
      call halo_swap(neighbors_, chunk, M, N, cart_comm)

      call cluster_step(chunk, M, N)

      call mpi_allreduce(sum(chunk(1:M, 1:N)), sum_, 1, &
        mpi_integer, mpi_sum, cart_comm%comm)

      if (sum_ == sum_old) exit

      sum_old = sum_

      call print_map_average(i, sum_, L, print_modulus)

      i = i + 1
    end do
  end


  subroutine cluster_step(chunk, M, N)
    integer, dimension(0:M + 1, 0:N + 1), intent(inout) :: chunk
    integer, intent(in) :: M, N

    integer :: i, j

    forall(j=1:N, i=1:M, chunk(i, j) /= 0)
      chunk(i, j) = max( chunk(i, j + 1) &
                       , chunk(i, j - 1) &
                       , chunk(i + 1, j) &
                       , chunk(i - 1, j) &
                       , chunk(i, j) )
    end forall
  end


  subroutine print_map_average(i, sum_, L, modulus)
    integer, intent(in) :: i, sum_, L, modulus

    if (mod(i, modulus) == 0 .and. rank == 0 ) then
      print *, "percolate: average cell value of map on step ", &
        i, " is ", float(sum_) / float(L ** 2)
    end if
  end


  subroutine percolates(map, L)
    integer, dimension(L, L), intent(in) :: map
    integer, intent(in) :: L

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

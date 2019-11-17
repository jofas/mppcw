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

  use mpi_f08

  implicit none

  type(CLIResults) :: cli

  integer :: rank, cart_rank, w_size

  integer, dimension(2) :: coords, dims

  integer :: L, M, N

  integer :: clustering_iter

  integer :: i, j
  integer :: chunk_elem_sum, big_elem_sum, big_elem_sum_old

  integer, dimension(:, :), allocatable :: big_map
  integer, dimension(:, :), allocatable :: map_chunk

  integer, dimension(4) :: pointer_to_chunk_in_big_map

  integer, dimension(:, :), allocatable :: mxn_by_rank

  type(MPI_Datatype), dimension(:), allocatable :: send_chunk

  integer :: n_left, n_right, n_upper, n_lower, source

  type(MPI_Status) :: stat
  type(MPI_Comm), save :: comm, comm_cart

  type(MPI_Datatype) :: column, row, inner_big, inner_chunk

  comm = mpi_comm_world

  cli = read_from_cli()

  L = cli%matrix_dimension

  call rinit(cli%seed)

  call mpi_init()

  call mpi_comm_rank(comm, rank)
  call mpi_comm_size(comm, w_size)

  dims(:) = 0
  call mpi_dims_create(w_size, 2, dims)

  ! this thing thinks, like dr henty, that a 2d array can
  ! be thought of as a 2d euclidean space (x-axis as hori-
  ! zontal) instead of a matrix, which is just confusing.
  call mpi_cart_create(comm, 2, dims, &
    [.true., .false.], .false., comm_cart)

  call mpi_cart_coords(comm_cart, rank, 2, coords)
  call mpi_cart_rank(comm_cart, coords, cart_rank)

  call mpi_cart_shift(comm_cart, 1, 1, n_left, n_right)
  call mpi_cart_shift(comm_cart, 0, 1, n_lower, n_upper)

  M = float(L) / float(dims(1))
  N = float(L) / float(dims(2))

  pointer_to_chunk_in_big_map(3) = coords(1) * M + 1
  pointer_to_chunk_in_big_map(4) = coords(2) * N + 1

  if (coords(1) == dims(1) - 1) then
    M = M + mod(L, dims(1))
  end if

  if (coords(2) == dims(2) - 1) then
    N = N + mod(L, dims(2))
  end if

  pointer_to_chunk_in_big_map(1) = M
  pointer_to_chunk_in_big_map(2) = N

  if (rank == 0) then
    allocate(mxn_by_rank(4, 0:w_size-1))
  end if

  call mpi_gather(pointer_to_chunk_in_big_map, 4, mpi_integer, &
    mxn_by_rank, 4, mpi_integer, 0, comm)

  if (rank == 0) then
    allocate(send_chunk(0:w_size-1))

    do i = 0, w_size - 1
      call mpi_type_vector(mxn_by_rank(2, i), mxn_by_rank(1, i), L, mpi_integer, send_chunk(i))
      call mpi_type_commit(send_chunk(i))
    end do
  end if

  call mpi_type_contiguous(M, mpi_integer, column)
  call mpi_type_vector(N, 1, M+2, mpi_integer, row)

  call mpi_type_vector(N, M, M+2, mpi_integer, inner_chunk)

  call mpi_type_commit(column)
  call mpi_type_commit(row)
  call mpi_type_commit(inner_chunk)

  if (rank == 0) then
    call init_big_map(big_map, L, cli%density_of_filled_cells)
  end if

  allocate(map_chunk(0:M + 1, 0:N + 1))

  map_chunk(:, :) = 0

  if (rank == 0) then
    do i = 1, w_size - 1
      call mpi_ssend( &
        big_map(mxn_by_rank(3, i), mxn_by_rank(4, i)), 1, send_chunk(i), &
        i, 0, comm_cart &
      )
    end do

    map_chunk(1:M, 1:N) = big_map(:M, :N)
  else
    call mpi_recv(map_chunk(1, 1), 1, inner_chunk, &
      0, 0, comm_cart, mpi_status_ignore)
  end if

  clustering_iter = 1
  big_elem_sum_old = 0
  do
    ! send right, receive left
    call mpi_sendrecv( &
      map_chunk(1, N), 1, column, n_right, 0, &
      map_chunk(1, 0), 1, column, n_left, 0, &
      comm_cart, mpi_status_ignore &
    )

    ! send left, receive right
    call mpi_sendrecv( &
      map_chunk(1, 1),     1, column, n_left, 0, &
      map_chunk(1, N + 1), 1, column, n_right, 0, &
      comm_cart, mpi_status_ignore &
    )

    ! send upper, receive lower
    call mpi_sendrecv( &
      map_chunk(M, 1), 1, row, n_upper, 0, &
      map_chunk(0, 1), 1, row, n_lower, 0, &
      comm_cart, mpi_status_ignore &
    )

    ! send lower, receive upper
    call mpi_sendrecv( &
      map_chunk(1, 1),     1, row, n_lower, 0, &
      map_chunk(M + 1, 1), 1, row, n_upper, 0, &
      comm_cart, mpi_status_ignore &
    )

    forall(j=1:N, i=1:M, map_chunk(i, j) /= 0)
      map_chunk(i, j) = max( map_chunk(i, j + 1) &
                           , map_chunk(i, j - 1) &
                           , map_chunk(i + 1, j) &
                           , map_chunk(i - 1, j) &
                           , map_chunk(i, j) )
    end forall

    chunk_elem_sum = sum(map_chunk(1:M, 1:N))

    call mpi_allreduce(chunk_elem_sum, big_elem_sum, 1, &
      mpi_integer, mpi_sum, comm)

    if (big_elem_sum == big_elem_sum_old) exit

    big_elem_sum_old = big_elem_sum

    if ( mod( clustering_iter &
            , int(L * cli%print_iter_factor) ) == 0 &
        .and. rank == 0 ) &
    then
      print *, "percolate: average cell value of map on step ", &
        clustering_iter, " is ", float(big_elem_sum) / float(L ** 2)
    end if
    clustering_iter = clustering_iter + 1
  end do

  if (rank == 0) then
    do i = 1, w_size - 1
      call mpi_recv( &
        big_map(mxn_by_rank(3, i), mxn_by_rank(4, i)), 1, send_chunk(i), &
        i, 0, comm_cart, mpi_status_ignore &
      )

    end do

    big_map(:M, :N) = map_chunk(1:M, 1:N)
  else
    call mpi_ssend( &
      map_chunk(1, 1), 1, inner_chunk, &
      0, 0, comm_cart &
    )
  end if

  if (rank == 0) then
    call pgm_write( &
      cli%pgm_file_path, big_map, cli%print_n_clusters &
    )
  end if
  call mpi_finalize()

contains

  subroutine init_big_map(big_map, L, d)
    integer, dimension(:, :), allocatable, intent(out) :: big_map
    integer, intent(in) :: L
    real, intent(in) :: d

    integer :: i, j, free_cell_count

    allocate(big_map(L, L))
    big_map(:, :) = 0

    free_cell_count = 0
    do i = 1, L
      do j = 1, L
        if (random_uniform() > d) then
          free_cell_count = free_cell_count + 1
          big_map(i, j) = free_cell_count
        end if
      end do
    end do
  end
end

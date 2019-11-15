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
  ! It should be noted here, that, while the percolation is
  ! searched from left to right, the output files genera-
  ! ted, display the matix bottom to top. That means, the
  ! last column is displayed as the first row.
  !

  use io
  use uni

  use mpi_f08

  implicit none

  type(CLIResults) :: cli

  !logical :: does_percolate

  integer(kind=mpi_address_kind) :: lb, extent

  integer :: rank, w_size

  integer :: ierr

  integer :: L

  integer :: i_upper, i_lower

  integer :: i, j, k, free_cell_count
  integer :: N
  integer :: changes, changes_sum, max_iter

  integer, dimension(:, :), allocatable :: big_map
  integer, dimension(:, :), allocatable :: map_chunk, old_map_chunk

  integer :: n_upper, n_lower, source
  integer, dimension(1) :: coords

  type(MPI_Status) :: stat
  type(MPI_Comm), save :: comm, comm_cart

  type(MPI_Datatype) :: column

  comm = mpi_comm_world

  cli = read_from_cli()

  L = cli%matrix_dimension

  call rinit(cli%seed)

  call mpi_init()

  call mpi_comm_rank(comm, rank)
  call mpi_comm_size(comm, w_size)

  if (rank == 0) then
    if (mod(L, w_size) /= 0) then
      print *, "Change either mpi world size or matrix dimension"
      call mpi_abort(comm, 0)
      stop
    end if

    call init_big_map(big_map, L, cli%density_of_filled_cells)
  end if

  ! TODO: 2d
  call mpi_cart_create(comm, 1, [w_size], [.false.], &
    .false., comm_cart)

  call mpi_cart_shift(comm_cart, 0,  1, source, n_upper)
  call mpi_cart_shift(comm_cart, 0, -1, source, n_lower)

  N = L / w_size

  allocate(map_chunk(L, 0:N + 1))
  allocate(old_map_chunk(L, 0:N + 1))

  map_chunk(:, :) = 0

  call mpi_scatter(big_map, L * N, mpi_integer, &
    map_chunk(1, 1), L * N, mpi_integer, 0, comm)

  call mpi_type_contiguous(L, mpi_integer, column, ierr)

  print *, ierr

  call mpi_type_get_extent(column, lb, extent)

  print *, lb, extent, L * sizeof(1), sizeof(map_chunk(:, N))

  if (rank == 0) then
    call mpi_ssend(map_chunk(1, N), L, mpi_integer, 1, 0, comm)
    print *, "successfully send"
  elseif (rank == 1) then
    call mpi_recv(map_chunk(1, 0), 1, column, 0, 0, comm, mpi_status_ignore)
    print *, "successfully received"
  end if

  call mpi_finalize()
  stop

  i = 1
  do
    ! send upwards, receive downward

    !call mpi_sendrecv(map_chunk(1, N), 1, column, & !L, mpi_integer, & ! column, &!L, mpi_integer, &
    !  n_upper, 0, &
    !  map_chunk(1, 0), 1, column, & !L, mpi_integer, &!1, column, & !L, mpi_integer,
    !  n_lower, 0, comm_cart, mpi_status_ignore)

    ! send downwards, receive upward
    call mpi_sendrecv(map_chunk(1, 1), L, mpi_integer, &
      n_lower, 0, &
      map_chunk(1, N + 1), L, mpi_integer, &
      n_upper, 0, comm_cart, mpi_status_ignore)

    old_map_chunk(:, :) = map_chunk(:, :)

    do j = 1, N
      do k = 1, L
        if(map_chunk(k, j) /= 0) then
          i_upper = modulo(k, L) + 1
          i_lower = modulo(k - 2, L) + 1

          ! columns with halo
          map_chunk(k, j) = max( map_chunk(k, j + 1) &
                          , map_chunk(k, j - 1) &
                          , map_chunk(i_upper, j) &
                          , map_chunk(i_lower, j) &
                          , map_chunk(k, j) )
        end if
      end do
    end do

    changes = count(map_chunk(:, :) - old_map_chunk(:, :) /= 0)

    if (mod(i, 100) == 0) then

      call mpi_allreduce(changes, changes_sum, 1, &
        mpi_integer, mpi_sum, comm)

      if (changes_sum == 0) exit

      if (rank == 0) print *, rank, i, changes_sum
    end if
    i = i + 1
  end do

  call mpi_gather(map_chunk(1, 1), L * N, mpi_integer, &
    big_map, L * N, mpi_integer, 0, comm)

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

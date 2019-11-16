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

  integer(kind=mpi_address_kind) :: lb, extent

  integer :: rank, cart_rank, w_size

  integer, dimension(2) :: coords, dims

  integer :: L, M, N

  integer :: clustering_iter

  integer :: i, j
  integer :: changes, changes_sum, max_iter

  integer, dimension(:, :), allocatable :: big_map
  integer, dimension(:, :), allocatable :: map_chunk
  integer, dimension(:, :), allocatable :: old_map_chunk

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
  if(rank == 0) print *, "dims", dims

  ! this thing thinks like david that a matrix can be
  ! thought of as a euclidean coordinate system. This is
  ! just rediculous and confusing.
  call mpi_cart_create(comm, 2, dims, &
    [.true., .false.], .false., comm_cart)

  !call mpi_cart_coords(comm_cart, rank, 2, coords)
  !call mpi_cart_rank(comm_cart, coords, cart_rank)

  call mpi_cart_shift(comm_cart, 1, 1, n_left, n_right)
  call mpi_cart_shift(comm_cart, 0, 1, n_lower, n_upper)

  M = float(L) / float(dims(1))
  N = float(L) / float(dims(2))

  print *, "mn", rank, M, N

  call mpi_finalize()
  stop

  call mpi_type_contiguous(M, mpi_integer, column)
  call mpi_type_vector(N, 1, M+2, mpi_integer, row)

  call mpi_type_vector(N, M, L, mpi_integer, inner_big)
  call mpi_type_vector(N, M, M+2, mpi_integer, inner_chunk)

  call mpi_type_commit(column)
  call mpi_type_commit(row)
  call mpi_type_commit(inner_big)
  call mpi_type_commit(inner_chunk)

  if (rank == 0) then

    ! TODO:
    !if (mod(L, w_size) /= 0) then
    !  print *, "Change either mpi world size or matrix dimension"
    !  call mpi_abort(comm, 0)
    !  stop
    !end if

    call init_big_map(big_map, L, cli%density_of_filled_cells)
  end if

  allocate(map_chunk(0:M + 1, 0:N + 1))
  allocate(old_map_chunk(0:M + 1, 0:N + 1))

  !print *, "neighbors", rank, "upper", n_upper, "left", n_left, "lower", n_lower, "right", n_right
  map_chunk(:, :) = 0

  if (rank == 0) then
    do i = 1, w_size - 1
      call mpi_cart_coords(comm_cart, i, 2, coords)
      coords(1) = coords(1) * M + 1
      coords(2) = coords(2) * N + 1

      call mpi_ssend( &
        big_map(coords(1), coords(2)), 1, inner_big, &
        i, 0, comm_cart &
      )
    end do

    map_chunk(1:M, 1:N) = big_map(:M, :N)
  else
    call mpi_recv(map_chunk(1, 1), 1, inner_chunk, &
      0, 0, comm_cart, mpi_status_ignore)
  end if

  !if (rank == 0) then
  !  do i=1,L
  !    print *, big_map(i, :)
  !  end do
  !  print *
  !end if


  clustering_iter = 1
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

    !if (rank == 1) print *, "lower recv", map_chunk(0,1:N)
    !if (rank == 1) print *, "upper recv", map_chunk(M+1,1:N)
    !call mpi_finalize()
    !stop

    old_map_chunk(:, :) = map_chunk(:, :)

    forall(j=1:N, i=1:M, map_chunk(i, j) /= 0)
      map_chunk(i, j) = max( map_chunk(i, j + 1) &
                           , map_chunk(i, j - 1) &
                           , map_chunk(i + 1, j) &
                           , map_chunk(i - 1, j) &
                           , map_chunk(i, j) )
    end forall

    changes = count(map_chunk(:, :) - old_map_chunk(:, :) /= 0)

    call mpi_allreduce(changes, changes_sum, 1, &
      mpi_integer, mpi_sum, comm)

    if (changes_sum == 0) exit

    if (mod(clustering_iter, int(L * 0.5)) == 0 &
      .and. rank == 0) &
    then
      print *, "changes", clustering_iter, changes_sum
    end if
    clustering_iter = clustering_iter + 1
  end do

  !print *, "mc", rank, "c", map_chunk

  if (rank == 0) then
    do i = 1, w_size - 1
      call mpi_cart_coords(comm_cart, i, 2, coords)
      coords(1) = coords(1) * M + 1
      coords(2) = coords(2) * N + 1

      call mpi_recv( &
        big_map(coords(1), coords(2)), 1, inner_big, &
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

  !if (rank == 0) then
  !  do i=1,L
  !    print *, big_map(i, :)
  !  end do
  !  print *
  !end if

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

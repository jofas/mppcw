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

  integer, dimension(:), allocatable :: changes_per_iteration
  integer :: cluster_num
  logical :: does_percolate

  integer :: rank, w_size

  integer :: L

  integer :: i_upper, i_lower

  integer :: i, j, k, free_cell_count
  integer :: m_dim, splitted_m_dim, inner_smap_size
  integer :: changes, changes_sum, max_iter

  integer, dimension(:, :), allocatable :: map_
  integer, dimension(:, :), allocatable :: smap, osmap

  integer :: n_upper, n_lower, source
  integer, dimension(1) :: coords

  type(MPI_Status) :: stat
  type(MPI_Comm), save :: comm, comm_cart
  comm = mpi_comm_world

  cli = read_from_cli()

  L = cli%matrix_dimension

  call rinit(cli%seed)

  call mpi_init()

  call mpi_comm_rank(comm, rank)
  call mpi_comm_size(comm, w_size)

  if (rank == 0) then
    if (mod(cli%matrix_dimension, w_size) /= 0) then
      print *, "Change either mpi world size or matrix dimension"
      call mpi_abort(comm, 0)
      stop
    end if

    allocate(map_(cli%matrix_dimension, cli%matrix_dimension))
    map_(:, :) = 0

    free_cell_count = 0
    do i = 1, cli%matrix_dimension
      do j = 1, cli%matrix_dimension
        if (random_uniform() > cli%density_of_filled_cells) then
          free_cell_count = free_cell_count + 1
          map_(i, j) = free_cell_count
        end if
      end do
    end do
  end if

  call mpi_cart_create(comm, 1, [w_size], [.false.], &
    .false., comm_cart)

  call mpi_cart_shift(comm_cart, 0,  1, source, n_upper)
  call mpi_cart_shift(comm_cart, 0, -1, source, n_lower)

  m_dim = cli%matrix_dimension
  splitted_m_dim = m_dim / w_size
  inner_smap_size = m_dim * splitted_m_dim

  allocate(smap(m_dim, 0:splitted_m_dim + 1))
  smap(:, :) = 0

  call mpi_scatter(map_, inner_smap_size, &
    mpi_integer, smap(:m_dim, 1:splitted_m_dim), &
    inner_smap_size, mpi_integer, 0, comm)

  i = 1
  do
    ! send upwards, receive downward
    call mpi_sendrecv(smap(:, splitted_m_dim), &
      m_dim, mpi_integer, n_upper, 0, &
      smap(:, 0), m_dim, mpi_integer, n_lower, 0, &
      comm_cart, stat)

    ! send downwards, receive upward
    call mpi_sendrecv(smap(:, 1), m_dim, &
      mpi_integer, n_lower, 0, &
      smap(:, splitted_m_dim + 1), m_dim, &
      mpi_integer, n_upper, 0, comm_cart, stat)

    osmap = smap

    do j = 1, splitted_m_dim
      do k = 1, m_dim
      if(smap(k, j) /= 0) then
        if (k == m_dim) then
          i_upper = 1
        else
          i_upper = k + 1
        end if

        if (k == 1) then
          i_lower = m_dim
        else
          i_lower = k - 1
        end if

        ! columns with halo
        smap(k, j) = max( smap(k, j + 1) &
                        , smap(k, j - 1) &
                        , smap(i_upper, j) &
                        , smap(i_lower, j) &
                        , smap(k, j) )
      end if
      end do
    end do

    changes = count(smap(:, :) - osmap(:, :) /= 0)

    if (mod(i, 100) == 0) then

      call mpi_allreduce(changes, changes_sum, 1, &
        mpi_integer, mpi_sum, comm)

      if (changes_sum == 0) exit

      if (rank == 0) print *, rank, i, changes_sum
    end if
    i = i + 1
  end do

  call mpi_gather(smap(:, 1:splitted_m_dim), &
    inner_smap_size, mpi_integer, map_, inner_smap_size, &
    mpi_integer, 0, comm)

  if (rank == 0) then
    call pgm_write( &
      cli%pgm_file_path, map_, cli%print_n_clusters &
    )
  end if
  call mpi_finalize()
end

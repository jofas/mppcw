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

  use map_class
  use sorted_clusters_class
  use color_map_class
  use io
  use uni
  use pgm_out

  use mpi_f08

  implicit none

  type(CLIResults) :: cli
  type(Map) :: m
  type(ColorMap) :: colors
  type(SortedClusters) :: sorted_clusters

  integer, dimension(:), allocatable :: changes_per_iteration
  integer :: cluster_num
  logical :: does_percolate

  integer :: rank, w_size

  integer :: i
  integer :: m_dim, splitted_m_dim, inner_smap_size
  integer :: changes, changes_sum, max_iter

  integer, dimension(:, :), allocatable :: map_
  integer, dimension(:, :), allocatable :: smap, osmap

  type(MPI_Comm), save :: comm
  comm = mpi_comm_world

  cli = read_from_cli()

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

    m = Map(cli%matrix_dimension, cli%density_of_filled_cells)
    map_ = m%map(1:cli%matrix_dimension, 1:cli%matrix_dimension)

    do i = 1, cli%matrix_dimension
      print *, map_(i, :)
    end do
    print *, ""
  end if

  m_dim = cli%matrix_dimension
  splitted_m_dim = m_dim / w_size
  inner_smap_size = m_dim * splitted_m_dim

  allocate(smap(0:m_dim + 1, 0:splitted_m_dim + 1))
  smap(:, :) = 0

  call mpi_scatter(map_, inner_smap_size, &
    mpi_integer, smap(1:m_dim, 1:splitted_m_dim), &
    inner_smap_size, mpi_integer, 0, comm)

  !print *, rank, smap

  max_iter = inner_smap_size * 2 - 1
  changes_sum = 0

  do i = 1, max_iter
    osmap = smap

    where (smap /= 0)
      smap = max( cshift(smap, shift=-1, dim=1) &
                , cshift(smap, shift= 1, dim=1) &
                , cshift(smap, shift=-1, dim=2) &
                , cshift(smap, shift= 1, dim=2) &
                , smap )
    end where

    changes = sum(smap - osmap)

    if (changes == 0) exit

    changes_sum = changes_sum + changes

    if (mod(i, 100) == 0) print *, i, changes_sum
  end do

  !changes_per_iteration = m%build_clusters()

  call mpi_gather(smap(1:m_dim, 1:splitted_m_dim), &
    inner_smap_size, mpi_integer, map_, inner_smap_size, &
    mpi_integer, 0, comm)

  if (rank == 0) then
    do i = 1, m_dim
      print *, map_(i, :)
    end do

    call pgm_write( &
      cli%pgm_file_path, map_, cli%print_n_clusters &
    )
  end if
  call mpi_finalize()

  !does_percolate = m%does_percolate_horizontically(cluster_num)

  !sorted_clusters = SortedClusters(m)

  !call reset_print_n_clusters_if_not_enough_clusters()

  !colors = ColorMap( &
  !  m, sorted_clusters%cluster_ids, cli%print_n_clusters &
  !)

  !call do_output()

contains

  subroutine reset_print_n_clusters_if_not_enough_clusters()
    if (cli%print_n_clusters > sorted_clusters%amount_of_clusters) then
      cli%print_n_clusters = sorted_clusters%amount_of_clusters
    end if
  end


  subroutine do_output()
    if (cli%verbose) call print_infos1()

    call write_data_file(cli%data_file_path, m%inner())

    if (cli%verbose) call print_infos2()

    call write_pgm_file( &
      cli%pgm_file_path, colors%color_map, cli%print_n_clusters &
    )
  end


  subroutine print_infos1()
    call print_params_and_actual_density( &
      cli%density_of_filled_cells, cli%matrix_dimension, &
      cli%seed, m%true_density &
    )

    call print_iterations(changes_per_iteration)

    call print_percolation_status(does_percolate, cluster_num)
  end


  subroutine print_infos2()
    call print_amount_of_clusters_and_size_of_biggest( &
      sorted_clusters%amount_of_clusters, &
      sorted_clusters%cluster_sizes(1)    &
    )

    call print_amount_of_displayed_clusters( &
      cli%print_n_clusters, sorted_clusters%amount_of_clusters &
    )
  end
end

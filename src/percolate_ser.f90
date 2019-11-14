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

  implicit none

  type(CLIResults) :: cli

  integer :: i, j, L, free_cell_count
  integer :: changes, changes_sum, max_iter

  integer, dimension(:, :), allocatable :: map, map_old

  cli = read_from_cli()

  call rinit(cli%seed)

  L = cli%matrix_dimension

  allocate(map(0:L+1, 0:L+1))
  allocate(map_old(0:L+1, 0:L+1))

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

  max_iter = L * 2 - 1
  changes_sum = 0

  do i = 1, max_iter
    map_old(:, :) = map(:, :)

    forall (i=1:L, j=1:L, map(i, j) /= 0)
      map(i, j) = max( map(i - 1, j), map(i + 1, j) &
                     , map(i, j - 1), map(i, j + 1) &
                     , map(i, j) )
    end forall

    changes = sum(map - map_old)

    if (changes == 0) exit

    changes_sum = changes_sum + changes

    if (mod(i, 100) == 0) print *, i, changes_sum
  end do

  call pgm_write( &
    cli%pgm_file_path, map(1:L, 1:L), cli%print_n_clusters &
  )
end

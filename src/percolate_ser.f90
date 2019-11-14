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

  integer :: L

  integer, dimension(:, :), allocatable :: map

  cli = read_from_cli()

  call rinit(cli%seed)

  L = cli%matrix_dimension

  print *, "percolate: params are L = ", L, " rho = ", &
    cli%density_of_filled_cells, " seed = ", cli%seed

  call init_map(map, L)

  call build_clusters(map, L)

  call percolates(map, L)

  call pgm_write( &
    cli%pgm_file_path, map(1:L, 1:L), cli%print_n_clusters &
  )
contains

  subroutine init_map(map, L)
    integer, dimension(:, :), allocatable, intent(out) :: map
    integer, intent(in) :: L

    integer :: i, j, free_cell_count

    allocate(map(L, 0:L+1))

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
      " actual density = ", float(L ** 2 - free_cell_count) / float(L ** 2)
  end


  subroutine build_clusters(map, L)
    integer, dimension(L, 0:L+1), intent(inout) :: map
    integer :: L

    integer, dimension(L, 0:L+1) :: map_old

    integer :: i, changes

    i = 1
    do
      map_old(:, :) = map(:, :)

      where (map /= 0)
        map = max( cshift(map, shift=-1, dim=1) &
                 , cshift(map, shift= 1, dim=1) &
                 , cshift(map, shift=-1, dim=2) &
                 , cshift(map, shift= 1, dim=2) &
                 , map )
      end where

      changes = count(map(:, :) - map_old(:, :) /= 0)

      if (changes == 0) exit

      if (mod(i, 100) == 0) &
        print *, "percolate: number of changes on step ", i, &
          " is ", changes
      i = i + 1
    end do
  end


  subroutine percolates(map, L)
    integer, dimension(L, 0:L), intent(in) :: map
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

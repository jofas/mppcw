module map_class
  use uni

  implicit none

  private

  type, public :: Map
    !
    ! Object implementing the matrix and the perculation
    ! operation based on it.
    !

    integer, dimension(:, :), allocatable :: map
    integer :: inner_size
    real    :: true_density
  contains
    procedure :: build_clusters
    procedure :: does_percolate_horizontically
    procedure :: inner
    procedure, private :: free_random_cells_and_set_true_density
    procedure, private :: set_cell_to_max_neighbor
  end type

  interface Map
    module procedure new
  end interface

contains

  type(Map) function new(inner_size, density) result(self)
    !
    ! Constructor function for a Map instance.
    !
    ! Populates itself with full and empty cells.
    !

    integer, intent(in) :: inner_size
    real,    intent(in) :: density

    allocate(self%map(0:inner_size + 1, 0:inner_size + 1))
    self%map(:, :) = 0

    self%inner_size = inner_size

    call self%free_random_cells_and_set_true_density(density)
  end

  subroutine free_random_cells_and_set_true_density(self, density)
    !
    ! Subroutine initializing self randomly with free cells.
    !

    class(Map), intent(inout) :: self
    real,       intent(in)    :: density

    integer :: free_cell_count, i, j

    free_cell_count = 0
    do i = 1, self%inner_size
      do j = 1, self%inner_size
        if (random_uniform() > density) then
          free_cell_count = free_cell_count + 1
          self%map(i, j) = free_cell_count
        end if
      end do
    end do

    self%true_density = float(free_cell_count) &
                      / float(self%inner_size ** 2)
  end

  function build_clusters(self) result(changes_per_iteration)
    !
    ! Subroutine for building the clusters.
    !
    ! Every empty cell is compared to its neighbors and is
    ! set to the max value of them/itself. This is done so
    ! many times till no more changes happe in an iteration.
    !
    ! Afterwards the cells with the same value build a
    ! cluster.
    !

    class(Map), intent(inout) :: self
    integer, dimension(:), allocatable :: changes_per_iteration

    integer, dimension(:),    allocatable :: changes_per_iteration_temp
    integer, dimension(:, :), allocatable :: old_map

    integer :: changes, i, max_iter

    allocate(old_map(self%inner_size + 2, self%inner_size + 2))

    ! manhattan distance from one end of the matrix diagonal to the
    ! other
    max_iter = self%inner_size * 2 - 1
    allocate(changes_per_iteration_temp(max_iter))

    do i = 1, max_iter
      changes = 0
      old_map = self%map

      call self%set_cell_to_max_neighbor()

      changes = sum(self%map - old_map)
      changes_per_iteration_temp(i) = changes

      if (changes == 0) exit
    end do

    allocate(changes_per_iteration(i))
    changes_per_iteration(:) = changes_per_iteration_temp(:i)

    deallocate(old_map)
    deallocate(changes_per_iteration_temp)
  end

  subroutine set_cell_to_max_neighbor(self)
    class(Map), intent(inout) :: self

    where (self%map /= 0)
      self%map = max( cshift(self%map, shift=-1, dim=1) &
                    , cshift(self%map, shift= 1, dim=1) &
                    , cshift(self%map, shift=-1, dim=2) &
                    , cshift(self%map, shift= 1, dim=2) &
                    , self%map )
    end where
  end

  logical function does_percolate_horizontically(self, cluster_num) &
      result(does_percolate)
    !
    ! Function testing if a cluster percolates from the
    ! left most column to the right most.
    !
    ! It is just checked, whether the first column contains
    ! any value from the last column. If so, the cluster
    ! with the value percolates.
    !

    class(Map), intent(in) :: self
    integer, optional, intent(out) :: cluster_num

    integer :: i, j, MAX_ITER
    integer, dimension(self%inner_size) :: first_col, last_col

    MAX_ITER       = self%inner_size
    cluster_num    = 0
    does_percolate = .false.

    first_col = self%map(:, 1)
    last_col  = self%map(:, MAX_ITER)

    do i = 1, MAX_ITER
      if (first_col(i) > 0) then
        do j = 1, MAX_ITER
          if (first_col(i) == last_col(j)) then
            cluster_num = first_col(i)
            does_percolate = .true.
            exit
          end if
        end do
      end if
    end do
  end

  function inner(self) result(inner_map)
    !
    ! Function returning the inner representation of the
    ! map.
    !
    ! This method is needed, because, for computational
    ! reasons, the matrix is embedded into a halo of full
    ! cells, which is removed invoking this function.
    !

    class(Map), intent(in) :: self
    integer, dimension(self%inner_size, self%inner_size) :: inner_map

    inner_map = self%map(1:self%inner_size, 1:self%inner_size)
  end
end

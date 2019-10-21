module sorted_clusters_class
  use map_class

  implicit none

  private

  type, public :: SortedClusters
    !
    ! Object containing the clusters of a Map instance,
    ! sorted by the cluster's size (amount of cells) in
    ! descending order.
    !

    integer, dimension(:), allocatable :: cluster_ids
    integer, dimension(:), allocatable :: cluster_sizes
    integer :: amount_of_clusters = 0

    integer, private :: size
  contains
    procedure, private :: init_clusters
    procedure, private :: set_amount_of_clusters
    procedure, private :: sort
    procedure, private :: partition
    procedure, private :: swap
  end type

  interface SortedClusters
    module procedure new
  end interface

contains

  type(SortedClusters) function new(m) result(self)
    !
    ! Constructor function for SortedClusters instance.
    !

    type(Map), intent(in) :: m

    call self%init_clusters(m)

    call self%sort(1, self%size)

    call self%set_amount_of_clusters()
  end

  subroutine init_clusters(self, m)
    !
    ! Subroutine initializing the clusters.
    !
    ! Since the map instance has already build the clusters
    ! (every cell is the number of the cell in the cluster
    ! with the highest value assigned in the beginning),
    ! this subroutine now computes the clusters' sizes
    ! (self%cluster_sizes) and associates every cluster
    ! with an id (self%cluster_ids).
    !
    ! The cluster id is needed for mapping the cluster size
    ! back to the cluster, once self%cluster_sizes is
    ! sorted.
    !

    class(SortedClusters), intent(inout) :: self
    type(Map), intent(in) :: m

    integer :: i, j, MAX_ITER

    MAX_ITER  = m%inner_size
    self%size = m%inner_size ** 2

    allocate(self%cluster_ids(self%size))
    allocate(self%cluster_sizes(self%size))

    self%cluster_sizes(:) = 0
    self%cluster_ids(:)   = (/ (i, i = 1, self%size) /)

    forall (i = 1:MAX_ITER, j = 1:MAX_ITER, m%map(i, j) > 0)
      self%cluster_sizes(m%map(i, j)) = &
        self%cluster_sizes(m%map(i, j)) + 1
    end forall
  end

  subroutine set_amount_of_clusters(self)
    !
    ! Subroutine counting the actual amount of clusters.
    !
    ! This is done after the clusters are sorted and
    ! self%cluster_sizes is just iterated till the cluster
    ! size equals zero.
    !

    class(SortedClusters), intent(inout) :: self
    integer :: i

    do i = 1, self%size
      if (self%cluster_sizes(i) > 0) then
        self%amount_of_clusters = i
      else
        exit
      end if
    end do
  end

  recursive subroutine sort(self, begin, end)
    !
    ! Subroutine implementing quicksort.
    !

    class(SortedClusters), intent(inout) :: self
    integer, intent(in) :: begin, end

    integer pivotindex
    pivotindex = begin + (end - begin) / 2

    if (end > begin) then
      call self%partition(begin, end, pivotindex)
      call self%sort(begin, pivotindex - 1)
      call self%sort(pivotindex + 1, end)
    end if
  end

  subroutine partition(self, begin, end, pivotindex)
    !
    ! Subroutine implementing partition operation used by
    ! the quicksort algorithm.
    !

    class(SortedClusters), intent(inout) :: self

    integer, intent(in) :: begin, end
    integer, intent(inout) :: pivotindex

    integer pivot_id, pivot_size, left, right, i, indexpoint

    pivot_id   = self%cluster_ids(pivotindex)
    pivot_size = self%cluster_sizes(pivotindex)

    call self%swap(end, pivotindex)

    left = begin
    right = end - 1
    indexpoint = left

    do i = left, right
      if (self%cluster_sizes(i) > pivot_size .or. &
          (self%cluster_sizes(i) == pivot_size .and. &
            self%cluster_ids(i) > pivot_id)) then

          call self%swap(i, indexpoint)
          indexpoint = indexpoint + 1
      end if
    end do

    call self%swap(indexpoint, end)

    pivotindex = indexpoint
  end

  subroutine swap(self, idx_1, idx_2)
    class(SortedClusters), intent(inout) :: self
    integer, intent(in) :: idx_1, idx_2

    integer :: temp_id, temp_size

    temp_id   = self%cluster_ids(idx_1)
    temp_size = self%cluster_sizes(idx_1)

    self%cluster_ids(idx_1)   = self%cluster_ids(idx_2)
    self%cluster_sizes(idx_1) = self%cluster_sizes(idx_2)

    self%cluster_ids(idx_2)   = temp_id
    self%cluster_sizes(idx_2) = temp_size
  end
end

module color_map_class
  use map_class

  implicit none

  private

  type, public :: ColorMap
    !
    ! Object mapping every inner cell of a Map instance to
    ! color used in the .pgm file.
    !

    integer, dimension(:, :), allocatable :: color_map
  contains
    procedure, private :: build_color_map
  end type

  interface ColorMap
    module procedure new
  end interface

contains

  type(ColorMap) function new( &
    m, cluster_sizes, n_biggest_clusters &
  ) result(self)
    !
    ! Constructor function for a ColorMap instance.
    !
    ! Takes a Map instance and first generates a ranking
    ! for the n_biggest_clusters amount of clusters (based
    ! on the cluster size in descending order).
    !
    ! Based on this ranking the colors are computed, this
    ! time in ascending order. That means the biggest
    ! cluster is always associated with 0, therefore being
    ! displayed as black in the .pgm file. The smaller the
    ! cluster, the lighter the color it is displayed with.
    !
    ! Therfore, every full cell and every cell not con-
    ! tained in one of the n_biggest_clusters is displayed
    ! with the lightest color (n_biggest_clusters).
    !

    type(Map), intent(in) :: m
    integer, intent(in)   :: n_biggest_clusters
    integer, dimension(:), intent(in) :: cluster_sizes

    integer, dimension(m%inner_size ** 2) :: cluster_ranking

    allocate(self%color_map(m%inner_size, m%inner_size))

    call build_cluster_ranking( &
      cluster_sizes, n_biggest_clusters, cluster_ranking &
    )

    call self%build_color_map(m, n_biggest_clusters, cluster_ranking)
  end


  subroutine build_cluster_ranking( &
    cluster_sizes, n_biggest_clusters, cluster_ranking &
  )
    !
    ! Subroutine generating the cluster_ranking.
    !
    ! Assigns the ranks from 0 .. n_biggest_clusters - 1.
    ! The rank is equivalent to the color a cell of the
    ! ColorMap instance will contain, if it is part of the
    ! cluster.
    !

    integer, dimension(:), intent(in)  :: cluster_sizes
    integer,               intent(in)  :: n_biggest_clusters
    integer, dimension(:), intent(out) :: cluster_ranking

    integer :: i

    do i = 1, n_biggest_clusters
      cluster_ranking(cluster_sizes(i)) = i - 1
    end do
  end


  subroutine build_color_map( &
    self, m, n_biggest_clusters, cluster_ranking &
  )
    !
    ! Subroutine generating the color_map.
    !
    ! Maps every inner cell to a color value, based on the
    ! cluster it is part of.
    !

    class(ColorMap), intent(inout) :: self
    type(Map), intent(in) :: m
    integer, intent(in) :: n_biggest_clusters
    integer, dimension(:), intent(in) :: cluster_ranking

    integer :: i, j

    self%color_map(:, :) = n_biggest_clusters

    forall ( &
      j = 1:m%inner_size, &
      i = 1:m%inner_size, &
      cluster_in_bound( &
        i, j, n_biggest_clusters, m, cluster_ranking &
      ) &
    )
      self%color_map(i, j) = cluster_ranking(m%map(i, j))
    end forall
  end


  pure logical function cluster_in_bound( &
    i, j, n_biggest_clusters, m, cluster_ranking &
  ) result(in_bound)
    !
    ! Function checking if a cell is part of one of the
    ! n_biggest_clusters or not.
    !

    integer, intent(in) :: i, j, n_biggest_clusters
    type(Map), intent(in) :: m
    integer, dimension(:), intent(in) :: cluster_ranking

    in_bound = m%map(i, j) > 0 .and. &
      cluster_ranking(m%map(i, j)) < n_biggest_clusters
  end
end

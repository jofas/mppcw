module cart_comm
  !
  ! Module containing a wrapper around a 2d cartesian topology
  ! communication handle (CartComm) and a neighbor type that
  ! holds the information for where and how to halo swap.
  !

  use mpi_f08

  implicit none

  private

  type Neighbors
    integer :: left, right, upper, lower
    type(MPI_Datatype) :: col_type, row_type
  end type

  type ScatterInfo
    integer, dimension(:, :), allocatable :: ptrs
    type(MPI_Datatype), dimension(:), allocatable :: types
    type(MPI_Datatype) :: chunk_type
  end type

  type CartComm
    type(MPI_Comm) :: comm
    integer, dimension(2) :: dims
    integer :: rank, w_size, M, N
    type(Neighbors) :: neighbors_
    type(ScatterInfo) :: scatter_info
  end type

  public :: CartComm
  public :: init_cart_comm
  public :: scatter
  public :: gather
  public :: halo_swap

contains

  subroutine init_cart_comm(self, L, comm)
    !
    ! Initializes the 2d cartesian communicator. The wrapper
    ! contains some more useful information in connection
    ! with the topology, namely the coordinates of the
    ! calling process and the dimensions of the topology.
    !
    type(CartComm), intent(inout) :: self
    integer, intent(in) :: L
    type(MPI_Comm), intent(in) :: comm

    call init_comm(self, comm)
    call init_scatter_info(self, L)
    call init_neighbors(self)
  end


  subroutine init_comm(self, comm)
    type(CartComm), intent(inout) :: self
    type(MPI_Comm), intent(in) :: comm

    ! rank for comm and the created cartesian comm are the
    ! same, because shuffle is set to false in
    ! mpi_cart_create.
    call mpi_comm_rank(comm, self%rank)
    call mpi_comm_size(comm, self%w_size)

    self%dims(:) = 0
    call mpi_dims_create(self%w_size, 2, self%dims)

    ! first axis is horizontal and second vertical (like
    ! coordinate system of a 2d euclidean space rather than
    ! a matrix
    call mpi_cart_create(comm, 2, self%dims, &
      [.true., .false.], .false., self%comm)
  end


  subroutine init_scatter_info(self, L)
    type(CartComm), intent(inout) :: self
    integer, intent(in) :: L

    integer, dimension(self%w_size - 1) :: Ms, Ns
    integer, dimension(4, 0:self%w_size - 1) :: gather_container

    integer :: M, N, M_PTR, N_PTR

    call get_chunk_info(self, M, N, M_PTR, N_PTR, L)

    self%M = M
    self%N = N

    call mpi_gather([M, N, M_PTR, N_PTR], 4, mpi_integer, &
      gather_container, 4, mpi_integer, 0, self%comm)

    if (self%rank == 0) then

      call init_scatter_info_gather_chunk_info(self, Ms, &
        Ns, gather_container)

      call init_scatter_info_types(self, Ms, Ns, L)

    else

      call mpi_type_vector(self%N, self%M, self%M+2, &
        mpi_integer, self%scatter_info%chunk_type)

      call mpi_type_commit(self%scatter_info%chunk_type)
    end if
  end


  subroutine init_neighbors(self)
    !
    ! Initializes the neighbor container (the four neighbor
    ! processes and types for swapping halos.
    !
    type(CartComm), intent(inout) :: self

    call mpi_cart_shift(self%comm, 1, 1, &
      self%neighbors_%left, self%neighbors_%right)

    call mpi_cart_shift(self%comm, 0, 1, &
      self%neighbors_%lower, self%neighbors_%upper)


    call mpi_type_contiguous(self%M, mpi_integer, &
      self%neighbors_%col_type)

    call mpi_type_vector(self%N, 1, self%M+2, mpi_integer, &
      self%neighbors_%row_type)


    call mpi_type_commit(self%neighbors_%col_type)
    call mpi_type_commit(self%neighbors_%row_type)
  end


  subroutine init_scatter_info_gather_chunk_info( &
    self, Ms, Ns, gather_container &
  )
    type(CartComm), intent(inout) :: self
    integer, dimension(self%w_size - 1), intent(out) :: Ms, Ns
    integer, dimension(4, 0:self%w_size - 1), intent(in) &
      :: gather_container

    allocate(self%scatter_info%ptrs(2, self%w_size - 1))
    self%scatter_info%ptrs(:, :) = gather_container(3:, 1:)

    Ms(:) = gather_container(1, 1:)
    Ns(:) = gather_container(2, 1:)
  end


  subroutine init_scatter_info_types(self, Ms, Ns, L)
    type(CartComm), intent(inout) :: self
    integer, dimension(self%w_size - 1), intent(in) :: Ms, Ns
    integer, intent(in) :: L

    integer :: i

    allocate(self%scatter_info%types(self%w_size - 1))

    do i = 1, self%w_size - 1
      call mpi_type_vector( &
        Ns(i), Ms(i), L, mpi_integer, self%scatter_info%types(i) &
      )

      call mpi_type_commit(self%scatter_info%types(i))
    end do
  end


  subroutine get_chunk_info(self, M, N, M_PTR, N_PTR, L)
    type(CartComm), intent(in)  :: self
    integer, intent(out) :: M, N, M_PTR, N_PTR
    integer, intent(in) :: L

    integer, dimension(2) :: coords

    call mpi_cart_coords(self%comm, self%rank, 2, coords)

    M = float(L) / float(self%dims(1))
    N = float(L) / float(self%dims(2))

    M_PTR = coords(1) * M + 1
    N_PTR = coords(2) * N + 1

    if (coords(1) == self%dims(1) - 1) then
      M = M + mod(L, self%dims(1))
    end if

    if (coords(2) == self%dims(2) - 1) then
      N = N + mod(L, self%dims(2))
    end if
  end


  subroutine gather(self, map, chunk)
    type(CartComm), intent(in) :: self
    integer, dimension(:, :), intent(inout) :: map
    integer, dimension(0:self%M + 1, 0:self%N + 1), &
      intent(in) :: chunk

    integer i

    if (self%rank == 0) then
      call recv_map(self, map)

      map(:self%M, :self%N) = chunk(1:self%M, 1:self%N)
    else
      call mpi_ssend( &
        chunk(1, 1), 1, self%scatter_info%chunk_type, &
        0, 0, self%comm &
      )
    end if
  end


  subroutine ssend_map(self, map)
    type(CartComm), intent(in) :: self
    integer, dimension(:, :), intent(in) :: map

    integer :: i

    do i = 1, self%w_size - 1
      call mpi_ssend( &
        map( self%scatter_info%ptrs(1, i)    &
           , self%scatter_info%ptrs(2, i) ), &
        1, self%scatter_info%types(i), &
        i, 0, self%comm &
      )
    end do
  end


  subroutine recv_map(self, map)
    type(CartComm), intent(in) :: self
    integer, dimension(:, :), intent(inout) :: map

    integer :: i

    do i = 1, self%w_size - 1
      call mpi_recv( &
        map( self%scatter_info%ptrs(1, i)    &
           , self%scatter_info%ptrs(2, i) ), &
        1, self%scatter_info%types(i), &
        i, 0, self%comm, mpi_status_ignore &
      )
    end do
  end


  subroutine halo_swap(self, chunk)
    type(CartComm), intent(in) :: self
    integer, dimension(0:self%M + 1, 0:self%N + 1), &
      intent(inout) :: chunk

    ! send right, receive left
    call mpi_sendrecv( &
      chunk(1, self%N), 1, self%neighbors_%col_type, &
      self%neighbors_%right, 0, &
      chunk(1, 0), 1, self%neighbors_%col_type, &
      self%neighbors_%left, 0, &
      self%comm, mpi_status_ignore &
    )

    ! send left, receive right
    call mpi_sendrecv( &
      chunk(1, 1), 1, self%neighbors_%col_type, &
      self%neighbors_%left, 0, &
      chunk(1, self%N + 1), 1, self%neighbors_%col_type, &
      self%neighbors_%right, 0, &
      self%comm, mpi_status_ignore &
    )

    ! send upper, receive lower
    call mpi_sendrecv( &
      chunk(self%M, 1), 1, self%neighbors_%row_type, &
      self%neighbors_%upper, 0, &
      chunk(0, 1), 1, self%neighbors_%row_type, &
      self%neighbors_%lower, 0, &
      self%comm, mpi_status_ignore &
    )

    ! send lower, receive upper
    call mpi_sendrecv( &
      chunk(1, 1),     1, self%neighbors_%row_type, &
      self%neighbors_%lower, 0, &
      chunk(self%M + 1, 1), 1, self%neighbors_%row_type, &
      self%neighbors_%upper, 0, &
      self%comm, mpi_status_ignore &
    )
  end


  subroutine scatter(self, map, chunk)
    type(CartComm), intent(in) :: self
    integer, dimension(:, :), intent(in) :: map
    integer, dimension(0:self%M + 1, 0:self%N + 1), &
      intent(inout) :: chunk

    if (self%rank == 0) then
      call ssend_map(self, map)

      chunk(1:self%M, 1:self%N) =  map(:self%M, :self%N)
    else
      call mpi_recv( &
        chunk(1, 1), 1, self%scatter_info%chunk_type, &
        0, 0, self%comm, mpi_status_ignore &
      )
    end if
  end
end
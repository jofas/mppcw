module cart_comm
  !
  ! Module containing a wrapper around a 2d cartesian topology
  ! communication handle (CartComm).
  !
  ! CartComm does all the message passing of percolate:
  !
  !   * scattering and gathering the map to/from its chunks
  !
  !   * halo swapping between neighboring processes
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
    ! Initializes the actual communicator, the information
    ! on how to scatter/gather the map to/from the chunks
    ! and the neighbors of the calling process used for
    ! halo swapping.
    !
    type(CartComm), intent(inout) :: self
    integer, intent(in) :: L
    type(MPI_Comm), intent(in) :: comm

    call init_comm(self, comm)
    call init_scatter_info(self, L)
    call init_neighbors(self)
  end


  subroutine init_comm(self, comm)
    !
    ! Builds the 2d cartesian communicator from comm and
    ! saves some more information about the environment
    ! (how the dimensions are split among processes, the
    ! amount of processes (w_size) and the rank of the
    ! calling process, which is needed later for identifi-
    ! cation of its neighbors and which chunk is passed to
    ! the process).
    !
    type(CartComm), intent(inout) :: self
    type(MPI_Comm), intent(in) :: comm

    ! rank for comm and the created cartesian comm are the
    ! same, because shuffle is set to false in
    ! mpi_cart_create.
    call mpi_comm_rank(comm, self%rank)
    call mpi_comm_size(comm, self%w_size)

    self%dims(:) = 0
    call mpi_dims_create(self%w_size, 2, self%dims)

    call mpi_cart_create(comm, 2, self%dims, &
      [.true., .false.], .false., self%comm)
  end


  subroutine init_scatter_info(self, L)
    !
    ! ScatterInfo contains information on how to scatter
    ! and gather the map to/from its chunks, which are
    ! distributed among the processes. It is a rather in-
    ! tricate data structure, because it differs depending
    ! on the calling processes. If the root process calls
    ! this routine, ScatterInfo contains information on
    ! how to send/receive data from/to the map, while all
    ! other processes contain information on how to send/
    ! receive data from/to chunks.
    !
    ! The root process ScatterInfo instance is more com-
    ! plex. For each process 1...|amount of processes| it
    ! contains a pointer to the first element of the chunk
    ! and a datatype (strided vector) for how to send the
    ! chunk. The other processes just have a datatype for
    ! receiving the chunk in its inner container (without
    ! the halos).
    !
    ! The chunks may differ between processes (when the
    ! amount of processes and L are not symmetric). This
    ! is the reason every process sends information about
    ! its chunk to the root process which then builds its
    ! ScatterInfo instance.
    !
    type(CartComm), intent(inout) :: self
    integer, intent(in) :: L

    integer, dimension(4, 0:self%w_size - 1) :: gather_container

    call init_chunk_info_and_gather(self, gather_container, L)

    if (self%rank == 0) then

      call init_scatter_info_root(self, gather_container, L)

    else

      call init_scatter_info_not_root(self)

    end if
  end


  subroutine init_neighbors(self)
    !
    ! Initializes the neighbor container (the four neighbor
    ! processes and types for swapping halos).
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


  subroutine init_chunk_info_and_gather(self, gather_container, L)
    !
    ! The chunks are distributed as equal as possible.
    ! If L and the amount of processes is not symmetric,
    ! the processes at the edges take the rest of
    ! L mod |amount of processes per dimension| as additi-
    ! onal chunk size.
    !
    ! This routine returns a gather container (only on the
    ! root process). For each process, it contains a row
    ! with the information:
    !
    !   1. M
    !
    !   2. N
    !
    !   3. M_PTR
    !
    !   4. N_PTR
    !
    ! The gather container is gathered to the root process
    ! which uses this information to build the pointers to
    ! the chunks. M_PTR, N_PTR point to the first element
    ! of the chunk and M, N are used to define a strided
    ! vector type for sending the chunk.
    !
    type(CartComm), intent(inout) :: self
    integer, dimension(4, 0:self%w_size - 1), intent(out) &
      :: gather_container
    integer, intent(in) :: L

    integer, dimension(2) :: coords
    integer :: M_PTR, N_PTR

    call mpi_cart_coords(self%comm, self%rank, 2, coords)

    self%M = float(L) / float(self%dims(1))
    self%N = float(L) / float(self%dims(2))

    M_PTR = coords(1) * self%M + 1
    N_PTR = coords(2) * self%N + 1

    if (coords(1) == self%dims(1) - 1) then
      self%M = self%M + mod(L, self%dims(1))
    end if

    if (coords(2) == self%dims(2) - 1) then
      self%N = self%N + mod(L, self%dims(2))
    end if

    call mpi_gather([self%M, self%N, M_PTR, N_PTR], &
      4, mpi_integer, gather_container, 4, mpi_integer, &
      0, self%comm)
  end


  subroutine init_scatter_info_root(self, gather_container, L)
    !
    ! Gather the information about the chunks and build the
    ! pointers to the first element of a chunk. The cor-
    ! responding types are initialized in
    ! init_scatter_info_types.
    !
    ! Only called from the root process.
    !
    type(CartComm), intent(inout) :: self
    integer, dimension(4, 0:self%w_size - 1), intent(in) &
      :: gather_container
    integer, intent(in) :: L

    allocate(self%scatter_info%ptrs(2, self%w_size - 1))

    self%scatter_info%ptrs(:, :) = gather_container(3:, 1:)

    call init_scatter_info_root_types(self, gather_container, L)
  end


  subroutine init_scatter_info_not_root(self)
    !
    ! Fills the scatter info with a strided vector type
    ! used for sending/receiving a chunk to/from the map.
    !
    ! Only called from non root processes.
    !
    type(CartComm), intent(inout) :: self

    call mpi_type_vector(self%N, self%M, self%M+2, &
      mpi_integer, self%scatter_info%chunk_type)

    call mpi_type_commit(self%scatter_info%chunk_type)
  end


  subroutine init_scatter_info_root_types(self, gather_container, L)
    !
    ! These are the strided vector types for scattering/
    ! gathering the chunks to/from the map.
    !
    ! Only called from the root process.
    !
    type(CartComm), intent(inout) :: self
    integer, dimension(4, 0:self%w_size - 1), intent(in) &
      :: gather_container
    integer, intent(in) :: L

    integer, dimension(self%w_size - 1) :: Ms, Ns
    integer :: i

    Ms(:) = gather_container(1, 1:)
    Ns(:) = gather_container(2, 1:)

    allocate(self%scatter_info%types(self%w_size - 1))

    do i = 1, self%w_size - 1
      call mpi_type_vector(Ns(i), Ms(i), L, mpi_integer, &
        self%scatter_info%types(i))

      call mpi_type_commit(self%scatter_info%types(i))
    end do
  end


  subroutine scatter(self, map, chunk)
    !
    ! Subroutine scattering map to the chunks. The root
    ! process synchronously sends the chunk to the cor-
    ! responding process.
    !
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


  subroutine gather(self, map, chunk)
    !
    ! Subroutine gathering map from the chunks. Every
    ! process synchronously sends its chunk to the root
    ! process.
    !
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
    !
    ! Send a chunk defined by its pointer and type.
    !
    ! Only called by the root process.
    !
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
    !
    ! Receive a chunk defined by its pointer and type.
    !
    ! Only called by the root process.
    !
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
    !
    ! Swap chunk halos between neighboring processes.
    !
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
end

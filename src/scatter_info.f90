module scatter_info

  use cart

  use mpi_f08

  implicit none

  private

  type ScatterInfo
    integer, dimension(:, :), allocatable :: ptrs
    type(MPI_Datatype), dimension(:), allocatable :: types
    integer :: M, N
  end type

  public :: ScatterInfo
  public :: init_scatter_info
  public :: scatter
  public :: gather

contains

  subroutine init_scatter_info(self, rank, w_size, L, cart_)
    type(ScatterInfo), intent(inout) :: self
    integer, intent(in) :: rank, w_size, L
    type(CartComm), intent(in) :: cart_
    integer, dimension(w_size - 1) :: Ms, Ns
    integer, dimension(4, 0:w_size - 1) :: gather_container
    integer :: M, N, M_PTR, N_PTR

    call get_chunk_info(M, N, M_PTR, N_PTR, L, cart_)

    call mpi_gather([M, N, M_PTR, N_PTR], 4, mpi_integer, &
      gather_container, 4, mpi_integer, 0, cart_%comm)

    if (rank == 0) then
      allocate(self%ptrs(2, w_size - 1))
      self%ptrs(:, :) = gather_container(3:, 1:)

      Ms(:) = gather_container(1, 1:)
      Ns(:) = gather_container(2, 1:)

      call scatter_info_init_types(self, Ms, Ns, w_size, L)
    end if

    self%M = M
    self%N = N
  end


  subroutine scatter(self, map, chunk, recv_type, rank, w_size, cart_)
    type(ScatterInfo), intent(in) :: self
    integer, dimension(:, :), intent(in) :: map
    integer, dimension(0:self%M + 1, 0:self%N + 1), intent(inout) :: chunk
    type(MPI_Datatype), intent(in) :: recv_type
    integer, intent(in) :: rank, w_size
    type(CartComm), intent(in) :: cart_

    if (rank == 0) then
      call ssend_map(self, map, w_size, cart_)
      chunk(1:self%M, 1:self%N) = map(:self%M, :self%N)
    else
      call mpi_recv(chunk(1, 1), 1, recv_type, &
        0, 0, cart_%comm, mpi_status_ignore)
    end if
  end


  subroutine gather(self, map, chunk, send_type, rank, w_size, cart_)
    type(ScatterInfo), intent(in) :: self
    integer, dimension(:, :), intent(inout) :: map
    integer, dimension(0:self%M + 1, 0:self%N + 1), intent(in) :: chunk
    type(MPI_Datatype), intent(in) :: send_type
    integer, intent(in) :: rank, w_size
    type(CartComm), intent(in) :: cart_

    integer i

    if (rank == 0) then
      call recv_map(self, map, w_size, cart_)
      map(:self%M, :self%N) = chunk(1:self%M, 1:self%N)
    else
      call mpi_ssend(chunk(1, 1), 1, send_type, &
        0, 0, cart_%comm)
    end if
  end


  subroutine ssend_map(self, map, w_size, cart_)
    type(ScatterInfo), intent(in) :: self
    integer, dimension(:, :), intent(in) :: map
    integer, intent(in) :: w_size
    type(CartComm), intent(in) :: cart_

    integer :: i

    do i = 1, w_size - 1
      call mpi_ssend( &
        map(self%ptrs(1, i), self%ptrs(2, i)), &
        1, self%types(i), i, 0, cart_%comm &
      )
    end do
  end


  subroutine recv_map(self, map, w_size, cart_)
    type(ScatterInfo), intent(in) :: self
    integer, dimension(:, :), intent(inout) :: map
    integer, intent(in) :: w_size
    type(CartComm), intent(in) :: cart_

    integer :: i
    do i = 1, w_size - 1
      call mpi_recv( &
        map(self%ptrs(1, i), self%ptrs(2, i)), &
        1, self%types(i), i, 0, cart_%comm, mpi_status_ignore &
      )
    end do
  end

  subroutine get_chunk_info(M, N, M_PTR, N_PTR, L, cart_)
    integer, intent(out) :: M, N, M_PTR, N_PTR
    integer, intent(in) :: L
    type(CartComm), intent(in)  :: cart_

    M = float(L) / float(cart_%dims(1))
    N = float(L) / float(cart_%dims(2))

    M_PTR = cart_%coords(1) * M + 1
    N_PTR = cart_%coords(2) * N + 1

    if (cart_%coords(1) == cart_%dims(1) - 1) then
      M = M + mod(L, cart_%dims(1))
    end if

    if (cart_%coords(2) == cart_%dims(2) - 1) then
      N = N + mod(L, cart_%dims(2))
    end if
  end


  subroutine scatter_info_init_types(self, Ms, Ns, w_size, L)
    type(ScatterInfo), intent(inout) :: self
    integer, dimension(w_size - 1), intent(in) :: Ms, Ns
    integer, intent(in) :: w_size, L

    integer :: i

    allocate(self%types(w_size - 1))

    do i = 1, w_size - 1
      call mpi_type_vector( &
        Ns(i), Ms(i), L, mpi_integer, self%types(i) &
      )

      call mpi_type_commit(self%types(i))
    end do
  end
end

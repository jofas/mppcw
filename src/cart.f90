module cart

  use mpi_f08

  implicit none

  public

  type CartComm
    type(MPI_Comm) :: comm
    integer, dimension(2) :: dims
    integer, dimension(2) :: coords
  end type

  type Neighbors
    integer :: left, right, upper, lower
    type(MPI_Datatype) :: col_type, row_type
  end type

  private :: get_neighbors

contains

  subroutine init_cart(self, rank, w_size, comm)
    type(CartComm), intent(inout) :: self
    integer, intent(in) :: rank, w_size
    type(MPI_Comm), intent(in) :: comm

    self%dims(:) = 0
    call mpi_dims_create(w_size, 2, self%dims)

    ! this thing thinks, like dr henty, that a 2d array can
    ! be thought of as a 2d euclidean space (x-axis as ho-
    ! rizontal) instead of a matrix, which is just confu-
    ! sing.
    call mpi_cart_create(comm, 2, self%dims, &
      [.true., .false.], .false., self%comm)

    call mpi_cart_coords(self%comm, rank, 2, self%coords)
  end


  subroutine get_neighbors(self, n)
    type(CartComm), intent(in) :: self
    type(Neighbors), intent(out) :: n

    call mpi_cart_shift(self%comm, 1, 1, n%left, n%right)
    call mpi_cart_shift(self%comm, 0, 1, n%lower, n%upper)
  end


  subroutine init_neighbors(self, M, N, cart_)
    type(Neighbors), intent(out) :: self
    integer, intent(in) :: M, N
    type(CartComm), intent(in) :: cart_

    call get_neighbors(cart_, self)

    call mpi_type_contiguous(M, mpi_integer, self%col_type)
    call mpi_type_vector(N, 1, M+2, mpi_integer, self%row_type)

    call mpi_type_commit(self%col_type)
    call mpi_type_commit(self%row_type)
  end


  subroutine halo_swap(self, chunk, M, N, cart_)
    type(Neighbors), intent(in) :: self
    integer, dimension(0:M + 1, 0:N + 1), intent(inout) :: chunk
    integer, intent(in) :: M, N
    type(CartComm), intent(in) :: cart_


    ! send right, receive left
    call mpi_sendrecv( &
      chunk(1, N), 1, self%col_type, self%right, 0, &
      chunk(1, 0), 1, self%col_type, self%left, 0, &
      cart_%comm, mpi_status_ignore &
    )

    ! send left, receive right
    call mpi_sendrecv( &
      chunk(1, 1),     1, self%col_type, self%left, 0, &
      chunk(1, N + 1), 1, self%col_type, self%right, 0, &
      cart_%comm, mpi_status_ignore &
    )

    ! send upper, receive lower
    call mpi_sendrecv( &
      chunk(M, 1), 1, self%row_type, self%upper, 0, &
      chunk(0, 1), 1, self%row_type, self%lower, 0, &
      cart_%comm, mpi_status_ignore &
    )

    ! send lower, receive upper
    call mpi_sendrecv( &
      chunk(1, 1),     1, self%row_type, self%lower, 0, &
      chunk(M + 1, 1), 1, self%row_type, self%upper, 0, &
      cart_%comm, mpi_status_ignore &
    )
  end
end

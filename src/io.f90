module io
  !
  ! Module for printing (stdout), writing (files)
  ! and the command line interface.
  !

  use cli_info

  implicit none

  private

  integer                  :: IOUNIT  = 12
  integer, parameter       :: STR_LEN = 32
  character(len = STR_LEN) :: FMT_STRING
  character(len = STR_LEN) :: CLI_OPT
  character(len = STR_LEN) :: CLI_IN

  type, public :: CLIResults
    !
    ! Object containing the values for the cli options.
    !

    integer :: matrix_dimension
    integer :: seed
    integer :: print_n_clusters
    real    :: density_of_filled_cells
    logical :: verbose
    character(len = STR_LEN) :: data_file_path
    character(len = STR_LEN) :: pgm_file_path
  contains
    procedure, private :: parse_argument
    procedure, private :: set_default
  end type

  public read_from_cli
  public write_data_file
  public write_pgm_file
  public print_params_and_actual_density
  public print_iterations
  public print_percolation_status
  public print_amount_of_clusters_and_size_of_biggest
  public print_amount_of_displayed_clusters

contains

  type(CLIResults) function read_from_cli() result(self)
    !
    ! Function which reads the command line arguments.
    !
    ! Use percolate -h for more infos.
    !
    ! If parsing the command line arguments fails, the
    ! execution is stopped and an error message is
    ! displayed.
    !

    integer :: i

    logical :: print_n_clusters_set = .false.

    call self%set_default()

    i = 1
    do
      call get_command_argument(i, CLI_OPT)

      call self%parse_argument(i, print_n_clusters_set)

      ! >=, to catch no commands
      if (i >= command_argument_count()) exit
      i = i + 1
    end do
  end


  subroutine parse_argument(self, i, print_n_clusters_set)
    !
    ! subroutine which parses the current command (an a
    ! corresponding value, if the command needs one).
    !
    ! If parsing the command line argument fails, the
    ! execution is stopped and an error message is dis-
    ! played.
    !

    class(CLIResults), intent(inout) :: self
    integer, intent(inout) :: i
    logical, intent(inout) :: print_n_clusters_set

    select case(CLI_OPT)

      case("--length", "-l")
        call read_command_value(i)

        self%matrix_dimension = parse_command_to_int()

        ! if print_n_clusters is not set by the command
        ! line interface, it is overridden and set to
        ! "all".
        if (.not. print_n_clusters_set) then
          self%print_n_clusters = self%matrix_dimension ** 2
        end if

      case("--density", "-d")
        call read_command_value(i)

        self%density_of_filled_cells = parse_command_to_real()

      case("--seed", "-s")
        call read_command_value(i)

        self%seed = parse_command_to_int()

      case("--print_n_clusters", "-p")
        call read_command_value(i)

        self%print_n_clusters = parse_command_to_int()

        print_n_clusters_set = .true.

      case("--verbose", "-v")
        self%verbose = .true.

      case("--data_file_path")
        call read_command_value(i)
        self%data_file_path = trim(cli_in)

      case("--pgm_file_path")
        call read_command_value(i)
        self%pgm_file_path = trim(cli_in)

      case("--help", "-h")
        call write_help_msg()
        stop

      case("--version")
        call write_version()
        stop

      case("")
        return

      case default
        write (*, *) "Command line arguments are wrong. &
                     See -h, --help for further information."
        stop
    end select
  end


  subroutine set_default(self)
    class(CLIResults), intent(out) :: self

    self%matrix_dimension        = 20
    self%seed                    = 1564
    self%print_n_clusters        = 20 ** 2
    self%density_of_filled_cells = .4
    self%verbose                 = .false.
    self%data_file_path          = "map.dat"
    self%pgm_file_path           = "map.pgm"
  end


  subroutine read_command_value(i)
    integer, intent(inout) :: i

    i = i + 1
    call get_command_argument(i, CLI_IN)
  end


  integer function parse_command_to_int() result(res)
    !
    ! Function parsing CLI_IN to an integer.
    !
    ! If the parsing fails, the program is stopped and an
    ! error message is displayed.
    !

    integer :: stat

    read (CLI_IN, *, iostat = stat) res

    if (stat /= 0) then
      write (*, *) "Could not parse: ", trim(CLI_OPT), &
                   ": ", CLI_IN
      stop
    end if
  end


  real function parse_command_to_real() result(res)
    !
    ! Function parsing CLI_IN to a real.
    !
    ! If the parsing fails, the program is stopped and an
    ! error message is displayed.
    !

    integer :: stat

    read (CLI_IN, *, iostat = stat) res

    if (stat /= 0) then
      write (*, *) "Could not parse: ", trim(CLI_OPT), &
                   ": ", CLI_IN
      stop
    end if
  end


  subroutine write_data_file(path, inner_map)
    !
    ! Subroutine writing the inner_map to file (path).
    !
    ! It writes the matrix column wise, beginning with the
    ! last column.
    !

    character(len = STR_LEN), intent(in) :: path
    integer, dimension(:, :), intent(in) :: inner_map

    integer :: inner_map_size, j

    inner_map_size = size(inner_map, dim=1)

    write (*, *) 'Opening file ', path
    open (unit=IOUNIT, file=path)
    write (*, *) 'Writing data ...'

    write (FMT_STRING, fmt='(''('', i3, ''(1x, I4))'')') &
      inner_map_size

    do j = inner_map_size, 1, -1
      write (IOUNIT, fmt=FMT_STRING) inner_map(:, j)
    end do

    write (*, *) '...done'
    close (IOUNIT)
    write (*, *) 'File closed'
  end


  subroutine write_pgm_file(path, colors, amount_of_clusters)
    !
    ! Subroutine writing the colors matrix to a .pgm file
    ! (path).
    !
    ! Like write_data_file, the colors matrix is written
    ! column wise, beginning with the last column.
    !

    character(len = STR_LEN), intent(in) :: path
    integer, dimension(:, :), intent(in) :: colors
    integer, intent(in) :: amount_of_clusters

    integer :: color_size, j

    color_size = size(colors, dim=1)

    write (*, *) 'Opening file ', path
    open (unit=IOUNIT, file=path)
    write (*, *) 'Writing data ...'

    call write_pgm_header(amount_of_clusters, color_size)

    do j = color_size, 1, -1
      write (IOUNIT, fmt=FMT_STRING) colors(:, j)
    end do

    write (*, *) '...done'
    close (IOUNIT)
    write (*, *) 'File closed'
  end


  subroutine write_pgm_header(amount_of_clusters, color_size)
    !
    ! Subroutine writing some header information to the
    ! .pgm file (handled by IOUNIT).
    !

    integer, intent(in) :: amount_of_clusters, color_size

    write (FMT_STRING, fmt='(''('', i3, ''(1x, I4))'')') &
      color_size

    write (IOUNIT, fmt='(''P2'')')

    write (IOUNIT, *) color_size, color_size

    if (amount_of_clusters > 0) then
      write (IOUNIT, *) amount_of_clusters
    else
      write (IOUNIT, *) 1
    end if
  end


  subroutine print_params_and_actual_density( &
    rho, L, seed, true_density &
  )
    real, intent(in)    :: rho, true_density
    integer, intent(in) :: L, seed

    write (*, *) 'Parameters are rho=', rho, ', L=', L, &
                 ', seed=', seed
    write (*, *) 'rho = ', rho, ', actual density = ', true_density
  end


  subroutine print_iterations(changes_per_iteration)
    integer, dimension(:), intent(in) :: changes_per_iteration

    integer :: i

    do i = 1, size(changes_per_iteration)
      write (*, *) 'Number of changes on loop ', i, ' is ', &
                   changes_per_iteration(i)
    end do
  end


  subroutine print_percolation_status(does_percolate, cluster_num)
    logical, intent(in) :: does_percolate
    integer, intent(in) :: cluster_num

    if (does_percolate) then
      write (*, *) 'Cluster DOES percolate. Cluster number: ', &
                   cluster_num
    else
      write (*, *) 'Cluster DOES NOT percolate'
    end if
  end


  subroutine print_amount_of_clusters_and_size_of_biggest( &
    amount, biggest_size &
  )
    integer, intent(in) :: amount, biggest_size

    write (*, *) 'Map has ', amount, &
      ' clusters, maximum cluster size is ', biggest_size
  end


  subroutine print_amount_of_displayed_clusters(print_n, actual_n)
    integer, intent(in) :: print_n, actual_n

    if (print_n == 1) then
      write (*, *) 'Displaying the largest cluster'
    else if (print_n == actual_n) then
      write (*, *) 'Displaying all clusters'
    else
      write (*, *) 'Displaying largest ', print_n, ' clusters'
    end if
  end
end

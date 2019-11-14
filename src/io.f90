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
    character(len = STR_LEN) :: pgm_file_path
  contains
    procedure, private :: parse_argument
    procedure, private :: set_default
  end type

  public read_from_cli
  public pgm_write

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


  subroutine pgm_write(percfile, map, ncluster)
    !
    ! Function to write a percolation map in greyscale Portable Grey
    ! Map (PGM) format. The largest "ncluster" clusters are identified
    ! and shown as shades of grey against a black background, with the
    ! largest cluster shown in white.
    !

    character (len = *), intent(in)     :: percfile
    integer, dimension(:,:), intent(in) :: map
    integer, intent(in)                 :: ncluster

    integer, parameter :: maxncluster = 9  ! Must be a single digit
    integer, parameter :: pixperline = 32  ! PGM limit 70 chars per line

    integer, dimension(pixperline)     :: pgmline
    integer, dimension(maxncluster)    :: foundcluster
    integer, allocatable, dimension(:) :: clustersize

    integer, parameter :: fmtlen = 64
    character (len = fmtlen)  :: fmtstring
    integer, parameter :: iounit = 12

    integer :: m, n

    integer :: i, j, npix, colour
    integer :: clusterid, icluster, lncluster, prevcluster, maxcluster

    m = size(map, 1)
    n = size(map, 2)

    allocate(clustersize(m*n))

    lncluster = ncluster

    if (lncluster > maxncluster) then

      write(*,*) "percwrite: WARNING ncluster too large, resetting to ", &
           maxncluster

      lncluster = maxncluster

   end if

   if (lncluster > 1) then

      write(*,*) "percwrite: visualising the largest ", lncluster, "clusters"

   else

      write(*,*) "percwrite: only visualising the largest cluster"

   end if

   ! Count up the size of each cluster

   clustersize(:) = 0

   do i = 1, m
      do j = 1, n

         clusterid = map(i,j)

         if (clusterid > 0) then

            clustersize(clusterid) = clustersize(clusterid) + 1

         end if

      end do
   end do

   !  Find the size of the "lncluster" largest clusters (by brute force!)

   prevcluster = m*n+1  ! Larger than the largest possible cluster

   do icluster = 1, lncluster

      maxcluster = 0

      do i = 1, m*n

         if (clustersize(i) > maxcluster .and. &
             clustersize(i) < prevcluster        ) then

            maxcluster = clustersize(i)

         end if

      end do

      foundcluster(icluster) = maxcluster
      prevcluster = maxcluster

   end do

   if (lncluster > 1) then
      write(*,*) "percwrite: cluster sizes are ", foundcluster(1:lncluster)
   else
      write(*,*) "percwrite: maximum cluster size is ", foundcluster(1)
   end if

   write(*,*) 'percwrite: opening file ', percfile

   open(unit=iounit, file=percfile)

   write(*,*) 'percwrite: writing data ...'

   write(fmtstring, fmt='(''('', ''i1,'',i5,''(1x, i1))'')') pixperline-1
   write(iounit,fmt='(''P2'')')
   write(iounit,*) m,  n

   write(iounit,*) lncluster

   npix = 0

   do j = n, 1, -1
      do i = 1, m

         clusterid = map(i,j)

         !  Write out the largest cluster(s), shading appropriately

         colour = 0

         if (clusterid > 0) then

            do icluster = 1, lncluster

               if (clustersize(clusterid) == foundcluster(icluster)) then

                  !  Largest (first) cluster is white

                  colour = lncluster - icluster + 1

               end if

            end do

         end if

         npix = npix + 1

         pgmline(npix) = colour

         if (npix == pixperline) then
            write(iounit,fmt=fmtstring) pgmline(1:npix)
            npix = 0
         end if

      end do
   end do

   if (npix /= 0) then
      write(iounit,fmt=fmtstring) pgmline(1:npix)
   end if

   write(*,*) 'percwrite: ... done'

   close(iounit)
   write(*,*) 'percwrite: file closed'
 end
end

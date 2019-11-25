module cli_info
  !
  ! Module containing some routines for displaying cli
  ! information.
  !

  private

  character(len=6), parameter :: VERSION = "v0.1.0"

  character(len=57), dimension(50) :: HELP_MSG = (/ &
  "percolate " // VERSION // &
  "                                        ", &
  "                                                        ", &
  "Program for computing, whether any cluster in a         ", &
  "randomly generated matrix percolates.                   ", &
  "                                                        ", &
  "Usage: percolate_[par|ser] [seed] [options ...]         ", &
  "                                                        ", &
  "                                                        ", &
  "seed: INT Seed for the random number generator.         ", &
  "          DEFAULT: 1564                                 ", &
  "                                                        ", &
  "                                                        ", &
  "options:                                                ", &
  "                                                        ", &
  "    -h, --help                    Prints this help mes- ", &
  "                                  sage.                 ", &
  "                                                        ", &
  "        --version                 Prints the version of ", &
  "                                  this program.         ", &
  "                                                        ", &
  "    -l, --length            INT   Sets the dimension of ", &
  "                                  the matrix.           ", &
  "                                  DEFAULT: 288.         ", &
  "                                                        ", &
  "    -d, --density           FLOAT Sets the density of   ", &
  "                                  the full cells.       ", &
  "                                  DEFAULT: 0.411        ", &
  "                                                        ", &
  "    -p, --print_n_clusters  INT   Sets the amount of    ", &
  "                                  clusters displayed in ", &
  "                                  the .pgm file.        ", &
  "                                  DEFAULT: max.         ", &
  "                                  Can display a maximum ", &
  "                                  of the 9 biggest      ", &
  "                                  clusters.             ", &
  "                                                        ", &
  "        --pgm_file_path     PATH  Sets the path for the ", &
  "                                  .pgm file.            ", &
  "                                  DEFAULT: map.pgm      ", &
  "                                                        ", &
  "        --print_iter_factor FLOAT Sets the factor of how", &
  "                                  often the average cell", &
  "                                  value during cluster- ", &
  "                                  ing is printed.       ", &
  "                                  iter % int(L * factor)", &
  "                                  DEFAULT: 0.5          ", &
  "                                  This flag is irrele-  ", &
  "                                  vant for the serial   ", &
  "                                  version of percolate. ", &
  "                                                        "  &
  /)

  public write_help_msg
  public write_version

contains

  subroutine write_help_msg()
    integer :: i

    do i = 1, size(HELP_MSG)
      write (*, *) HELP_MSG(i)
    end do
  end


  subroutine write_version()
    write (*, *) "percolate ", VERSION
  end
end

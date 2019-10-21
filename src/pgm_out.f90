module pgm_out

  implicit none

contains

  !  Function to write a percolation map in greyscale Portable Grey
  !  Map (PGM) format. The largest "ncluster" clusters are identified
  !  and shown as shades of grey against a black background, with the
  !  largest cluster shown in white.

  subroutine pgm_write(percfile, map, ncluster)

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

   return

 end
end

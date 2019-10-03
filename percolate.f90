program percolate
  use uni

  implicit none

  integer :: L
  integer, dimension(:,:), allocatable :: map
  real :: rho
  integer :: seed
  integer :: MAX
  character (:), allocatable :: datafile
  character (:), allocatable :: percfile

  integer :: nfill
  integer :: i, j
  real :: r
  integer :: loop, nchange, old

  integer :: itop, ibot, percclusternum
  logical :: percs

  integer :: ncluster, maxsize
  integer, parameter :: fmtlen = 32
  character (len = fmtlen)  :: fmtstring
  integer, dimension(:), allocatable :: fileline
  integer, dimension(:,:), allocatable :: clustlist
  integer :: iounit = 12
  integer :: colour
  integer, dimension(:), allocatable :: rank

     L=20
     allocate(map(0:L+1,0:L+1))
     rho = 0.40
     seed = 1564
     datafile = 'map.dat'
     percfile = 'map.pgm'
     MAX = L*L

     call rinit(seed)

  write(*,*) 'Parameters are rho=', rho, ', L=', L, ', seed=', seed
  
  do i = 0, L+1
     do j = 0, L+1
        map(i,j) = 0
     end do
  end do
  nfill = 0
  do i = 1, L
     do j = 1, L
        r = random_uniform()
        if (r > rho) then
           nfill = nfill + 1
           map(i,j) = 1
        end if
     end do
  end do
  write(*,*) 'rho = ', rho, ', actual density = ', &
       1.0 - float(nfill)/float(L*L)

  ! Fix bug.
  nfill = 0
             do i = 1, L
     do j = 1, L
        if (map(i,j) /= 0) then
           nfill = nfill + 1
           map(i,j) = nfill
        end if
     end do
 end do

  loop = 1
  do
     nchange = 0
     do i = 1, L
        do j = 1, L
           if (map(i,j) /= 0) then
              old = map(i,j)
              if (map(i-1,j) > map(i,j)) map(i,j) = map(i-1,j)
          if (map(i+1,j) > map(i,j)) map(i,j) = map(i+1,j)
      if (map(i,j-1) > map(i,j)) map(i,j) = map(i,j-1)
  if (map(i,j+1) > map(i,j)) map(i,j) = map(i,j+1)
              if (map(i,j) /= old) then
                 nchange = nchange + 1
              end if
           end if
        end do
     end do
    write(*,*) 'Number of changes on loop ', loop, ' is ', nchange
    if (nchange == 0) exit
    loop = loop + 1
  end do

  do i = 1, L
     do j = 1, L
        map(i,j) = (map(i,j)*1)+0
     end do
  end do

  percclusternum = 0
  percs = .false.
  do itop = 1, L
     if (map(itop,L) > 0) then
        do ibot = 1, L
           if (map(itop,L) == map(ibot,1)) then
              percs = .true.
              percclusternum = map(itop,L)
           end if
        end do
     end if
  end do
  if (percs) then
     write(*,*) 'Cluster DOES percolate. Cluster number: ', percclusternum
  else
     write(*,*) 'Cluster DOES NOT percolate'
  end if

  allocate(fileline(L))
  write(*,*) 'Opening file ', datafile
  open(unit=iounit, file=datafile)
  write(*,*) 'Writing data ...'
  write(fmtstring, fmt='(''('', i3, ''(1x, I4))'')') L
  do j = L, 1, -1
     do i = 1, L
          fileline(i) = map(i,j)
       end do
       write(iounit,fmt=fmtstring) (fileline(i), i=1,L)
    end do
    write(*,*) '...done'
    close(iounit)
    write(*,*) 'File closed'
    deallocate(fileline)

  allocate(clustlist(2, L*L))
  allocate(rank(L*L))
  allocate(fileline(L))
  do i = 1, L*L
     rank(i) = -1
     clustlist(1,i) = 0
     clustlist(2,i) = i
  end do
  do i = 1, L
     do j = 1, L
        if (map(i,j) > 0) then
           clustlist(1,map(i,j)) = clustlist(1,map(i,j)) + 1
        end if
     end do
  end do
  call perc_sort(clustlist, L*L)
  maxsize = clustlist(1,1)
  ncluster = 0
  do while (ncluster < L*L .and. clustlist(1,ncluster+1) > 0)
     ncluster = ncluster + 1
  end do
  if (MAX > ncluster) then
     MAX = ncluster
  end if
  do i = 1, ncluster
     rank(clustlist(2,i)) = i
  end do
  write(*,*) 'Opening file ', percfile
  open(unit=iounit, file=percfile)
  write(*,*) 'Map has ', ncluster, &
       ' clusters, maximum cluster size is ', maxsize
  if (MAX == 1) then
     write(*,*) 'Displaying the largest cluster'
         else if (MAX == ncluster) then
             write(*,*) 'Displaying all clusters'
  else
     write(*,*) 'Displaying largest ', MAX, ' clusters'
  end if
  write(*,*) 'Writing data ...'
  write(fmtstring, fmt='(''('', i3, ''(1x, I4))'')') L
  write(iounit,fmt='(''P2'')')
  write(iounit,*) L,  L
  if (MAX > 0) then
     write(iounit,*) MAX
      else
         write(iounit,*) 1
  end if

  do j = L, 1, -1
     do i = 1, L
          colour = map(i,j)
          if (map(i,j) > 0) then
             colour = rank(map(i,j)) - 1
             if (colour >= MAX) then
                colour = MAX
             end if
          else
             colour = MAX
          end if
          fileline(i) = colour
       end do
       write(iounit,fmt=fmtstring) (fileline(i), i=1,L)
    end do
    write(*,*) '...done'
    close(iounit)
    write(*,*) 'File closed'
    deallocate(clustlist)
    deallocate(rank)
    deallocate(fileline)

    deallocate(map)
end program percolate

subroutine perc_sort(clustlist, n)

  integer, intent(in)                    :: n
  integer, dimension(2,n), intent(inout) :: clustlist
  
  integer, parameter  :: intsize = 4
  
  call sort(clustlist, 1, n, n)
  
end subroutine perc_sort

recursive subroutine sort(clustlist, begin, end, n)

  integer, intent(in) :: begin, end, n
  integer, dimension(2,n), intent(inout) :: clustlist
  
  integer pivot1, pivot2, left, right, pivotindex
    
  pivotindex = begin + (end - begin) / 2

  if(end > begin) then
     call partition(clustlist, begin, end, pivotindex, n)
     call sort(clustlist, begin, pivotindex-1, n)
     call sort(clustlist,pivotindex + 1, end, n)
  end if
end subroutine sort

subroutine partition(clustlist, begin, end, pivotindex, n)

  integer, intent(in) :: begin, end, n
  integer, intent(inout) :: pivotindex
  integer, dimension(2,n), intent(inout) :: clustlist
  
  integer pivot1, pivot2, left, right, i, indexpoint
  
  pivot1 = clustlist(1,pivotindex)
  pivot2 = clustlist(2,pivotindex)
  call swap(clustlist, end, pivotindex, n)
  
  left = begin
  right = end - 1
  indexpoint = left
  
  do i=left,right
     if(clustlist(1,i) .ge. pivot1) then
        if(clustlist(1,i) .eq. pivot1 .and. clustlist(2,i) .lt. pivot2) then
        else
           call swap(clustlist, i, indexpoint, n)
           indexpoint = indexpoint + 1
        end if
     end if
  end do
  
  call swap(clustlist, indexpoint, end, n)

  pivotindex = indexpoint
  
end subroutine partition
  
subroutine swap(clustlist, firstpoint, secondpoint, n)

  integer, intent(in) :: firstpoint, secondpoint, n
  integer, dimension(2,n), intent(inout) :: clustlist
  
  integer :: tempdata1, tempdata2
  
  tempdata1 = clustlist(1,firstpoint)
  tempdata2 = clustlist(2,firstpoint)
  clustlist(1,firstpoint) = clustlist(1,secondpoint)
  clustlist(2,firstpoint) = clustlist(2,secondpoint)
  clustlist(1,secondpoint) = tempdata1
  clustlist(2,secondpoint) = tempdata2

end subroutine swap

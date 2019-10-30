! moritz siegel 191016
! sudoku solver via simulated annealing & heat bath / metropolis

! QUESTIONS:
! deviation from equilibral energy of 55*9*3 = 1215


program sadoku
  
  implicit none
  
  integer, parameter :: imax=9, jmax=9, nn=8
  integer, dimension(9), parameter :: im=[2,5,8,2,5,8,2,5,8] ! koords of middle items
  integer, dimension(9), parameter :: jm=[2,2,2,5,5,5,8,8,8]
  integer, dimension(9), parameter :: inn=[1,-1,0,0,1,1,-1,-1] ! nearest & next nearest neighbour coordinates
  integer, dimension(9), parameter :: jnn=[0,0,1,-1,1,-1,1,-1]

!  integer, dimension (imax,jmax) :: mr, seg, cor, err, erc, ers
!  integer :: i,j,k,c, idummy,jdummy,inew,jnew, segnew
!  integer :: flips, accept, errtot, errata(5), correct(5)
!  real*8 :: r, t, enew,eold, perc, dist
!  real*8, parameter :: jj=0.99d0, lambda=0.9d0, ti=10000d0, tf=0.001d0
!  real*8, parameter :: far=0.20d0 ! coupling factor of the next nearest neighbours
  logical, dimension(imax,jmax) :: mask
  
  ! init
  t = ti
  c = 1
  flips = 0
  err(:,:) = 0
  erc(:,:) = 0
  ers(:,:) = 0
  errata(:) = 0
  correct(:) = 0
  accept = 0
  
  call srand(time()) ! initialise random numbers

  ! read sudoku, define mask (empty slots), random-seed those
  open(unit=10, file='sudoku_in.dat')
  do i = 1,imax
     do j = 1,jmax
        read() sudoku(i,j)
        if ( sudoku_in(i,j) .eq. 0 ) then ! def mask := free slots
           mask(i,j) = .true.
           sudoku(i,j)  = int(9.d0*rand()) + 1
        else
           mask(i,j) = .false.
        end if
     end do
  end do
  close(unit=10)

  eold = energy(sudoku)

  ! divide & rule:  count errata
  rows(:,:) = 0
  cols(:,:) = 0
  patch(:,:) = 0
  do i = 1, imax
     do j = 1, jmax
        ! row sums
        if ( rows(i, sudoku(i,j)) .eq. 0 ) then ! ask if number is still free
           rows(i, sudoku(i,j)) = 1
        else
           errata(i,j) = errata(i,j) + 1 ! increment count if number is already taken
        end if
        
        ! col sums
        if ( cols(j, sudoku(i,j)) .eq. 0 ) then
           cols(j, sudoku(i,j)) = 1
        else
           errata(i,j) = errata(i,j) + 1
        end if

        ! patch sums, using i & j as patch & next neighbours
        imm = im(i)
        jmm = jm(i)
        inew = imm + inn(j)
        jnew = jmm + jnn(j)
        if ( patch(i, sudoku(inew,jnew)) .eq. 0 ) then
           patch(i, sudoku(inew,jnew)) = 1
        else
           errata(i,j) = errata(i,j) + 1
        end if
     end do
  end do
  
  numbers 
     colsum = sum(sudoku(m,:)

     ! patch sums
     patchsum = sudoku(i,j)
     imm = im(m)
     jmm = jm(m)
     nn : do k = 1, nn ! loop over nearest neighbors
        inew = i + inn(k)
        jnew = j + jnn(k)
        patchsum = patchsum + sudoku(inew,jnew)
     end do nn
     
     ! change in energy
     energy = abs(rowsum + colsum + patchsum - de)
  end do divide
  
              
  ! sweeeeep
  open(unit=20, file='flips.dat')
  cool : do while (t .gt. tf)
     row : do i = 1, imax
        col : do j = 1, jmax
           
           if ( mask(i,j) .eqv. .false. ) break ! skip if the number is fixed parameter
           
           ! flip pixel to random except old value
           guess(i,j) = mod(sudoku(i,j) + int(8.d0*rand()) , 9) + 1
                    
           ! change in energy
              energy = abs(energy - 3*sudoku(i,j) + + colsum + patchsum - de)
              
              ! check for floating point overflow & define r
              if ( de/t .gt. 20.d0 ) then
                 r = 0.d0
              else if ( de/t .lt. 0.001d0 ) then
                 r = 1.d0
              else
                 r = exp(-de/t)
              end if
              
              ! accept or deny & metropolis algorithm
              if ( rand() .lt. r ) then
                 sudoku(i,j) = guess
                 flips = flips + 1
              end if
           end do divide
        end do col
     end do row
     
     ! write & reset flips
     write(*,*) 'run: ',c, 'temperature: ',t ,'flips: ',flips
     write(20,*) c, t, flips
     accept = accept + flips
     flips = 0
     
     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + 1
     
  end do cool
  close(unit=20)
  
  ! write image & compare
  open(unit=30, file='sudoku_out.dat')
  do i=1,imax
     do j=1,jmax
        write(30,'(3i5)') seg(i,j)
        write(*,*) seg(i,j)
     end do
     write(*,*)
  end do
  close(unit=30)
  
  write(*,*) 'runs: ',c, 'total flips: ',accept
  
end program sadoku

! moritz siegel 191016
! sudoku solver via simulated annealing & heat bath / metropolis
! ij : 11 12 13 ...
!      21 22 
!      31    ...
!      ...

! OPEN ISSUES:

program sadoku
  
  implicit none
  
  logical, dimension(9,9) :: mask ! predefined sudoku values
  integer, dimension(9), parameter :: im=[2,2,2,5,5,5,8,8,8] ! middlemen
  integer, dimension(9), parameter :: jm=[2,5,8,2,5,8,2,5,8]
  integer, dimension(9), parameter :: in=[-1,-1,-1,0,0,0,1,1,1] ! nearest neighbours
  integer, dimension(9), parameter :: jn=[-1,0,1,-1,0,1,-1,0,1]
  integer, dimension (9,9) :: sudoku, row, col, prt
  integer :: i, j, m, n
  integer :: flips, accept, energy
  real*8 :: r, t
  

  call srand(time()) ! initialise random numbers

  
  ! read sudoku, define mask of empty slots, random-seed those
  open(unit=10, file='sudoku.dat')
  do i = 1,9
     do j = 1,9
        read(10) sudoku(i,j)
        if ( sudoku(i,j) .eq. 0 ) then ! mask 1 := free slot
           mask(i,j) = .true.
           sudoku(i,j)  = int(9.d0*rand()) + 1
        else
           mask(i,j) = .false.
        end if
     end do
  end do
  
  
  ! divide & rule:  count energy
  row(:,:) = -1 ! minus 1 := free slot
  col(:,:) = -1
  prt(:,:) = -1
  patch : do m = 1, 9 ! select patch middlemen
     element : do n = 1, 9 ! cycle through neighbours within selected patch
        i = im(m) + in(n)
        j = jm(m) + jn(n)
        row(i,sudoku(i,j)) = row(i,sudoku(i,j)) + 1
        col(j,sudoku(i,j)) = col(j,sudoku(i,j)) + 1
        prt(m,sudoku(i,j)) = prt(m,sudoku(i,j)) + 1
     end do element
  end do patch
  
  
  ! sweeeeep
  t = ti
  c = 1
  flips = 0
  accept = 0
  open(unit=20, file='flips.dat')
  open(unit=30, file='energy.dat')
  cool : do while ( t .gt. tf )
     patch : do m = 1, 9 ! select patch
        element : do n = 1, 9 ! cycle through selected patch
           i = im(m) + in(j)
           j = jm(m) + jn(j)
           
           if ( mask(i,j) .eqv. .false. ) break ! skip if current number is fixed parameter
           
           ! flip field to random except old value
           guess = mod(sudoku(i,j) + int(8.d0*rand()) , 9) + 1
           
           ! update rows & cols & patch & energy
           drow = - abs(row(i,sudoku(i,j))) + abs(row(i,sudoku(i,j)) - 1) - abs(row(i,guess)) + abs(row(i,guess) + 1)
           dcol = - abs(col(j,sudoku(i,j))) + abs(col(j,sudoku(i,j)) - 1) - abs(col(j,guess)) + abs(col(j,guess) + 1)
           dprt = - abs(prt(m,sudoku(i,j))) + abs(prt(m,sudoku(i,j)) - 1) - abs(prt(m,guess)) + abs(prt(m,guess) + 1)
           de = drow + dcol + dprt
           
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
              row(i,sudoku(i,j)) = row(i,sudoku(i,j)) - 1
              row(i,guess) = row(i,guess) + 1
              col(j,sudoku(i,j)) = col(j,sudoku(i,j)) - 1
              col(j,guess) = col(j,guess) + 1
              prt(m,sudoku(i,j)) = prt(m,sudoku(i,j)) - 1
              prt(m,guess) = prt(m,guess) + 1
              flips = flips + 1
              energy = energy + de
           end if
        end do element
     end do patch
     
     ! write & reset flips
     write(*,*) 'run: ',c, 'temperature: ',t ,'flips: ',flips, 'energy: ',energy
     write(20,*) c, t, flips, energy
     accept = accept + flips
     flips = 0
     
     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + 1
     
  end do cool

  
  ! write image & compare
  open(unit=40, file='sadoku.dat')
  do i=1,9
     do j=1,9
        write(*,*) sudoku(i,j)
        write(40,'(3i5)') sudoku(i,j)
     end do
     write(*,*)
  end do

  
  write(*,*) 'runs: ',c, 'total flips: ',accept, 'final energy: ',energy

  
end program sadoku

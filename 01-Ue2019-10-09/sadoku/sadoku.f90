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
  integer, dimension (9,9) :: sudokuin, sudoku, col, row, prt
  integer :: i, j, m, n, c
  integer :: guess, de, dcol, drow, dprt, flips, accept, energy
  real*8, parameter ::   tf=0.1, ti=60, lambda=0.99
  real*8 :: r, t
  

  call srand(time()) ! initialise random numbers

  ! read sudoku, define mask of empty slots, random-seed those
  open(unit=10, file='sudoku.dat') ! assuming its a matrix of space separated lines of integers in ascii
  read(10,*) sudokuin
  sudoku = sudokuin
  do i = 1,9
     do j = 1,9
        if ( sudoku(i,j) .eq. 0 ) then ! mask 1 := free slot
           mask(i,j) = .true.
           sudoku(i,j)  = int(9.d0*rand()) + 1
        else
           mask(i,j) = .false.
        end if
     end do
  end do
  write(*,*) "sudokuin"
  write(*,*) ( ( sudokuin(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "mask"
  write(*,*) ( ( mask(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "initial sudoku"
  write(*,*) ( ( sudoku(i,j), i=1,9 ), new_line('A'), j=1,9 )

  ! divide & rule:  count energy
  col(:,:) = -1 ! minus 1 := free slot
  row(:,:) = -1
  prt(:,:) = -1
  do m = 1, 9 ! select patch middlemen
     do n = 1, 9 ! cycle through neighbours within selected patch
        i = im(m) + in(n)
        j = jm(m) + jn(n)
        col(i,sudoku(i,j)) = col(i,sudoku(i,j)) + 1
        row(j,sudoku(i,j)) = row(j,sudoku(i,j)) + 1
        prt(m,sudoku(i,j)) = prt(m,sudoku(i,j)) + 1
     end do
  end do
  energy = sum(abs(col(:,:))) + sum(abs(row(:,:))) + sum(abs(prt(:,:)))
  write(*,*) "col 1--9, numbers 1|9"
  write(*,*) ( ( col(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "row 1--9, numbers 1|9"
  write(*,*) ( ( row(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "prt 1--9, numbers 1|9"
  write(*,*) ( ( prt(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "initial energy", energy
  
  ! sweeeeep
  t = ti
  c = 1
  flips = 0
  accept = 0
  open(unit=20, file='stats.dat')
  open(unit=30, file='energy.dat')
  cool : do while ( t .gt. tf )
     patch : do m = 1, 9 ! select patch
        element : do n = 1, 9 ! cycle through selected patch
           i = im(m) + in(j)
           j = jm(m) + jn(j)
           
           if ( mask(i,j) .eqv. .false. ) cycle ! skip if current number is fixed parameter
           
           ! flip field to random except old value
           guess = mod(sudoku(i,j) + int(8.d0*rand()) , 9) + 1
           
           ! update cols & rows & patch & energy
           dcol = - abs(col(i,sudoku(i,j))) + abs(col(i,sudoku(i,j)) - 1) - abs(col(i,guess)) + abs(col(i,guess) + 1)
           drow = - abs(row(j,sudoku(i,j))) + abs(row(j,sudoku(i,j)) - 1) - abs(row(j,guess)) + abs(row(j,guess) + 1)
           dprt = - abs(prt(m,sudoku(i,j))) + abs(prt(m,sudoku(i,j)) - 1) - abs(prt(m,guess)) + abs(prt(m,guess) + 1)
           de = dcol + drow + dprt ! -6 < de < 6
           write(*,*) "de:", de
           
           ! check for floating point overflow & define r
           if ( de/t .gt. 20.d0 ) then
              r = 0.d0
           else if ( de/t .lt. 0.001d0 ) then
              r = 1.d0
           else
              r = exp(-de/t)
           end if
           write(*,*) "r:",r, "-de/t", -de/t
          
           ! accept or deny & metropolis algorithm
           if ( rand() .lt. r ) then
              col(i,sudoku(i,j)) = col(i,sudoku(i,j)) - 1
              col(i,guess) = col(i,guess) + 1
              row(j,sudoku(i,j)) = row(j,sudoku(i,j)) - 1
              row(j,guess) = row(j,guess) + 1
              prt(m,sudoku(i,j)) = prt(m,sudoku(i,j)) - 1
              prt(m,guess) = prt(m,guess) + 1
              sudoku(i,j) = guess
              flips = flips + 1
              energy = energy + de
              write(*,*) "energy:", energy         
           end if
        end do element
     end do patch
     
     ! write & reset flips
     write(*,*) 'run: ',c, 'temperature: ',t ,'flips/cell: ',flips/dble(81), 'energy: ',energy
     write(20,*) c, t, flips, energy
     accept = accept + flips
     flips = 0
     
     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + 1
     
  end do cool
  
  write(*,*) "final sudoku"
  write(*,*) ( ( sudoku(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "col 1--9, numbers 1|9"
  write(*,*) ( ( col(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "row 1--9, numbers 1|9"
  write(*,*) ( ( row(i,j), i=1,9 ), new_line('A'), j=1,9 )
  write(*,*) "prt 1--9, numbers 1|9"
  write(*,*) ( ( prt(i,j), i=1,9 ), new_line('A'), j=1,9 )
  
  write(*,*) 'runs: ',c, 'total flips: ',accept, 'final energy: ',energy

  
end program sadoku

! moritz siegel 191016
! sudoku solver via simulated annealing & heat bath / metropolis
! ij : 11 12 13 ...
!      21 22 
!      31    ...
!      ...

! OPEN ISSUES:

program sadoku

  implicit none
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  real*8, parameter :: tf=0.001, ti=90, lambda=0.99
  ! ////////////////////////////////////////////////////////////////////////////
  logical, dimension(9,9) :: mask, num ! predefined sudoku values
  integer, dimension(9), parameter :: im=[2,2,2,5,5,5,8,8,8] ! middlemen
  integer, dimension(9), parameter :: jm=[2,5,8,2,5,8,2,5,8]
  integer, dimension(9), parameter :: in=[-1,-1,-1,0,0,0,1,1,1] ! nearest neighbours
  integer, dimension(9), parameter :: jn=[-1,0,1,-1,0,1,-1,0,1]
  integer, dimension (9,9) :: sudokuin, sudoku, col, row, prt, colte, rowte, prtte, sudote, coltemp, rowtemp
  integer :: i, j, m, n, ai, aj, bi, bj, c, dummy, o, p, same, chains
  integer :: a, b, de, dcol, drow, flips, accept, energy, eold, reheat
  real*8 :: r, t

  call srand(time()) ! initialise random numbers
  
  call seedoku()

  ! sweeeeep
  t = ti
  c = 0
  flips = 0
  accept = 0
  reheat = 0
  same = 0
  eold = energy
  open(unit=20, file='stats.dat')
  cool : do while ( t .gt. tf )
     markov : do n = 1,chains
        ! select patch
        m = int(8.d0*rand()) + 1
        
        ! select two different fields in that patch
        a = 0
        do ! if number is fixed state, guess again
           a = mod(a + int(8.d0*rand()) , 9) + 1
           ai = im(m) + in(a)
           aj = jm(m) + jn(a)
           if ( mask(ai,aj) .eqv. .false. ) exit
        end do
        do ! if number is fixed state, guess again
           b = mod(a + int(8.d0*rand()) , 9) + 1
           bi = im(m) + in(b)
           bj = jm(m) + jn(b)
           if ( mask(bi,bj) .eqv. .false. ) exit
        end do

        ! update cols & rows & energy
        coltemp = col
        rowtemp = row
        coltemp(ai,sudoku(ai,aj)) = coltemp(ai,sudoku(ai,aj)) -1
        coltemp(ai,sudoku(bi,bj)) = coltemp(ai,sudoku(bi,bj)) +1
        coltemp(bi,sudoku(ai,aj)) = coltemp(bi,sudoku(ai,aj)) +1
        coltemp(bi,sudoku(bi,bj)) = coltemp(bi,sudoku(bi,bj)) -1
        rowtemp(aj,sudoku(ai,aj)) = rowtemp(aj,sudoku(ai,aj)) -1
        rowtemp(aj,sudoku(bi,bj)) = rowtemp(aj,sudoku(bi,bj)) +1
        rowtemp(bj,sudoku(ai,aj)) = rowtemp(bj,sudoku(ai,aj)) +1
        rowtemp(bj,sudoku(bi,bj)) = rowtemp(bj,sudoku(bi,bj)) -1
        de = sum(abs(coltemp(:,:))) + sum(abs(rowtemp(:,:))) &
             - sum(abs(col(:,:))) - sum(abs(row(:,:)))
        
        write(*,*) "de soll", de
        
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
           ! update row & col
           col = coltemp
           row = rowtemp
           ! swap the two numbers
           dummy = sudoku(ai,aj)
           sudoku(ai,aj) = sudoku(bi,bj)
           sudoku(bi,bj) = dummy
           ! count flips & update energy
           flips = flips + 1
           energy = sum(abs(col(:,:))) + sum(abs(row(:,:)))
           write(*,*) "energy", energy

           ! if solved, do not solve on
           if (energy .eq. 0 ) then
              write(*,*) "solved sudoku"
              write(*,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
              write(*,*) "col 1--9, numbers 1|9"
              write(*,'(9i2)') ( ( col(i,j), i=1,9 ), j=1,9 )
              write(*,*) "row 1--9, numbers 1|9"
              write(*,'(9i2)') ( ( row(i,j), i=1,9 ), j=1,9 )
              write(*,*) 'runs: ',c*chains+n, 'total flips: ',accept, 'final energy: ',energy

              stop

           end if
        end if

     end do markov
     
     ! count number of chains that do not change in energy
     if ( energy .eq. eold  ) then
        same = same + 1
     else
        same = 0
     end if
     eold = energy

     ! reheat if stuck / solution not obtained
     if ( same .gt. 50 .and. reheat .lt. 10) then
        call seedoku() ! reseed
        same = 0
        reheat = reheat + 1
        t = ti
        write(*,*) 'chain: ',n, 'reheating'
        !write(20,*) ! seperate blocks of data
        cycle cool
     end if

     write(*,*) "sudoku"
     write(*,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
     write(*,*) "col 1--9, numbers 1|9"
     write(*,'(9i2)') ( ( col(i,j), i=1,9 ), j=1,9 )
     write(*,*) "row 1--9, numbers 1|9"
     write(*,'(9i2)') ( ( row(i,j), i=1,9 ), j=1,9 )

     ! write & reset flips
     write(*,*) 'run: ',c*chains+n, 'temperature: ',t ,'flips: ',flips, 'energy: ',energy
     write(20,*) c, t, flips, energy
     accept = accept + flips
     flips = 0

     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + 1

  end do cool


  write(*,*) "final sudoku"
  write(*,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
  write(*,*) "col 1--9, numbers 1|9"
  write(*,'(9i2)') ( ( col(i,j), i=1,9 ), j=1,9 )
  write(*,*) "row 1--9, numbers 1|9"
  write(*,'(9i2)') ( ( row(i,j), i=1,9 ), j=1,9 )
  write(*,*) 'chains: ',chains, 'total runs: ',c*chains+n, 'total reheats:', reheat, 'total flips: ',accept, 'final energy: ',energy
  
  
contains


  subroutine seedoku() ! read sudoku & define mask of empty slots & random-seed sudoku
    implicit none
    open(unit=10, file='sudoku.dat') ! assuming its a matrix of space separated lines of integers in ascii
    read(10,*) sudokuin
    sudoku = sudokuin
    num(:,:) = .false.
    chains = 0
    do m = 1, 9 ! select patch middlemen
       do n = 1, 9 ! cycle through neighbours within selected patch
          i = im(m) + in(n)
          j = jm(m) + jn(n)

          ! read sudoku & discard fixed state numbers
          if ( sudoku(i,j) .eq. 0 ) then ! mask 0 := free slot
             mask(i,j) = .false.
          else
             mask(i,j) = .true.
             num(m,sudoku(i,j)) = .true.
             chains = chains + 1
          end if
       end do
       do n = 1, 9 ! cycle through neighbours within selected patch
          i = im(m) + in(n)
          j = jm(m) + jn(n)
          ! write empty sudoku cells
          if ( mask(i,j) .eqv. .false. ) then
             do
                sudoku(i,j)  = int(9.d0*rand()) + 1
                if ( num(m,sudoku(i,j)) .eqv. .true. ) then ! number is already present, guess again
                   cycle 
                else ! take that
                   num(m,sudoku(i,j)) = .true.
                   exit
                end if
             end do
          end if
       end do
    end do
    write(*,*) "sudokuin"
    write(*,'(9i2)') ( ( sudokuin(i,j), i=1,9 ), j=1,9 )
    write(*,*) "mask"
    write(*,'(9l2)') ( ( mask(i,j), i=1,9 ), j=1,9 )
    write(*,*) "initial sudoku"
    write(*,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
    
    chains = chains*chains
    
    ! divide & rule:  count energy
    col(:,:) = -1
    row(:,:) = -1
    do i = 1, 9
       do j = 1, 9
          col(i,sudoku(i,j)) = col(i,sudoku(i,j)) + 1
          row(j,sudoku(i,j)) = row(j,sudoku(i,j)) + 1
       end do
    end do
    energy = sum(abs(col(:,:))) + sum(abs(row(:,:)))
    write(*,*) "initial energy", energy
    write(*,*) "col 1--9, numbers 1|9"
    write(*,'(9i2)') ( ( col(i,j), i=1,9 ), j=1,9 )
    write(*,*) "row 1--9, numbers 1|9"
    write(*,'(9i2)') ( ( row(i,j), i=1,9 ), j=1,9 )
    
    close(unit=10)
  end subroutine seedoku


end program

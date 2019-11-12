! moritz siegel 191112
! sudoku solver via simulated annealing & heat bath / metropolis & reheats & gradually slows

program sadoku

  implicit none
  ! \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
  integer, parameter :: maxiter = 10000000 ! absolute last limit, please do not exceed 2147483647
  real*8, parameter :: frozen = 10 ! reheat if energy does not change after frozen chains steps
  real*8, parameter :: rehot = 20 ! limit reheats
  real*8 :: hot = 15 ! initial temperatur, is later automatically set to ~80% initial acceptance 
  real*8 :: lambda = 0.5 ! define initial cooling speed, set to 0.1<x<0.99
  real*8 :: chains = 1 ! adjust markov chain length, theoretic optimum is 1
  integer :: prints = 10000 ! print interval, recomended 100<x<10000
  ! ////////////////////////////////////////////////////////////////////////////
  logical, dimension(9,9) :: mask, num ! predefined sudoku values
  integer, dimension(9), parameter :: im=[2,2,2,5,5,5,8,8,8] ! middlemen
  integer, dimension(9), parameter :: jm=[2,5,8,2,5,8,2,5,8]
  integer, dimension(9), parameter :: in=[-1,-1,-1,0,0,0,1,1,1] ! nearest neighbours
  integer, dimension(9), parameter :: jn=[-1,0,1,-1,0,1,-1,0,1]
  integer, dimension (9,9) :: sudokuin, sudoku, col, row, coltemp, rowtemp
  integer :: i, j, m, n, ai, aj, bi, bj, dummy, o, p, freeze, links
  integer :: a, b, c, de, flips, accept, energy, reheat, eold
  real*8, dimension (100) :: des
  real*8 :: r, t, u, cold, desmean, sigma 
  logical :: est
  
  ! init
  open(unit=20, file='stats.dat')
  open(unit=30, file='sadoku.dat')
  t = hot
  c = 0
  flips = 0
  accept = 0
  reheat = 0
  freeze = 0
  eold = energy
  des = 0
  est = .true.
  cold = 0
  
  call random_seed() ! initialise random numbers
  
  call seedoku()

  ! sweeeeep
  cool : do while ( c .lt. maxiter )
     accept = accept + flips
     flips = 0
     markov : do n = 1,links

        ! select patch & select two different fields in that patch
        call random_number(u)
        m = int(9*u) + 1
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

        ! estimate mean & deviation of initial de & set hot accordingly & restart
        if ( est .eqv. .true. ) then 
           des(n) = de
           if ( n .gt. 81 ) then ! only first 100 cells of first markov chain
              desmean = sum(des)*0.01
              sigma = sqrt( sum( (des-desmean)**2 ) * 0.01 )
              hot = - (abs(desmean) + sigma) / log(0.8) ! 80% of the initial cells (mean+deviation) should flip
              cold = (abs(desmean) + sigma) / log(links*chains)
              t = hot
              est = .false.
              write(*,*) 'estimating temperatur poles: hot=', hot, 'cold=', cold, '>> restarting'
              cycle cool
           end if
        end if
        
        ! check for floating point overflow & define r
        if ( de/t .gt. 10 ) then
           r = 0
        else if ( de/t .lt. 0.001 ) then
           r = 1
        else
           r = exp(-de/t)
        end if

        ! accept or deny & metropolis algorithm
        call random_number(u)
        if ( u .lt. r ) then
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

           ! if solved, do not solve on
           if (energy .eq. 0 ) then
              call printoku('sudoku solved!')
              write(30,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
              stop
           end if
        end if
        
        ! print stats / whole sudoku snapshots & stats
        if ( mod(n, prints) .eq. 0 ) then
           write(20,*) c + n, t, flips/real(prints), energy
           accept = accept + flips
           flips = 0
        end if
        !call printoku('snapshot sudoku')
     
     end do markov

     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + links
     
     ! count number of chains that do not change in energy
     if ( energy .eq. eold  ) then
        freeze = freeze + 1
     else
        freeze = 0
     end if
     eold = energy
     
     ! stuck / no solution obtained so far / out of bounds
     if ( freeze .gt. frozen .and. reheat .lt. rehot ) then
        freeze = 0
        call heatoku('stuck >> reseeding & reheating')
        cycle cool
     else if ( t .lt. cold .and. reheat .lt. rehot ) then
        call heatoku('no solution obtained >> reseeding & reheating')
        cycle cool
     else if ( t .lt. cold .and. reheat .eq. rehot ) then
        write(*,*) 'no solution obtained within maximum reheats >> aborting'
        call printoku('final sudoku')
        stop
     end if
     
  end do cool
  
  
contains


  subroutine seedoku() ! read sudoku & define mask of empty slots & random-seed sudoku
    implicit none
    open(unit=10, file='sudoku.dat') ! assuming its a matrix of space separated lines of integers in ascii
    read(10,*) sudokuin
    sudoku = sudokuin
    num(:,:) = .false.
    links = 0
    do m = 1, 9 ! select patch middlemen
       do n = 1, 9 ! cycle through neighbours within selected patch
          i = im(m) + in(n)
          j = jm(m) + jn(n)

          ! read sudoku & discard fixed state numbers
          if ( sudoku(i,j) .eq. 0 ) then ! mask 0 := free slot
             mask(i,j) = .false.
             links = links + 1
          else
             mask(i,j) = .true.
             num(m,sudoku(i,j)) = .true.
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
    write(*,*) "input sudoku"
    write(*,'(9i2)') ( ( sudokuin(i,j), i=1,9 ), j=1,9 )
    write(*,*) "mask sudoku"
    write(*,'(9l2)') ( ( mask(i,j), i=1,9 ), j=1,9 )

    ! calculate solution space ~freecells**2
    links = int(links*links*chains)

    ! set printer interval to maximum if necessary
    if ( prints .gt. links ) prints = links
    
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
    close(unit=10)

    call printoku('initial sudoku')
    
  end subroutine seedoku

  
  subroutine printoku(text)
    implicit none
    character(len=*), intent(in) :: text
    
    write(*,*) text
    write(*,'(9i2)') ( ( sudoku(i,j), i=1,9 ), j=1,9 )
    write(*,*) "col 1--9, numbers 1|9"
    write(*,'(9i2)') ( ( col(i,j), i=1,9 ), j=1,9 )
    write(*,*) "row 1--9, numbers 1|9"
    write(*,'(9i2)') ( ( row(i,j), i=1,9 ), j=1,9 )
    write(*,*) 'iteration: ',c, ', temperatur:',t, ', energy: ',energy
    write(*,*) 'links: ',links, ', reheats:', reheat, ', accepts: ',accept
    write(*,*) 'parameters: hot=', hot, ', cold=', cold, ', lambda=', lambda
    write(*,*) 'parameters: frozen=', frozen, ', rehot=', rehot, ', chains=', chains

    write(20,*) c, t, flips/real(links), energy
        
  end subroutine printoku

  subroutine heatoku(text)
    implicit none
    character(len=*), intent(in) :: text
    
    lambda = 0.9*lambda + 0.1 ! lower speed
    reheat = reheat + 1
    t = hot
    est = .true.
    write(*,*) text
    write(20,*) new_line('a') ! seperate blocks of data
    
    call seedoku() ! reseed
    
  end subroutine heatoku
  
end program

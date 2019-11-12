! moritz siegel 181007
! ising model with periodic boundaries using monte carlo simulation with metropolis algorithm:
! segment mr image using bayes probability, in nearest and next nearest neighbour configuration (apparantly not better! why?).

program samr
  
  implicit none

  integer, parameter :: imax=254, jmax=333, nn=8 ! only nearest neighbours: set to 4
  integer, dimension(8), parameter :: inn=[1,-1,0,0,1,1,-1,-1] ! nearest & next nearest neighbour coordinates
  integer, dimension(8), parameter :: jnn=[0,0,1,-1,1,-1,1,-1]
  integer, dimension(5), parameter :: mv=[30,426,602,1223,167] ! mean values & standard deviation, [bg,wm,gm,csf,sb]
  integer, dimension(5), parameter :: sigma=[30,59,102,307,69] 
  integer, dimension (imax,jmax) :: mr, seg, cor, err, erc, ers
  integer :: i,j,k,c, idummy,jdummy,inew,jnew, segnew
  integer :: flips, accept, errtot, errata(5), correct(5)
  real*8 :: r, t, enew,eold, perc, dist
  real*8, parameter :: jj=0.99d0, lambda=0.9d0, ti=10000d0, tf=0.001d0
  real*8, parameter :: far=0.20d0 ! coupling factor of the next nearest neighbours
  
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

  ! read mr & random-seed seg
  open(unit=10, file='mr.dat')
  do i = 1,imax
     do j = 1,jmax
        read(10,'(3i5)') idummy,jdummy,mr(i,j)
        seg(i,j)  = int(5.d0*rand()) + 1
     end do
     read(10,*)
  end do
  close(unit=10)
  
  ! sweeeeep
  open(unit=20, file='flips.dat')
  do while (t .gt. tf)
     do i = 1, imax
        do j = 1, jmax
           ! calculate energy of current setting
           eold = 0.5d0*( ((mr(i,j)) - mv(seg(i,j))) &
                / dble(sigma(seg(i,j))) )**2.d0 + log(dble(sigma(seg(i,j))))
           
           ! flip pixel to random except old value
           segnew = mod(seg(i,j) + int(4.d0*rand()) , 5) + 1
           
           ! calculate energy of the new setting
           enew = 0.5d0*( ((mr(i,j) - mv(segnew))) &
                / dble(sigma(segnew)) )**2.d0 + log(dble(sigma(segnew)))
           
           ! loop over nearest neighbors
           do k = 1, nn
              dist = 1.d0
              inew = i + inn(k)
              jnew = j + jnn(k)
                            
              ! check periodic boundary conditions
              if (inew .le. 0) then
                 inew = imax
              else if(inew .gt. imax) then
                 inew = 1
              end if
              if (jnew .le. 0) then
                 jnew = jmax
              else if(jnew .gt. jmax) then
                 jnew = 1
              end if
              
              ! define distance of current to center position
              if ( inn(k) .ne. 0 .and. jnn(k) .ne. 0 ) dist = far
              !dist = norm2([dble(inew),dble(jnew)]) ! correct vector norm for more than 8 neighbours
              
              ! calculate energy to nn of old setting
              if ( seg(i,j) == seg(inew,jnew) ) eold = eold - jj * dist
              
              ! calculate energy to nn of new setting
              if ( segnew == seg(inew,jnew) ) enew = enew - jj * dist
           end do ! nn
           
           ! check for floating point overflow & define r
           if ( (enew-eold)/t .gt. 20.d0 ) then
              r = 0.d0
           else if ( (enew-eold)/t .lt. 0.001d0 ) then
              r = 1.d0
           else
              r = exp(-(enew-eold)/t)
           end if
         
           ! accept or deny & metropolis algorithm
           if ( rand() .lt. r ) then
              seg(i,j) = segnew
              flips = flips + 1
           end if

        end do ! j
     end do ! i

     ! write & reset flips
     write(*,*) 'run: ',c, 'temperature: ',t ,'flips: ',flips
     write(20,*) c, t, flips
     accept = accept + flips
     flips = 0
     
     ! decrease temperature & count sweeps
     t = t * lambda
     c = c + 1

  end do ! t
  close(unit=20)
    
  ! write image & compare
  open(unit=30, file='seg.dat')
  open(unit=40, file='cor.dat')
  open(unit=50, file='err.dat')
  open(unit=60, file='ers.dat')
  open(unit=70, file='erc.dat')
  do i=1,imax
     do j=1,jmax
        ! write seg image
        write(30,'(3i5)') i,j,seg(i,j)
        
        ! read in correct image & compare
        read(40,'(3i5)') idummy,jdummy,cor(i,j)
        correct(cor(i,j)) = correct(cor(i,j)) + 1 ! generate count of correct tissue types
        
        ! generate spatial error matrix & count error
        if ( cor(i,j) .ne. seg(i,j) ) then
           err(i,j) = 1 ! generate spatial error matrix
           ers(i,j) = seg(i,j) ! wrong segmented pixels types...
           erc(i,j) = cor(i,j) ! ...should be pixels types
           errata(seg(i,j)) = errata(seg(i,j)) + 1 ! count error based on type
        endif
      
        ! write error image & type specific error images
        write(50,'(3i5)') i,j,err(i,j)
        write(60,'(3i5)') i,j,ers(i,j)
        write(70,'(3i5)') i,j,erc(i,j)
     end do ! j
     write(30,*)
     read(40,*)
     write(50,*)
     write(60,*)
     write(70,*)
  end do ! i
  close(unit=30)
  close(unit=40)
  close(unit=50)
  close(unit=60)
  close(unit=70) 
  
  ! print count & percentage of correct designations
  errtot = sum(err)
  perc = errtot/dble(imax*jmax)
  
  write(*,*) 'runs: ',c, 'total flips: ',accept, 'errtot: ',errtot, 'err %:',perc
  write(*,*) 'error of tissue types (type, total, errata, err %):'
  write(*,*) ( i, correct(i), errata(i), errata(i)/dble(correct(i)), new_line('A'), i=1,5 )
  
end program samr

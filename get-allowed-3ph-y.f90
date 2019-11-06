program get_allowed
!
! A program to calculate the allowed 3-phonon interactions given
! a set of omega,q (read from a file). The code checks for transitions
! that conserve energy and quasi-momentum. For each mode lambda (where
! lambda combines q point and mode frequency), all pairs of modes that 
! correspond to phonon decays or collisions are found.
!
! Following Togo et al PRB 91 094306 (2015), this code aims to find
! the allowed transitions that contribute to the phonon linewidth, 
! Eq. 11. In this way, it is possible to see what scattering processes
! should be important in a material without doing the full 3-phonon
! force constants calculation.
!
! *******************
! This version differs in that the mesh.yaml file is read in directly
! to get the q points and frequencies, rather than having it converted
! first using the phonopy-mesh-convert-yaml.sh script
! *******************
!
! For each i, all pairs j,k are found that consist of allowed transitions.
! The output is:
! (1) a 3D set of coordinates of omega_i, omega_j, omega_k (each point
! is an allowed transition)
! (2) The 2 phonon joint DOS (Eq. 19 of Togo's paper)
! (3) The 2 phonon JDOS with occupations (Eqs. 23 and 24 of the paper)
!     For this calc a temperature is needed and is requested from the
!     user.
!
! J. Buckeridge July 2019
!
  implicit none
!
  integer, parameter :: nmax = 1e5
  integer, parameter :: nsmaller = 1e3
  integer, parameter :: maxang = 360
  real*8, parameter :: tiny = 1.d-12
  real*8, parameter :: big = 1.d12
  real*8, parameter :: cmtomev = 0.12398d0
  real*8, parameter :: small = 1.d-5
  real*8, parameter :: kboltz = 8.617343d-5 ! Boltzmann constant (eV/K)
  integer i, j, k, inew1, inew2, jnew1, jnew2, m, n, p, q, r
  integer stat, nqpts, natom, nmode, nlambda, countd, countc
  integer ngamma, nzero, idist1, idist2
  integer ilbrack, icomma, len1, mesh(3), meshqpts
  integer jmin, kmin, jlammin, nomin, omin(nmax,nsmaller)
  integer nequiv(nmax), totequiv, equivqpt(nmax,nmax)
  integer nunique
  integer iqchk, ijval
  real*8 pi, factor, lattvec(3,3), reclatt(3,3), vol
  real*8 qpt(3,nmax), qpt_c(3,nmax), omega(nmax,nsmaller), gvec(3), distg(nmax)
  real*8 omegamax, omegamin, distmin, dist, maxqcpt, enertest, qptest
  real*8 etestc1, etestc2, etestd
  real*8 tmpom1(3), tmpom2(3)
  character*80 str1
  character*1 lbrack, comma, colon
  logical checkg, checkq(3), checkqtot
!
  pi = acos(-1.d0)
!
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
  write(*,'(a)') "     333      PPP  H" 
  write(*,'(a)') "        3     P  P H" 
  write(*,'(a)') "       3   -- PPP  HHH   OOO  NNN   OOO   NNN " 
  write(*,'(a)') "        3     P    H  H O   O N  N O   O  N  N" 
  write(*,'(a)') "     333      P    H  H  OOO  N  N  OOO   N  N"
  write(*,'(a)')
  write(*,'(a)') "    TTTTT                  SSS     T                SSS" 
  write(*,'(a)') "      T                   S        T               S"
  write(*,'(a)') "      T   RRR   AAA  NNN   SSS  I TTT I  OOO  NNN   SSS"
  write(*,'(a)') "      T   R    A   A N  N     S I  T  I O   O N  N     S"
  write(*,'(a)') "      T   R     AAAA N  N  SSS  I  TT I  OOO  N  N  SSS" 
  write(*,'(a)') 
  write(*,'(a)') 
  write(*,'(a)') "j.buckeridge@ucl.ac.uk 2019"
  write(*,'(a)') "**************************************************************"
  write(*,'(a)') 
!
  write(*,*)
  write(*,'(a)') "Opening mesh.yaml to read input..."
  open(unit=11,file="mesh.yaml",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "ERROR: cannot find mesh.yaml file!!"
     stop
  else
     write(*,'(a)') "Found mesh.yaml. Checking q mesh..."
  endif
!
! To read in data from yaml, need to extract substrings and remove trailing
! commas. For this we need to define two useful substrings:
!
  lbrack = "["
  comma = ","
  colon = ":"
!
! Read in mesh dimensions - need to separate out strings in order to do
! this and avoid trailing commas
!
  read(11,'(a)',iostat=stat) str1
  if(stat /= 0) call error_read("mesh  ",0,0,stat)
  len1 = len(trim(str1))
  ilbrack = index(str1,lbrack) + 1
  str1 = str1(ilbrack:len1)
  icomma = index(str1,comma) - 1
  read(str1(1:icomma),*,iostat=stat) mesh(1)
  if(stat /= 0) call error_read("mesh 1",0,0,stat)
  len1 = len(trim(str1))
  ilbrack = icomma + 2
  str1 = str1(ilbrack:len1)
  icomma = index(str1,comma) - 1
  read(str1(1:icomma),*,iostat=stat) mesh(2)
  if(stat /= 0) call error_read("mesh 2",0,0,stat)
  len1 = len(trim(str1)) - 1
  ilbrack = icomma + 2
  str1 = str1(ilbrack:len1)
  read(str1,*,iostat=stat) mesh(3)
  if(stat /= 0) call error_read("mesh 3",0,0,stat)
  write(*,'(a16,3(x,i3))') "Mesh dimensions:", (mesh(i),i=1,3)
!
! Now read in number of qpts from file
!
  read(11,'(a)') str1
  if(stat /= 0) call error_read("qpts  ",0,0,stat)
  ilbrack = index(str1,colon) + 1
  len1 = len(trim(str1))
  str1 = str1(ilbrack:len1)
  read(str1,*,iostat=stat) nqpts
  if(stat /= 0) call error_read("nqpts ",0,0,stat)
!
! Check that the mesh dimensions and number of qpts are compatible.
! If there are fewer qpts, then symmetry has been used, which this code
! cannot deal with. If the number of qpts is too large, f knows what 
! happened but it's bad
!
  meshqpts = 1
  do i=1,3
     meshqpts = meshqpts * mesh(i)
  enddo
  if(nqpts < meshqpts) then
     write(*,*)
     write(*,'(a)') "ERROR: symmetry has been used to reduce number of qpts!!"
     write(*,'(a)') "       I can't deal with this situation - bye!!"
     stop
  endif
  if(nqpts > meshqpts) then
     write(*,*)
     write(*,'(a)') "ERROR: number of qpts is greater than what should be in the mesh!!"
     write(*,'(a)') "       this is too weird - bye!!"
     stop
  endif
!
! Next read in reciprocal lattice vectors. These must be multiplied 
! by 2pi to get the correct units
!
  write(*,'(a)') "Reading reciprocal lattice vectors..."
  read(11,*,iostat=stat)
  if(stat /= 0) call error_read("a line",0,0,stat)
  do i=1,3
     read(11,'(a)',iostat=stat) str1
     if(stat /= 0) call error_read("rlatt ",0,0,stat)
     len1 = len(trim(str1))
     ilbrack = index(str1,lbrack) + 1
     str1 = str1(ilbrack:len1)
     icomma = index(str1,comma) - 1
     read(str1(1:icomma),*,iostat=stat) reclatt(i,1)
     if(stat /= 0) call error_read("rlatt1",i,0,stat)
     len1 = len(trim(str1))
     ilbrack = icomma + 2
     str1 = str1(ilbrack:len1)
     icomma = index(str1,comma) - 1
     read(str1(1:icomma),*,iostat=stat) reclatt(i,2)
     if(stat /= 0) call error_read("rlatt2",i,0,stat)
     len1 = len(trim(str1)) - 6
     ilbrack = icomma + 2
     str1 = str1(ilbrack:len1)
     read(str1,*,iostat=stat) reclatt(i,3)
     if(stat /= 0) call error_read("rlatt3",i,0,stat)
  enddo
!
  reclatt(:,:) = reclatt(:,:) * 2.d0 * pi
!
  write(*,'(a)') "...done"
!
! Now read in the number of atoms N. There will be 3N modes at each qpt
!
  write(*,'(a)') "Checking number of atoms..."
  read(11,'(a)',iostat=stat) str1
  if(stat /= 0) call error_read("natom ",0,0,stat)
  ilbrack = index(str1,colon) + 1
  len1 = len(trim(str1))
  str1 = str1(ilbrack:len1)
  read(str1,*,iostat=stat) natom
  if(stat /= 0) call error_read("natom ",0,0,stat)
!
! Report the number of atoms
!
  write(*,'(a31,x,i3)') "Number of atoms in your system:", natom
!
! Set the number of modes = 3 * natom
!
  nmode = 3 * natom
  nlambda = nqpts * nmode
  write(*,'(a9,x,i8,x,a16)') "There are", nlambda, "(q,omega) points"
!
! Skip several lines in the yaml file where the lattice type and info on 
! the atoms is given. We don't need this info
!
  do i=1,5
     read(11,*,iostat=stat)
     if(stat /= 0) call error_read("a line",0,0,stat)
  enddo
  do i=1,natom
     do j=1,3
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("a line",0,0,stat)
     enddo
  enddo
!
! Look for "phonon" in the next line as sometimes there are blank lines
! before this one
!
  do
     read(11,'(a)',iostat=stat) str1
     if(stat /= 0) call error_read("a line",0,0,stat)
     if(str1(1:2) == "ph") exit
  enddo
!
! Now begin reading in qpt coordinates and phonon frequencies
!
  write(*,*)
  write(*,'(a)') "Reading in list of q-points and frequencies ..."
  write(*,'(a)') "WARNING: frequencies assumed to be in cm-1!!"
  do i=1,nqpts
!
! First read in qpt coordinates
!
     read(11,'(a)',iostat=stat) str1
     if(stat /= 0) call error_read("qptsl ",i,0,stat)
     len1 = len(trim(str1))
     ilbrack = index(str1,lbrack) + 1
     str1 = str1(ilbrack:len1)
     icomma = index(str1,comma) - 1
     read(str1(1:icomma),*,iostat=stat) qpt(1,i)
     if(stat /= 0) call error_read("qpts1 ",i,0,stat)
     len1 = len(trim(str1))
     ilbrack = icomma + 2
     str1 = str1(ilbrack:len1)
     icomma = index(str1,comma) - 1
     read(str1(1:icomma),*,iostat=stat) qpt(2,i)
     if(stat /= 0) call error_read("qpts2 ",i,0,stat)
     len1 = len(trim(str1)) - 1
     ilbrack = icomma + 2
     str1 = str1(ilbrack:len1)
     read(str1,*,iostat=stat) qpt(3,i)
     if(stat /= 0) call error_read("qpts3 ",i,0,stat)
!
! Convert qpts from fractional to cartesian
!
     do j=1,3
        qpt_c(j,i) = 0.d0
        do k=1,3
           qpt_c(j,i) = qpt_c(j,i) + qpt(k,i) * reclatt(k,j)
        enddo
     enddo
!
! Read in distance from Gamma for each qpt
!
     read(11,'(a)',iostat=stat) str1
     if(stat /= 0) call error_read("distl ",i,0,stat)
     ilbrack = index(str1,colon) + 1
     len1 = len(trim(str1))
     str1 = str1(ilbrack:len1)
     read(str1,*,iostat=stat) distg(i)
     if(stat /= 0) call error_read("dist  ",i,0,stat)
!
! Now read in frequencies. Skip a couple of lines first
!
     do j=1,2
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("a line",i,0,stat)
     enddo
     do j=1,nmode
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("a line",i,0,stat)
        read(11,'(a)',iostat=stat) str1
        if(stat /= 0) call error_read("omegal",i,j,stat)
        ilbrack = index(str1,colon) + 1
        len1 = len(trim(str1))
        str1 = str1(ilbrack:len1)
        read(str1,*,iostat=stat) omega(i,j)
        if(stat /= 0) call error_read("omega ",i,j,stat)
!
! Convert to eV
!
        omega(i,j) = omega(i,j) * cmtomev * 1.d-3
     enddo
     if( i < nqpts) then
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("a line",i,0,stat)
     endif
  enddo
!
  close(unit=11)
!
! Now we need to analyse the qpts list - we need to determine what to use as our
! criteria to decide when energies and wave vectors are close enough together so
! that the delta fns in Eq. 11 (Togo) are equal to one. For wave vectors, we find 
! the minimum distance between all pairs and take half that distance as the criterion.
!
! For energies it's a bit trickier. Ideally we would use the value of the three
! acoustic modes at Gamma, but Gamma may not be in our list of qpts. Also, if there
! are imaginary modes at Gamma, it's very difficult to find the zero acoustic modes.
! Also also, those modes need to be set to zero so they won't be included in the
! transitions (as they would have infinite occupation)
!
! What we will do is find the closest qpt to Gamma, find the lowest non-negative
! frequency there, and set all modes with freqs below this at Gamma to zero. We will
! then find the maximum frequency, and set the criterion as small * this frequency  
!
  omegamax = 0.d0
  omegamin = big
  ngamma = 0
  write(*,*)
  write(*,'(a)') "Checking data to establish values for energy and momentum checks..."
!
! Scan through all (q,omega) and check if Gamma is present, find the max omega and
! find the min omega at points away from Gamma
!
  do i=1,nqpts
     if( abs(distg(i)) < tiny ) then
        ngamma = i
     else
        do j=1,nmode
           if(omega(i,j) > tiny .and. omega(i,j) < omegamin) omegamin = omega(i,j)
        enddo
     endif
     do j=1,nmode
        if(omega(i,j) > omegamax) omegamax = omega(i,j)
     enddo
  enddo
  write(*,'(a45,x,f13.7,x,a4)') "Minimum frequency in range (away from Gamma):", &
       omegamin / (cmtomev * 1.d-3), "cm-1"
  write(*,'(a45,x,f13.7,x,a4)') "Maximum frequency in range                  :", &
       omegamax / (cmtomev * 1.d-3), "cm-1"
  write(*,'(a)') "Energy differences less than 10^-5 times this max frequency will be deemed to be zero..."
  write(*,*)
!
! If Gamma is present, set all frequencies below omegamin to zero. Report on this,
! and if any more modes beyond the three acoustic were set to zero
!
  nzero = 0
  if(ngamma > 0) then
     write(*,'(a33,x,i6)') "Found Gamma point in list at qpt:", ngamma 
     write(*,'(a)') "We need to set acoustic modes here to zero..."
     do i=1,nmode
        if(omega(ngamma,i) < omegamin) then
           omega(ngamma,i) = 0.d0
           nzero = nzero + 1
        endif
     enddo
     if(nzero == 0) then
        write(*,'(a)') "WARNING: Found very unusual dispersion at Gamma. No modes with lower"
        write(*,'(a)') "         frequency than lowest (non-negative) mode at next closest q pt!!"
        write(*,'(a)') "         I will continue with the calculation but something is up with"
        write(*,'(a)') "         your Gamma point frequencies..."
     else
        write(*,'(a11,x,i3,x,a32)') "Setting the", nzero, "lowest modes to zero at Gamma..."
        if(nzero > 3) write(*,'(a23,x,i3,x,a25)') "(i.e. you have at least", nzero - 3, &
             "imaginary modes at Gamma)"
     endif
  endif
!
! Now check for imaginary modes and set them to zero so they will not be included
! in the check for transitions
!
  write(*,'(a)') "Now checking for imaginary modes..."
  k = 0
  do i=1,nqpts
     do j=1,nmode
        if(omega(i,j) < 0.d0) then
           omega(i,j) = 0.d0
           k = k+1
           write(*,'(a16,x,i6,x,i6,x,a45)') "(q,omega) number", i, j, &
                "has negative frequency and will be omitted..."
        endif
     enddo
  enddo
  write(*,'(a5,x,i6,x,a18)') "Found", k, "imaginary modes..."
  write(*,'(a)') "...done"
!
! Next look for minimum distance between qpts. When it is found, get the maximum
! difference between components of the two qpts involved. Half of this will be used
! as the check to see if momentum is conserved
!
  write(*,*)
  write(*,'(a)') "Finding minimum distance between qpts..."
  distmin = 1.d6
  do i=1,nqpts
     do j=1,nqpts
        if(abs(qpt_c(1,i) - qpt_c(1,j)) < tiny .and. abs(qpt_c(2,i) - qpt_c(2,j)) < tiny &
             .and. abs(qpt_c(3,i) - qpt_c(3,j)) < tiny) cycle
        dist = 0
        do k=1,3
           dist = dist + ( qpt_c(k,i) - qpt_c(k,j) )**2
        enddo
        dist = sqrt(dist)
        if(dist < distmin) then
           idist1 = i
           idist2 = j
           distmin = dist
        endif
     enddo
  enddo
  maxqcpt = 0.d0
  do i=1,3
     if(abs(qpt_c(i,idist1) - qpt_c(i,idist2)) > maxqcpt) then
        maxqcpt = abs(qpt_c(i,idist1) - qpt_c(i,idist2))
     endif
  enddo
  write(*,'(a35,x,i6,x,a3,x,i6,x,a2,xf13.7,x,a3)') "Found minimum distance between qpts", &
       idist1, "and", idist2, "of", distmin, "A-1"
  write(*,'(a)') "Any 2 qpts where the differences between all components are less than half that"
  write(*,'(a)') "of these two will be deemed to be equal..."
!
! Set the criteria for energy and momentum
!
  enertest = omegamax * small
  qptest = 0.5d0 * maxqcpt
!
  write(*,*)
  write(*,'(a)') "Now checking for allowed transitions..."
!
! Open scratch file for output
!
  open(unit=12,status="scratch")
!
  countd = 0
  countc = 0
!
! First check q vectors for quasimomentum conservation. Calculate G vector
! and then check if q - q' - q" = G. We only check for decays as, if energy
! and momentum conservation is ok for a decay, then it is ok for a collision
! too. Check each component first, then check if all are true together.
! iqchk keeps track of how many cases of quasimomentum conservation are found
!
  iqchk = 0
  do i=1,nqpts
     ijval = 0
     do j=1,nqpts
        ijval = ijval+1
        inner_qloop: do k=ijval,nqpts
!
! We sum k from ijval in order to avoid checking cases where j,k = k',j', where
! j',k' is a previously checked case. In this way, we avoid double checking, as
! switching q' and q" in q - q' - q" = G makes no difference. We could also check
! for cases where q'=-q" etc (see labbook 04-11-19), but that would be expensive
! and wouldn't actually save time
!
!
! Initialise checks for q vector components
!
           checkqtot = .false.
           checkq(:) = .false.
!
! Search through possible G vectors (max value of m or n or p is 3). When a
! viable G vector is found, exit the search to save time. Keep count of number
! of cases found where q vector conservation applies
!
           outer_gloop: do m=-3,3
              do n=-3,3
                 do p=-3,3
                    do q=1,3
                       gvec(q) = real(m) * reclatt(1,q) + real(n) * reclatt(2,q)&
                            + real(p) * reclatt(3,q)
!
! Check for momentum conservation for decay processes
!
                       if(abs(qpt_c(q,i) - qpt_c(q,j) - qpt_c(q,k) - gvec(q)) <= qptest) then
                          checkq(q) = .true.
                       endif
                    enddo
                    if(checkq(1) .and. checkq(2) .and. checkq(3)) then
                       checkqtot = .true.
                       iqchk = iqchk + 1
                       exit outer_gloop
                    endif
                 enddo
              enddo
           enddo outer_gloop
!
! If momentum is conserved for decays, check energy conservation:
! omega - omega' - omega" = 0
! Exclude modes with zero frequency from these checks. Also exclude
! cases where omega < omega' or omega", as if this is true then the
! energy conservation condition cannot hold
!
           if(checkqtot) then
              do m=1,nmode
                 if(abs(omega(i,m)) <= tiny) cycle
                 do n=1,nmode
                    if(abs(omega(j,n)) <= tiny .or. omega(i,m) < omega(j,n)) cycle
                    do p=1,nmode
                       if(abs(omega(k,p)) <= tiny .or. omega(i,m) < omega(k,p)) cycle
                       etestd = omega(i,m) - omega(j,n) - omega(k,p)
                       if ( abs(etestd) <= enertest ) then
!
! Now we have an allowed phonon decay. Compute relevant occupations 
! and write frequencies in decay to scratch. Also write to scratch the
! allowed collision process q' + q" - q = G (omega' + omega" - omega = 0)
!
                          write(12,'(3(e21.13,3x))') omega(i,m) / (cmtomev * 1.d-3), &
                               omega(j,n) / (cmtomev * 1.d-3), omega(k,p) / (cmtomev * 1.d-3)
                          countd = countd + 1
                          countc = countc + 1
                       endif
                    enddo
                 enddo
              enddo
           endif
!
! Close loops for qpt checks
!
        enddo inner_qloop
     enddo
  enddo
!
! Say how many transitions were found
!
  write(*,'(a28,x,i12)') "Total number of conserved q:", iqchk
  write(*,'(a28,x,i12)') "Total number of transitions:", countd + countc
  write(*,'(a18,x,i12)') "Phonon decays    :", countd
  write(*,'(a18,x,i12)') "Phonon collisions:", countc 
  write(*,'(a)') "...done"
!
! Now remove redundancies from set of transitions. Open file to write 
! results in, and scan through scratch file to find unique transitions
!
  write(*,*)
  write(*,'(a)') "Removing non-unique transitions from output for plotting.."
!
  open(unit=11,file="allowed-omega.dat",status="replace",iostat=stat)
  nunique = 0
  rewind(unit=12)
  do i=1,countd
     read(12,*) (tmpom1(j),j=1,3)
!
! Change units so comparisons below are ok
!
     tmpom1(:) = tmpom1(:) * (cmtomev * 1.d-3)
     checkqtot = .true.
!
! Scan through values already written to output and check if current omega
! read from scratch is already contained there. Remember that points with 
! omega_2 and omega_3 swapped are equivalent and don't need to be outputted
! separately 
!
     if(nunique > 0) then
        do j=1,nunique
           read(11,*) (tmpom2(k),k=1,3)
!
! Change units so comparisons below are ok
!
           tmpom2(:) = tmpom2(:) * (cmtomev * 1.d-3)
           checkq(:) = .false.
!
! Check if all omega values are equal
!
           do k=1,3
              if(abs(tmpom1(k) - tmpom2(k)) < enertest) checkq(k) = .true.
           enddo
!
! Check if omega_2 and omega_3 are the same but just swapped
!
           if((checkq(2) .eqv. .false.) .and. (checkq(3) .eqv. .false.)) then
              if(abs(tmpom1(2) - tmpom2(3)) < enertest .and. abs(tmpom1(3) - tmpom2(2)) < enertest) then
                 checkq(2) = .true.
                 checkq(3) = .true.
              endif
           endif
           if(checkq(1) .and. checkq(2) .and. checkq(3)) checkqtot = .false.
        enddo
     endif
!
! If they aren't all equivalent, write to output (for the first point in the
! list, write straight to output). To include collisions in this output, 
! permute omega_1 and omega_2
!
     if(checkqtot) then
!
! Change units back
!
        write(11,'(3(e21.13,3x))') (tmpom1(j) / (cmtomev * 1.d-3) ,j=1,3)
        nunique = nunique+1
     endif
     rewind(unit=11)
  enddo
!
! Tell user how many unique transitions written to file
!
  write(*,'(a37,x,i7)') "Number of unique transitions to plot:", nunique
  write(*,'(a)') "...done"
!
  do i=11,12
     close(unit=i)
  enddo
!
end program get_allowed
!
! Subroutine to check for errors when reading in data
!
subroutine error_read(name, num1, num2,stat)
  implicit none
  character*6 name
  integer num1, num2, stat
  if(stat == 0) return
  if(stat > 0) then
     if( num1 == 0 .and. num2 == 0 ) then
        write(*,'(a22,x,a6)') "ERROR: failed to read:", name
        stop
     elseif( num2 == 0 ) then
        write(*,'(a22,x,a6,x,i5)') "ERROR: failed to read:", name, num1
        stop
     else
        write(*,'(a22,x,a6,x,i5)') "ERROR: failed to read:", name, num1, num2
        stop
     endif
  else
     if( num1 == 0 .and. num2 == 0 ) then
        write(*,'(a31,x,a6)') "ERROR: EOF when trying to read:", name
        stop
     elseif( num2 == 0 ) then
        write(*,'(a31,x,a6,x,i5)') "ERROR: EOF when trying to read:", name, num1
        stop
     else
        write(*,'(a31,x,a6,x,i5)') "ERROR: EOF when trying to read:", name, num1, num2
        stop
     endif
  endif
end subroutine error_read

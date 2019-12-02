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
! Moreover, if a phonon band structure calculation has been performed,
! so that the file band.yaml exists, the code checks if any of the 
! transitions correspond to (q,omega) points present in the band calc.
! If so, a file containing lines that can be plotted to show the allowed
! transitions is produced
!
! J. Buckeridge November 2019
!
  implicit none
!
  integer, parameter :: nmax = 1e5
  integer, parameter :: nsmaller = 1e3
  integer, parameter :: maxang = 360
  real*8, parameter :: tiny = 1.d-12
  real*8, parameter :: big = 1.d12
  real*8, parameter :: small = 1.d-5
  real*8, parameter :: kboltz = 8.617343d-5 ! Boltzmann constant (eV/K)
  integer i, j, k, inew1, inew2, jnew1, jnew2, m, n, p, q
  integer stat, nqpts, natom, nmode, nlambda, countc, countd
  integer ngamma, nzero, idist1, idist2
  integer ilbrack, icomma, len1, mesh(3), meshqpts
  integer jmin, kmin, jlammin, nomin
  integer bznqpts, bznatom, nbzgamma, totequiv
  integer memtest, nunique, ibzequiv
  integer iqchk, ijval, tmpqpt(3)
  integer, allocatable :: bzngamma(:), omin(:,:), nequiv(:), equivqpt(:,:)
  integer, allocatable :: idegen(:)
  real*8 pi, scale, factor, lattvec(3,3), reclatt(3,3), gvec(3), vol
  real*8 brlatt(3,3)
  real*8 omegamax, omegamin, distmin, dist, maxqcpt, enertest, qptest
  real*8 etestc1, etestc2, etestd
  real*8 temperature, occ1, occ2, occ3
  real*8 tmpom1(3), tmpom2(3)
  real*8, allocatable :: qpt(:,:), qpt_c(:,:), omega(:,:), distg(:)
  real*8, allocatable :: bzqpt(:,:), bzqpt_c(:,:), bzomega(:,:), bzdistg(:)
  real*8, allocatable :: d21dos(:,:), d22dos(:,:)
  real*8, allocatable :: n21dos(:,:), n22dos(:,:)
  character*80 str1
  character*5 unit
  character*1 lbrack, comma, colon
  logical checkg, checkq(3), checkqtot, checkomin, checkjdos
  logical checkbz, checkrec, checkgam, checkequiv
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
! First check for unit.dat file in which to read in the chosen units.
! If it is not present, default to cm-1
!
  write(*,*)
  write(*,'(a)') "Checking for unit.dat file to read in units..."
  open(unit=11,file="unit.dat",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "WARNING: unit.dat not found. Assuming cm-1..."
     scale = 1.d0 / 8065.54429d0
  else
!
! Found the file. Change the scale variable accordingly
!
     write(*,'(a)') "Found unit.dat file - attempting to read units..."
     do
        read(11,'(a)') str1
        if(str1(1:1) == "#") cycle
        read(str1,*,iostat=stat) unit
        if(unit(1:1) == "c" .or. unit(1:1) == "C") then
           write(*,'(a)') "Found units of cm-1..."
           scale = 1.d0 / 8065.54429d0
        elseif(unit(1:1) == "m" .or. unit(1:1) == "M") then
           write(*,'(a)') "Found units of meV..."
           scale = 1.d-3
        elseif(unit(1:1) == "t" .or. unit(1:1) == "T") then
           write(*,'(a)') "Found units of THz..."
           scale = 1.d0 / 241.8d0
        elseif(unit(1:1) == "e" .or. unit(1:1) == "E") then
           write(*,'(a)') "Found units of eV..."
           scale = 1.d0
        else
           write(*,'(a)') "ERROR: failed to read in units from unit.dat!"
           write(*,'(a)') "       I will assume cm-1..."
           scale = 1.d0 / 8065.54429d0
        endif
        exit
     enddo
     close(unit=11)
  endif
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
! If there are less qpts, then symmetry has been used, which this code
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
!
! Allocate memory to store the qpt and omega arrays
!
  allocate( qpt(3,nqpts), stat=memtest)
  call mem_test("qpt    ",memtest)
  allocate( qpt_c(3,nqpts), stat=memtest)
  call mem_test("qpt_c  ",memtest)
  allocate( distg(nqpts), stat=memtest)
  call mem_test("distg  ",memtest)
  allocate( omega(nqpts,nmode), stat=memtest)
  call mem_test("omega  ",memtest)
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
        omega(i,j) = omega(i,j) * scale
     enddo
     if( i < nqpts) then
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("a line",i,0,stat)
     endif
  enddo
!
  close(unit=11)
!
! Check for presence of temperature.dat file. If it is present, then
! the 2 phonon JDOS will be written to file
!
  checkjdos = .false.
  open(unit=15,file="temperature.dat",status="old",iostat=stat)
  if(stat == 0) then
     write(*,'(a)') "Found temperature.dat file - attempting to read temperature..."
     do
        read(15,'(a)') str1
        if(str1(1:1) == "#") cycle
        read(str1,*,iostat=stat) temperature
        if(stat /= 0) then
           write(*,'(a)') "WARNING: failed to read in temperature correctly. Skipping JDOS calculation..."
           exit
        else     
           write(*,'(a19,x,f13.6,x,a1)') "Found temperature =", temperature, "K"
           write(*,'(a)') "JDOS will be calculated and written to file jdos.dat"
           write(*,'(a)') "NOTE: JDOS calculated here should be scaled by 1/N"
           write(*,*)         
           checkjdos = .true.
           exit
        endif
     enddo
     close(unit=15)
  else
     write(*,'(a)') "No temperature.dat file present - JDOS not calculated"
     write(*,*)
  endif
!
! We can now deallocate qpt as it is no longer needed, unless the JDOS is 
! calculated, as we will then need it for the output
!
  if(checkjdos .eqv. .false.) then
     deallocate(qpt,stat=memtest)
     call dealloc_test("qpt    ",memtest)
  endif
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
  write(*,'(a45,x,f13.7,x,a5)') "Minimum frequency in range (away from Gamma):", &
       omegamin / (scale), unit
  write(*,'(a45,x,f13.7,x,a5)') "Maximum frequency in range                  :", &
       omegamax / (scale), unit
  write(*,'(a)') "Energy differences less than 10^-5 times this max frequency will be deemed to be zero..."
  write(*,*)
!
! We can now deallocate distg as it is no longer needed
!
  deallocate(distg,stat=memtest)
  call dealloc_test("distg  ",memtest)
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
  do i=1,nqpts
     do j=1,nmode
        if(omega(i,j) < 0.d0) then
           omega(i,j) = 0.d0
           write(*,'(a16,x,i6,x,i6,x,a45)') "(q,omega) number", i, j, &
                "has negative frequency and will be omitted..."
        endif
     enddo
  enddo
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
! Now open band.yaml to read in data
!
  checkbz = .false.
  write(*,*)
  write(*,'(a)') "Looking for band.yaml to read input..."
  open(unit=11,file="band.yaml",status="old",iostat=stat)
!
  if(stat /= 0) then
     write(*,'(a)') "...band.yaml not found..."
  else
     write(*,'(a)') "Found band.yaml. Reading in data..."
     checkbz = .true.
  endif
!
  if(checkbz) then
!
! Read in number of qpts in BZ path from file
!
     read(11,'(a)') str1
     if(stat /= 0) call error_read("bzqpts",0,0,stat)
     ilbrack = index(str1,colon) + 1
     len1 = len(trim(str1))
     str1 = str1(ilbrack:len1)
     read(str1,*,iostat=stat) bznqpts
     if(stat /= 0) call error_read("bznqpt",0,0,stat)
!
! Report the number of qpts in BZ path
!
     write(*,'(a26,x,i6)') "Number of qpts in BZ path:", bznqpts
!
! Next read in reciprocal lattice vectors for the bandstructure calc. 
! These must be multiplied by 2pi to get the correct units and should
! be the same as those used for the mesh calc, but if not we can still
! use the data. In some cases the reciprocal lattice may not be in the
! band.yaml file, so check for this
!
     write(*,'(a)') "Reading reciprocal lattice vectors from band.yaml..."
     i = 1
     do
        i = i+1
        read(11,'(a)',iostat=stat) str1
        if(stat < 0) call error_read("line  ",i,0,stat)
        if(str1(1:5) == "recip" .or. stat < 0) exit
     enddo
     if (stat >= 0) then
        do i=1,3
           read(11,'(a)',iostat=stat) str1
           if(stat /= 0) call error_read("blatt ",0,0,stat)
           len1 = len(trim(str1))
           ilbrack = index(str1,lbrack) + 1
           str1 = str1(ilbrack:len1)
           icomma = index(str1,comma) - 1
           read(str1(1:icomma),*,iostat=stat) brlatt(i,1)
           if(stat /= 0) call error_read("blatt1",i,0,stat)
           len1 = len(trim(str1))
           ilbrack = icomma + 2
           str1 = str1(ilbrack:len1)
           icomma = index(str1,comma) - 1
           read(str1(1:icomma),*,iostat=stat) brlatt(i,2)
           if(stat /= 0) call error_read("blatt2",i,0,stat)
           len1 = len(trim(str1)) - 6
           ilbrack = icomma + 2
           str1 = str1(ilbrack:len1)
           read(str1,*,iostat=stat) brlatt(i,3)
           if(stat /= 0) call error_read("blatt3",i,0,stat)
        enddo
!
        brlatt(:,:) = brlatt(:,:) * 2.d0 * pi
!
        write(*,'(a)') "...done"
!
! Check if the reciprocal lattice vectors are different to those in mesh.yaml
!
        checkrec = .false.
        do i=1,3
           do j=1,3
              if(abs(brlatt(i,j) - reclatt(i,j)) > small) checkrec = .true.
           enddo
        enddo
        if(checkrec) then
           write(*,'(a)') "WARNING: reciprocal lattices in band.yaml and mesh.yaml are different"
           write(*,*)
        endif
     else
        write(*,'(a)') "...not found. Using reciprocal lattice vectors from mesh.yaml..."
        brlatt(:,:) = reclatt(:,:)
        rewind(unit=11)
     endif
!
! Look for line starting with "natom" and read the number of atoms
!
     do
        read(11,'(a)',iostat=stat) str1
        if(stat /= 0) call error_read("natom ",0,0,stat)
        if(str1(1:5) == "natom") then       
           ilbrack = index(str1,colon) + 1
           len1 = len(trim(str1))
           str1 = str1(ilbrack:len1)
           read(str1,*,iostat=stat) bznatom
           if(stat /= 0) call error_read("natomb",0,0,stat)
           exit
        endif
     enddo
!
! Check that the number of atoms is the same as that in mesh.yaml. If not
! then complain to the user and quit
!
     if(natom /= bznatom) then
        write(*,*)
        write(*,'(a)') "ERROR: the number of atoms in band.yaml is incompatible with that in mesh.yaml!!"
        write(*,'(a)') "       I can't continue :("
        stop
     endif
!
! Look for line starting with "phonon"
!
     do
        read(11,'(a)',iostat=stat) str1
        if(stat /= 0) call error_read("phonon",0,0,stat)
        if(str1(1:2) == "ph") exit
     enddo
!
! Now begin reading in qpt coordinates and phonon frequencies
!
     write(*,*)
     write(*,'(a)') "Reading in BZ path q-points and frequencies ..."
!
! Allocate memory to store the qpt and omega arrays
!
     allocate( bzqpt(3,bznqpts), stat=memtest)
     call mem_test("bzqpt  ",memtest)
     allocate( bzqpt_c(3,bznqpts), stat=memtest)
     call mem_test("bzqpt_c",memtest)
     allocate( bzdistg(bznqpts), stat=memtest)
     call mem_test("bzdistg",memtest)
     allocate( bzomega(bznqpts,nmode), stat=memtest)
     call mem_test("bzomega",memtest)
     do i=1,bznqpts
!
! First read in qpt coordinates
!
        read(11,'(a)',iostat=stat) str1
        if(stat /= 0) call error_read("bqptsl",i,0,stat)
        len1 = len(trim(str1))
        ilbrack = index(str1,lbrack) + 1
        str1 = str1(ilbrack:len1)
        icomma = index(str1,comma) - 1
        read(str1(1:icomma),*,iostat=stat) bzqpt(1,i)
        if(stat /= 0) call error_read("bqpts1",i,0,stat)
        len1 = len(trim(str1))
        ilbrack = icomma + 2
        str1 = str1(ilbrack:len1)
        icomma = index(str1,comma) - 1
        read(str1(1:icomma),*,iostat=stat) bzqpt(2,i)
        if(stat /= 0) call error_read("bqpts2",i,0,stat)
        len1 = len(trim(str1)) - 1
        ilbrack = icomma + 2
        str1 = str1(ilbrack:len1)
        read(str1,*,iostat=stat) bzqpt(3,i)
        if(stat /= 0) call error_read("bqpts3",i,0,stat)
!
! Convert qpts from fractional to cartesian
!
        do j=1,3
           bzqpt_c(j,i) = 0.d0
           do k=1,3
              bzqpt_c(j,i) = bzqpt_c(j,i) + bzqpt(k,i) * brlatt(k,j)
           enddo
        enddo
!
! Read in distance from Gamma for each qpt
!
        read(11,'(a)',iostat=stat) str1
        if(stat /= 0) call error_read("bdistl",i,0,stat)
        ilbrack = index(str1,colon) + 1
        len1 = len(trim(str1))
        str1 = str1(ilbrack:len1)
        read(str1,*,iostat=stat) bzdistg(i)
        if(stat /= 0) call error_read("bdist ",i,0,stat)
!
! Now read in frequencies. Skip the "band" line first
!
        read(11,*,iostat=stat)
        if(stat /= 0) call error_read("band  ",i,0,stat)
        do j=1,nmode
           read(11,*,iostat=stat)
           if(stat /= 0) call error_read("bline ",i,0,stat)
           read(11,'(a)',iostat=stat) str1
           if(stat /= 0) call error_read("omegal",i,j,stat)
           ilbrack = index(str1,colon) + 1
           len1 = len(trim(str1))
           str1 = str1(ilbrack:len1)
           read(str1,*,iostat=stat) bzomega(i,j)
           if(stat /= 0) call error_read("omega ",i,j,stat)
!
! Convert to eV
!
           bzomega(i,j) = bzomega(i,j) * scale
        enddo
        if( i < bznqpts) then
           read(11,*,iostat=stat)
           if(stat /= 0) call error_read("bline ",i,0,stat)
        endif
     enddo
!
! We can now deallocate bzqpt as it is no longer needed
!
     deallocate(bzqpt,stat=memtest)
     call dealloc_test("bzqpt  ",memtest)
!
! Close file
!
     write(*,'(a)') "...done"
     close(unit=11)
!
! Check that both mesh and band contain the Gamma point. If so, then
! check if any frequencies are different - if they are, this will be
! due to the inclusion of the non-analytical correction for TO-LO 
! splitting. Then change the frequencies in the mesh gamma point to
! those in the band file (as these have the TO-LO splitting interpolated
! at Gamma, while the mesh values do not)
!
! First allocate memory for array to store location of gamma points in
! BZ list
!
     allocate( bzngamma(bznqpts), stat=memtest)
     call mem_test("bzngam ",memtest)
!
! Of course, if there are only 3 modes (i.e. one atom), then this check
! is not necessary, but we still need to save the locations of the Gamma
! points
!
     bzngamma = 0
     if(ngamma > 0) then
        write(*,*)
        write(*,'(a)') "Checking for Gamma point in BZ path in band.yaml..."
        nbzgamma = 0
!
! Now try to find all occurrence of Gamma in BZ list and save their locations
!
        do i=1,bznqpts
           checkq(:) = .false.
           do j=1,3
              if(abs(bzqpt_c(j,i)) <= small) checkq(j) = .true.
           enddo
           if(checkq(1) .and. checkq(2) .and. checkq(3)) then
              nbzgamma = nbzgamma + 1
              bzngamma(nbzgamma) = i
              exit
           endif
        enddo
        if(nbzgamma > 0) then
!
! Found it. Check for differences in frequency, but leave out the first few modes
! as these will probably be different even though they are meant to be zero or 
! imaginary. nzero was used above to determine how many modes at Gamma were set
! to zero. If all modes at Gamma were set to zero, which could be the case if the
! system has just one atom for example, then skip this check 
!
           write(*,'(a52,x,i6)') "Gamma point found in BZ list, 1st occurrence at qpt:", &
                bzngamma(1)
           checkgam = .false.
           if(nzero < nmode) then
              do i=nzero+1,nmode
                 if( abs(bzomega(bzngamma(1),i) - omega(ngamma,i)) > enertest ) then
                    checkgam = .true.
                    write(*,'(a66,x,i4)') "Different frequency at Gamma b/n band.yaml &
                         and mesh.yaml for mode:", i
!
! The frequencies that differ must be set the same
!
                    omega(ngamma,i) = bzomega(bzngamma(1),i)
                 endif
              enddo
           endif
           if(checkgam) then
              write(*,'(a)') "(I am assuming that non-analytical correction has been used!"
              write(*,'(a)') "Mesh calcs do not get the LO-TO splitting at Gamma - so we"
              write(*,'(a)') "will make the Gamma point frequencies for the mesh equal to"
              write(*,'(a)') "those of the band calculation...)"
           else
              write(*,'(a)') "Frequencies in band.yaml and mesh.yaml are the same at Gamma..."
           endif
        else
           write(*,'(a)') "...Gamma point not found..."
        endif
     endif
!
! We can now deallocate bzqpt_c as it is no longer needed
!
     deallocate(bzqpt_c,stat=memtest)
     call dealloc_test("bzqpt_c",memtest)
!
! Go through list of qpts in band.yaml and find which qpts in mesh they are equivalent
! to according to their frequencies. Watch out for Gamma, as the lowest freq modes may
! differ slightly but should still be considered the same
!
! First allocate the arrays
!
     allocate( nequiv(nqpts), stat=memtest)
     call mem_test("nequiv ",memtest)
     allocate( equivqpt(nqpts,nqpts), stat=memtest)
     call mem_test("equivqp",memtest)
!
     write(*,*)
     write(*,'(a)') "Checking for equivalent qpts between BZ list and mesh qpts..."
     write(*,'(a)') "(by comparing frequencies, watch out if you have BZ directions"
     write(*,'(a)') "that have zero dispersion for ALL bands)"
!
     equivqpt = 0
     totequiv = 0
     do i=1,nqpts
        nequiv(i) = 0
!
! First check if this is the Gamma point. If so, then assign it correctly
! if the Gamma point is also in the mesh
!
        if(i == ngamma) then
           nequiv(i) = nbzgamma
           if(nbzgamma > 0) then
              do j=1,nequiv(i)
                 equivqpt(i,j) = bzngamma(j)
              enddo
           endif
        else
!
! This is not the Gamma point, so check for equivalent qpts in full list
! by comparing frequencies
!  
           do j=1,bznqpts
              checkequiv = .true.
              do k=1,nmode
                 if( abs(omega(i,k) - bzomega(j,k)) > enertest ) then
                    checkequiv = .false.
                    exit
                 endif
              enddo
!
! If frequencies do not differ, or if we are comparing Gamma points, we have two
! equivalent qpts
!
              if(checkequiv) then
                 nequiv(i) = nequiv(i) + 1
                 equivqpt(i,nequiv(i)) = j
              endif
           enddo
        endif
        totequiv = totequiv + nequiv(i)
     enddo
     write(*,'(a49,x,i7)') "Total equivalences between BZ qpts and mesh qpts:", totequiv 
!
! We can now deallocate bzomega and bzngamma as they are no longer needed
!
     deallocate(bzomega,stat=memtest)
     call dealloc_test("bzomega",memtest)
     deallocate(bzngamma,stat=memtest)
     call dealloc_test("bzngam ",memtest)
!
! Close check for presence of band.yaml file
!
  endif
!
  write(*,*)
!
  write(*,'(a)') "Now checking for allowed transitions..."
!
! Open file for output
!
  open(unit=20,status="scratch")
  if(checkjdos) then
!
! Allocate arrays and set initial contributions to JDOS to zero
!
     allocate( d21dos(nqpts,nmode), stat=memtest)
     call mem_test("d21dos ",memtest)
     allocate( d22dos(nqpts,nmode), stat=memtest)
     call mem_test("d22dos ",memtest)
     d21dos = 0.d0
     d22dos = 0.d0
     allocate( n21dos(nqpts,nmode), stat=memtest)
     call mem_test("n21dos ",memtest)
     allocate( n22dos(nqpts,nmode), stat=memtest)
     call mem_test("n22dos ",memtest)
     n21dos = 0.d0
     n22dos = 0.d0
  endif
!
! First check q vectors for quasimomentum conservation. Calculate G vector
! and then check if q - q' - q" = G. If this is true, then so is the collision
! process q' + q" - q = -G. Check each component first, then check if all 
! are true together
!
  countd = 0
  countc = 0
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
           checkq(:) = .false.
           checkqtot = .false.
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
                          write(20,'(3(e21.13,3x),3(i10,3x))') omega(i,m) / (scale), &
                               omega(j,n) / (scale), omega(k,p) / (scale), &
                               i, j, k
                          countd = countd + 1
                          countc = countc + 1
!
! Add to summation for phonon JDOS (Eq. 19 of Togo's paper)
!
                          if(checkjdos) then
                             d22dos(i,m) = d22dos(i,m) + 1.d0
                             d21dos(j,n) = d21dos(j,n) + 1.d0
                             d21dos(k,p) = d21dos(k,p) + 1.d0
!
! Use occupations to compute N2 JDOS (Eq. 24 of Togo's paper)
!
                             occ1 = 1.d0 / (exp( omega(i,m) / (kboltz * temperature) ) - 1.d0)
                             occ2 = 1.d0 / (exp( omega(j,n) / (kboltz * temperature) ) - 1.d0)
                             occ3 = 1.d0 / (exp( omega(k,p) / (kboltz * temperature) ) - 1.d0)
                             n22dos(i,m) = n22dos(i,m) + occ2 + occ3 + 1.d0
                             n21dos(j,n) = n21dos(j,n) + abs(occ1 - occ3)
                             n21dos(k,p) = n21dos(k,p) + abs(occ1 - occ2)
                          endif
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
! Write output for the joint DOS
!
  if(checkjdos) then
     write(*,*) 
     write(*,'(a)') "Now writing JDOS to file (this could be a big one!)..."
     open(unit=12,file="jdos.dat",status="replace",iostat=stat)
     write(12,'(a1,x,3(a2,22x),a5,19x,2(a5,19x,2(a6,18x)))') "#", "qx", "qy", "qz", "omega", &
          &"d2dos", "d2dos1", "d2dos2", "n2dos", "n2dos1", "n2dos2"
     do i=1,nqpts
        do j=1,nmode
           write(12,'(10(e21.13,3x))') (qpt(k,i),k=1,3), omega(i,j) / (scale), &
                &d21dos(i,j) + d22dos(i,j), d21dos(i,j), d22dos(i,j), &
                &n21dos(i,j) + n22dos(i,j), n21dos(i,j), n22dos(i,j)
        enddo
     enddo
     close(unit=12)
     deallocate(qpt,stat=memtest)
     call dealloc_test("qpt    ",memtest)
     deallocate(d21dos,stat=memtest)
     call dealloc_test("d21dos ",memtest)
     deallocate(d22dos,stat=memtest)
     call dealloc_test("d22dos ",memtest)
     deallocate(n21dos,stat=memtest)
     call dealloc_test("n21dos ",memtest)
     deallocate(n22dos,stat=memtest)
     call dealloc_test("n22dos ",memtest)
     write(*,'(a)') "...done"
  endif
!
! Deallocate arrays not needed anymore
!
  deallocate(omega,stat=memtest)
  call dealloc_test("omega  ",memtest)
  deallocate(qpt_c,stat=memtest)
  call dealloc_test("qpt_c  ",memtest)
!
! Now remove redundancies from set of transitions. Open a new scratch file
! to write results in first, and scan through scratch file to find unique 
! transitions. The scratch file is used so that we can count the degeneracy
! of each unique transition. We write the output to file afterwards
!
  write(*,*)
  write(*,'(a)') "Removing non-unique transitions from output for plotting.."
!
  open(unit=21,status="scratch")
!
! Allocate array to store degeneracies
!
  allocate(idegen(countd),stat=memtest)
  call mem_test("idegen ",memtest)
!
! Inititialise variables
!
  nunique = 0
  ibzequiv = 0
  idegen = 0
!
  rewind(unit=20)
  do i=1,countd
     read(20,*) (tmpom1(j),j=1,3), (tmpqpt(j),j=1,3)
!
! Change units so comparisons below are ok
!
     tmpom1(:) = tmpom1(:) * (scale)
     checkqtot = .true.
!
! Scan through values already written to output and check if current omega
! read from scratch is already contained there. Remember that points with 
! omega_2 and omega_3 swapped are equivalent and don't need to be outputted
! separately 
!
     if(nunique > 0) then
        do j=1,nunique
           read(21,*) (tmpom2(k),k=1,3)
!
! Change units so comparisons below are ok
!                                                                                                                 
           tmpom2(:) = tmpom2(:) * (scale)
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
           if(checkq(1) .and. checkq(2) .and. checkq(3)) then
              checkqtot = .false.
              idegen(j) = idegen(j) + 1
           endif
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
        write(21,'(3(e21.13,3x))') (tmpom1(j) / (scale) ,j=1,3)
        nunique = nunique+1
        idegen(nunique) = 1
!
! Now check if the non-unique transition is contained within the band data.
! If so, create a plot file so that it can be plotted as two lines, starting
! at the decaying phonon's position in the band diagram and ending at the 
! positions of the resultant phonons. Look for multiple occurances, and if
! they arise, then plot the lines so that they join to the nearest points
! on the band diagram
!
        if(checkbz) then
           checkequiv = .true.
           do j=1,3
              if(nequiv(tmpqpt(j)) <= 0) checkequiv = .false.
           enddo
           if(checkequiv) then
              if(ibzequiv == 0) open(unit=14,file="bz-transitions.dat",status="replace",iostat=stat)
              ibzequiv = ibzequiv + 1
              do m=1,nequiv(tmpqpt(1))
                 do n=1,nequiv(tmpqpt(2))
                    do p=1,nequiv(tmpqpt(3))
                       write(14,'(2(e21.13,3x))') bzdistg(equivqpt(tmpqpt(2),n)), tmpom1(2) / (scale)
                       write(14,'(2(e21.13,3x))') bzdistg(equivqpt(tmpqpt(1),m)), tmpom1(1) / (scale)
                       write(14,*)
                       write(14,'(2(e21.13,3x))') bzdistg(equivqpt(tmpqpt(3),p)), tmpom1(3) / (scale)
                       write(14,'(2(e21.13,3x))') bzdistg(equivqpt(tmpqpt(1),m)), tmpom1(1) / (scale)
                       write(14,*)
                    enddo
                 enddo
              enddo
           endif
        endif
     endif
     rewind(unit=21)
  enddo
!
! Tell user how many unique transitions written to file
!
  write(*,'(a37,x,i7)') "Number of unique transitions to plot:", nunique
!
! Tell user how many transitions will be plotted in band plot and close file
!
  if(checkbz) then 
     write(*,'(a37,x,i7)') "Number of transitions in bz plot    :", ibzequiv
     if(ibzequiv > 0) close(unit=14)
  endif
!
! Now write to file the unique allowed transitions plus their degeneracies.
! The degeneracies are multiplied by 3, as each allowed transition corresponds
! to one decay and two collisions
!
  open(unit=11,file="allowed-omega.dat",status="replace",iostat=stat)
  write(11,'(a1,x,3(a6,18x),a10)') "#", "omega1", "omega2", "omega3", "degeneracy"
  do i=1,nunique
     read(21,*) (tmpom1(j),j=1,3)
     write(11,'(3(e21.13,3x),i10)') (tmpom1(j),j=1,3), idegen(i)
  enddo
!
! Deallocate degeneracy array
!
  deallocate(idegen,stat=memtest)
  call dealloc_test("idegen ",memtest)
!
  write(*,'(a)') "...done"
!
  close(unit=11)
  close(unit=20)
  close(unit=21)
!
! Deallocate remaining arrays from BZ calc
!
  if(checkbz) then
     deallocate(bzdistg,stat=memtest)
     call dealloc_test("bzdistg",memtest)
     deallocate(equivqpt,stat=memtest)
     call dealloc_test("equivqp",memtest)
     deallocate(nequiv,stat=memtest)
     call dealloc_test("nequiv ",memtest)
  endif
!
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
!
subroutine mem_test(name,stat)
  implicit none
  character*7 name
  integer stat
  if(stat == 0) then
     return
  else
     write(*,'(a21,x,a7)') "ERROR: OOM for array:", name
     stop
  endif
end subroutine mem_test
!
subroutine dealloc_test(name,stat)
  implicit none
  character*7 name
  integer stat
  if(stat == 0) then
     return
  else
     write(*,'(a28,x,a7)') "ERROR: failed to deallocate:", name
     stop
  endif
end subroutine dealloc_test

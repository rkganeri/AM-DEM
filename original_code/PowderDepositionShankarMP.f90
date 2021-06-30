program PowderDepositionShankarMP
use RDistributions
use omp_lib
implicit none

integer ::numparticles, i, j, k, N, Nlength, NWidth, NHeight, NTotal, counter, layer, material, numlayers, numbinsx, numbinsy, numbinsz, totalnumbins, totalnumparticles
double precision, dimension(:,:), allocatable :: CenterCoords, RandNumbers, Vnew, Vold, PositionOld, PositionNew, PsiTot1, PsiTot2, PsiTot3, PsiTot4, Y1, Y2, Y3, Y4, mean, stddev
double precision, dimension(:), allocatable :: Radius, Mass, Temp, volume
integer, dimension(:), allocatable :: Inside, linkedbinlist
double precision, dimension(3) :: RandNumVec, Vwall
double precision :: xdist, ydist, zdist, distance, length, width, height, r1, r2, r3, t, tstart, tend, tstep, finalheight, m1, m2, m3, rho1, rho2, rho3, Vol1, Vol2, Vol3, meanradius, stddevradius, BINLENGTH, binwidth, binheight, minRadius, maxRadius, rho, nu, Eyoung
integer, dimension(:,:,:), allocatable :: bins
integer, dimension(:,:), allocatable :: particlebin



!!OPEN TECPLOT OUTPUT FILE
Open(unit=100, file='PowderDepositionShankarMP2800.dat')
write(unit=100,fmt='(A)') 'Title="DEMParticles"'
write(unit=100,fmt='(A)') 'Variables="X","Y","Z","R","In","T"'
70 format(F13.8,2x,F13.8,2x,F13.8,2x,F13.8,2x,I10,2x,F13.8)

Open(unit=200, file='PowderDepositionShankarMPFinal2800.dat')

! Generate a set of particles and randomly place them in a box (as small as possible to produce 1 mm height)
length = 0.00040D0 ! [m]
width = 0.00080D0
height = 0.00350D0

meanradius = 24.5e-6 ! particle diameter evenly distributed from 45 - 53 microns
!stddevradius = 4.0e-6
minradius = 22.5e-6
maxradius = 26.5e-6
numlayers = 1
numparticles = 2800
totalnumparticles = numparticles*numlayers

allocate(CenterCoords(numparticles,3))
allocate(RandNumbers(numparticles,4))
allocate(Radius(numparticles))
allocate(Mass(numparticles))
allocate(Volume(numparticles))
allocate(Vnew(numparticles,3))
allocate(Vold(numparticles,3))
allocate(PositionOld(numparticles,3))
allocate(PositionNew(numparticles,3))
allocate(psitot1(numparticles,3))
allocate(psitot2(numparticles,3))
allocate(psitot3(numparticles,3))
allocate(psitot4(numparticles,3))
allocate(Y1(numparticles,6))
allocate(Y2(numparticles,6))
allocate(Y3(numparticles,6))
allocate(Y4(numparticles,6))
allocate(Inside(numparticles))
allocate(Temp(numparticles))
Inside(:) = 0
Temp(:) = 373.0


tstart = 0.0D0
tend = 0.08D0

t = tstart
tstep = 5.0e-8
Vold(:,:) = 0.0D0 ! Velocity of particles


!! 316L SS material properties at 373 K
rho = 7919.0D0 ! kg/m^3
nu = 0.260D0  ! Poisson's ratio
Eyoung = 193.0e9  ! Pa


!!$FinalHeight = 0.02*height
!!$Vwall(1) = 0.0D0
!!$Vwall(2) = 0.0D0
!!$Vwall(3) = (FinalHeight-height)/(0.8*tend-tstart)


Call init_Random_Seed()
Call Random_Number(RandNumbers)
!! Determine mean and std deviation of particle radii for Gaussian distribution
!mean = 2.5e-5
!stddev = 1.0e-5

! Run loop to sinter 6 layers (2 for material 1, 2 for material 2, 2 for material 3)
do layer = 1,1

   ! Generate random coordinates of particle centers (z-coord is constrained to top 1/3rd of box), and radius sizes according to Gaussian distribution
   do i = 1,numparticles 
!!$      Radius(i) = rand_normal(meanRadius,stddevRadius)
!!$      if (Radius(i) < minRadius) then  !Ensure that minimum particle radius is 10 microns and max is 50 microns (average is 25 microns, stdev is 10 microns)
!!$         Radius(i) = minRadius
!!$      elseif (Radius(i) > maxRadius) then
!!$         Radius(i) = maxRadius
!!$      end if

      Radius(i) = minRadius + (maxRadius-minRadius)*RandNumbers(i,4)

      Volume(i) = 4.0D0/3.0D0*PI*Radius(i)**3
      CenterCoords(i,1) = Radius(i) + (length-2*Radius(i))*RandNumbers(i,1)
      CenterCoords(i,2) = Radius(i) + (width-2*Radius(i))*RandNumbers(i,2)
      CenterCoords(i,3) = Height - Radius(i) - 0.8*height*RandNumbers(i,3)
   end do



   ! Ensuring no particles (spheres) overlap
   i = numparticles - 1
   do while (i > 0)
      N = i + 1
      do while (N <= numparticles)
         xdist = CenterCoords(i,1)-CenterCoords(N,1)
         ydist = CenterCoords(i,2)-CenterCoords(N,2)
         zdist = CenterCoords(i,3)-CenterCoords(N,3)
         distance = (xdist**2 + ydist**2 + zdist**2)**(0.5)

         if (distance<(Radius(i)+Radius(N))) then
            Call Random_Number(RandNumVec)
            CenterCoords(i,1) = Radius(i) + (length-2*Radius(i))*RandNumVec(1)
            CenterCoords(i,2) = Radius(i) + (width-2*Radius(i))*RandNumVec(2)
            CenterCoords(i,3) = Height -Radius(i) - 0.50*height*RandNumVec(3)
            N = i+1
         else
            N = N+1
         end if
      end do
      i = i - 1
   end do

!!$! Calculate mass of each particle, from 0-1 mm is material 1, 1-1.5 mm is mat 2, 1.5-3 mm is mat3
!!$   do i = 1,numparticles
!!$      if (CenterCoords(i,2) <= 0.001) then
!!$         Mass(i) = Rho1*Volume(i)
!!$      elseif ((CenterCoords(i,2) > 0.001) .AND. (CenterCoords(i,2) <= 0.0015)) then
!!$         Mass(i) = Rho2*Volume(i)
!!$      else
!!$         Mass(i) = Rho3*Volume(i)
!!$      end if
!!$   end do

! Define mass of each particle to be that of HX SLM material
   do i = 1,numparticles
      Mass(i) = Rho*Volume(i)
   end do
      

   PositionOld(:,:) = CenterCoords(:,:)

   ! Define number of bins
   numbinsx = floor(length/(4.0*maxval(Radius)))
   numbinsy = floor(width/(4.0*maxval(Radius)))
   numbinsz = floor(height/(4.0*maxval(Radius)))
   totalnumbins = numbinsx*numbinsy*numbinsz
   binlength = length/numbinsx
   binwidth = width/numbinsy
   binheight = height/numbinsz

   ! Check x, y, and z coordinates of each particle and specify which bin it is in
   allocate(bins(numbinsx,numbinsy,numbinsz))  ! 3D grid array specifying first particle in a bin grid
   allocate(particlebin(totalnumparticles,3))  ! 1D vector specifying bin number of each particle center
   allocate(linkedbinlist(totalnumparticles))


!!!!!BEGIN TIME INCREMENTATION AND PARTICLE/ WALL INTERACTIONS!!!!!
   counter = 1
   k = 0
   do while (t < tend)

!!!!!!!!BINNING STUFF!!!!!!!!!!! (rebin each time step) 
      do i = 1,totalnumparticles
         particlebin(i,1) = ceiling(CenterCoords(i,1)/binlength) ! assign bin/ grid coordinates to each particle
         particlebin(i,2) = ceiling(CenterCoords(i,2)/binwidth) 
         particlebin(i,3) = ceiling(CenterCoords(i,3)/binheight)
      end do

      !create bin grid and linked list according to Alejandro's method
      bins(:,:,:) = 0
      linkedbinlist(:) = 0
      do i = totalnumparticles,1,-1
         linkedbinlist(i) = bins(particlebin(i,1),particlebin(i,2),particlebin(i,3))
         bins(particlebin(i,1),particlebin(i,2),particlebin(i,3)) = i
      end do

      





!!!!!WALL MOVES DOWN ---EDIT THIS LATER---  
!!$      if (t < 0.8*tend) then
!!$         height = height + Vwall(3)*tstep
!!$      end if

!!!!! DETERMINE FORCES ACTING ON EACH PARTICLE

!!$!!!!! UPDATE POSITION AND VELOCITY USING FWD EULER SCHEME
!!$      call Force(numparticles,PositionOld,VOld,Radius,Vwall,t,tend,tstart,tstep,height,length,width,PsiTot1)
!!$
!!$      Vnew(:,1) = Vold(:,1) + 1.0/mass(:)*Psitot1(:,1)*tstep
!!$      Vnew(:,2) = Vold(:,2) + 1.0/mass(:)*Psitot1(:,2)*tstep
!!$      Vnew(:,3) = Vold(:,3) + 1.0/mass(:)*Psitot1(:,3)*tstep
!!$
!!$      PositionNew(:,1) = PositionOld(:,1) + Vold(:,1)*tstep
!!$      PositionNew(:,2) = PositionOld(:,2) + Vold(:,2)*tstep
!!$      PositionNew(:,3) = PositionOld(:,3) + Vold(:,3)*tstep


!!!!! UPDATE POSITION AND VELOCITY USING EXPLICIT RK-4 METHOD
      call Force(numparticles,PositionOld,VOld,Radius,Mass,t,tend,tstart,tstep,height,length,width,PsiTot1,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),rho)
      Y1(:,1:3) = Vold(:,1:3)
      Y1(:,4:6) = PositionOld(:,1:3)

      Y2(:,1) = Vold(:,1) + tstep/2.0*1.0/mass(:)*Psitot1(:,1)
      Y2(:,2) = Vold(:,2) + tstep/2.0*1.0/mass(:)*Psitot1(:,2)
      Y2(:,3) = Vold(:,3) + tstep/2.0*1.0/mass(:)*Psitot1(:,3)
      Y2(:,4) = PositionOld(:,1) + tstep/2.0*Y1(:,1)
      Y2(:,5) = PositionOld(:,2) + tstep/2.0*Y1(:,2)
      Y2(:,6) = PositionOld(:,3) + tstep/2.0*Y1(:,3)


!!!!!! CONSTRAIN PARTICLES FROM LEAVING BOX
!!$      do i = 1,numparticles
!!$         if (Y2(i,1) < 0 + Radius(i)) then
!!$            Y2(i,1) = Radius(i)
!!$         elseif (Y2(i,1) > length - Radius(i)) then
!!$            Y2(i,1) = length-Radius(i)
!!$         end if
!!$         if (Y2(i,2) < 0 + Radius(i)) then
!!$            Y2(i,2) = Radius(i)
!!$         elseif (Y2(i,2) > width - Radius(i)) then
!!$            Y2(i,2) = width-Radius(i)
!!$         end if
!!$         if (Y2(i,3) < 0 + Radius(i)) then
!!$            Y2(i,3) = Radius(i)
!!$         elseif (Y2(i,3) > height - Radius(i)) then
!!$            Y2(i,3) = height - Radius(i)
!!$         end if
!!$      end do
!!!!!!!!!!!

      call Force(numparticles,Y2(:,4:6),Y2(:,1:3),Radius,Mass,t+tstep/2.0,tend,tstart,tstep,height,length,width,PsiTot2,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),rho) !Psi(Vold=Y2(:,1:3),t=t+tstep/2)
      Y3(:,1) = Vold(:,1) + tstep/2.0*1.0/mass(:)*Psitot2(:,1)  
      Y3(:,2) = Vold(:,2) + tstep/2.0*1.0/mass(:)*Psitot2(:,2)  
      Y3(:,3) = Vold(:,3) + tstep/2.0*1.0/mass(:)*Psitot2(:,3) 
      Y3(:,4) = PositionOld(:,1) + tstep/2.0*Y2(:,1)
      Y3(:,5) = PositionOld(:,2) + tstep/2.0*Y2(:,2)
      Y3(:,6) = PositionOld(:,3) + tstep/2.0*Y2(:,3)


!!!!!! CONSTRAIN PARTICLES FROM LEAVING BOX
!!$      do i = 1,numparticles
!!$         if (Y3(i,1) < 0 + Radius(i)) then
!!$            Y3(i,1) = Radius(i)
!!$         elseif (Y2(i,1) > length - Radius(i)) then
!!$            Y3(i,1) = length-Radius(i)
!!$         end if
!!$         if (Y3(i,2) < 0 + Radius(i)) then
!!$            Y3(i,2) = Radius(i)
!!$         elseif (Y3(i,2) > width - Radius(i)) then
!!$            Y3(i,2) = width-Radius(i)
!!$         end if
!!$         if (Y3(i,3) < 0 + Radius(i)) then
!!$            Y3(i,3) = Radius(i)
!!$         elseif (Y3(i,3) > height - Radius(i)) then
!!$            Y3(i,3) = height - Radius(i)
!!$         end if
!!$      end do
!!!!!!!!!!!


      call Force(numparticles,Y3(:,4:6),Y3(:,1:3),Radius,Mass,t+tstep/2.0,tend,tstart,tstep,height,length,width,PsiTot3,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),rho) !Psi(Vold=Y3(:,1:3),t=t+tstep/2)
      Y4(:,1) = Vold(:,1) + tstep*1.0/mass(:)*Psitot3(:,1)  
      Y4(:,2) = Vold(:,2) + tstep*1.0/mass(:)*Psitot3(:,2)  
      Y4(:,3) = Vold(:,3) + tstep*1.0/mass(:)*Psitot3(:,3)  
      Y4(:,4) = PositionOld(:,1) + tstep*Y3(:,1)
      Y4(:,5) = PositionOld(:,2) + tstep*Y3(:,2)
      Y4(:,6) = PositionOld(:,3) + tstep*Y3(:,3)


!!!!!! CONSTRAIN PARTICLES FROM LEAVING BOX
!!$      do i = 1,numparticles
!!$         if (Y4(i,1) < 0 + Radius(i)) then
!!$            Y4(i,1) = Radius(i)
!!$         elseif (Y4(i,1) > length - Radius(i)) then
!!$            Y4(i,1) = length-Radius(i)
!!$         end if
!!$         if (Y4(i,2) < 0 + Radius(i)) then
!!$            Y4(i,2) = Radius(i)
!!$         elseif (Y4(i,2) > width - Radius(i)) then
!!$            Y4(i,2) = width-Radius(i)
!!$         end if
!!$         if (Y4(i,3) < 0 + Radius(i)) then
!!$            Y4(i,3) = Radius(i)
!!$         elseif (Y4(i,3) > height - Radius(i)) then
!!$            Y4(i,3) = height - Radius(i)
!!$         end if
!!$      end do
!!!!!!!!!!!

      call Force(numparticles,Y4(:,4:6),Y4(:,1:3),Radius,Mass,t+tstep,tend,tstart,tstep,height,length,width,PsiTot4,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),rho) !Psi(Vold=Y4(:,1:3),t=t+tstep)
      Vnew(:,1) = Vold(:,1) + tstep/6.0*(1.0/mass(:)*Psitot1(:,1) + 2.0/mass(:)*Psitot2(:,1) + 2.0/mass(:)*Psitot3(:,1) + 1.0/mass(:)*Psitot4(:,1)) 
      Vnew(:,2) = Vold(:,2) + tstep/6.0*(1.0/mass(:)*Psitot1(:,2) + 2.0/mass(:)*Psitot2(:,2) + 2.0/mass(:)*Psitot3(:,2) + 1.0/mass(:)*Psitot4(:,2)) 
      Vnew(:,3) = Vold(:,3) + tstep/6.0*(1.0/mass(:)*Psitot1(:,3) + 2.0/mass(:)*Psitot2(:,3) + 2.0/mass(:)*Psitot3(:,3) + 1.0/mass(:)*Psitot4(:,3)) 

  !    Y1(:,4) = PositionOld(:,1)
 !     Y2(:,4) = PositionOld(:,1) + tstep/2.0*Y1(:,1)
   !   Y3(:,4) = PositionOld(:,1) + tstep/2.0*Y2(:,1)
!      Y4(:,4) = PositionOld(:,1) + tstep*Y3(:,1)

 !     Y1(:,5) = PositionOld(:,2)
 !     Y2(:,5) = PositionOld(:,2) + tstep/2.0*Y1(:,2)
  !    Y3(:,5) = PositionOld(:,2) + tstep/2.0*Y2(:,2)
  !    Y4(:,5) = PositionOld(:,2) + tstep*Y3(:,2)

  !    Y1(:,6) = PositionOld(:,3)
  !    Y2(:,6) = PositionOld(:,3) + tstep/2.0*Y1(:,3)
  !    Y3(:,6) = PositionOld(:,3) + tstep/2.0*Y2(:,3)
   !   Y4(:,6) = PositionOld(:,3) + tstep*Y3(:,3)

      PositionNew(:,1) = PositionOld(:,1) + tstep/6.0*(Y1(:,1) + 2.0*Y2(:,1) + 2.0*Y3(:,1) + Y4(:,1))
      PositionNew(:,2) = PositionOld(:,2) + tstep/6.0*(Y1(:,2) + 2.0*Y2(:,2) + 2.0*Y3(:,2) + Y4(:,2))
      PositionNew(:,3) = PositionOld(:,3) + tstep/6.0*(Y1(:,3) + 2.0*Y2(:,3) + 2.0*Y3(:,3) + Y4(:,3))



!!$!!!!CONSTRAIN PARTICLES FROM LEAVING BOX USING ELASTIC COLLISION AT WALLS
!!$   do i=1,numparticles
!!$
!!$      if (PositionNew(i,1) < 0 + Radius(i)) then
!!$         PositionNew(i,1) = Radius(i)
!!$         !Vold(i,1) = - 0.7*Vold(i,1)
!!$         !Psitot(i,1) = - Psitot(i,1)
!!$      elseif (PositionNew(i,1) > length - Radius(i)) then
!!$         PositionNew(i,1) = length-Radius(i)
!!$         !Vold(i,1) = - 0.7*Vold(i,1)
!!$         !Psitot(i,1) = - Psitot(i,1)
!!$      end if
!!$
!!$      if (PositionNew(i,2) < 0 + Radius(i)) then
!!$         PositionNew(i,2) = Radius(i)
!!$         !Vold(i,2) = - 0.7*Vold(i,2)
!!$         !Psitot(i,2) = - Psitot(i,2)
!!$      elseif (PositionNew(i,2) > width - Radius(i)) then
!!$         PositionNew(i,2) = width-Radius(i)
!!$         !Vold(i,2) = -0.7*Vold(i,2)
!!$         !Psitot(i,2) = - Psitot(i,2)
!!$      end if
!!$
!!$      if (PositionNew(i,3) < 0 + Radius(i)) then
!!$         PositionOld(i,3) = Radius(i)
!!$         !Vold(i,3) = - 0.7*Vold(i,3)
!!$         !Psitot(i,3) = - Psitot(i,3)
!!$      elseif (PositionNew(i,3) > height - Radius(i)) then
!!$         PositionNew(i,3) = height - Radius(i)
!!$         !Vold(i,3) = - 0.7*Vold(i,3)
!!$         !Psitot(i,3) = - Psitot(i,3)
!!$      end if
!!$   
!!$   end do



!!!OUTPUT FILES TO TECPLOT every 10000 time steps (160 total)
      if ((counter==10000*k) .OR. (counter==1)) then
         write(unit=100,fmt='(A)') 'Zone'
         do i = 1,numParticles 
            write(unit=100,fmt=70) 1E3*CenterCoords(i,1), 1E3*CenterCoords(i,2), 1E3*CenterCoords(i,3), 1E3*Radius(i), Inside(i), Temp(i)
         end do
         print *, k
         k = k + 1
      end if



      CenterCoords(:,:) = PositionNew(:,:)
      PositionOld(:,:) = PositionNew(:,:)
      Vold(:,:) = Vnew(:,:)

      counter = counter + 1
      t = t + tstep

   end do

end do

write(unit=200,fmt='(I10)') numparticles

do i = 1,numParticles 
   write(unit=200,fmt=70) 1e3*CenterCoords(i,1), 1e3*CenterCoords(i,2), 1e3*CenterCoords(i,3), 1e3*Radius(i)
end do

close(unit=100)
close(unit=200)

end program PowderDepositionShankarMP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! FORCE CALCULATION SUBROUTINE !!!!!!
SUBROUTINE Force(numparticles,PositionOld,VOld,Radius,Mass,t,tend,tstart,tstep,height,length,width,PsiTot,numbinsx,numbinsy,numbinsz,bins,linkedbinlist,ParticleBin,rho)
  implicit none

  integer, intent(in) :: numparticles, numbinsx, numbinsy, numbinsz
  integer, dimension(numparticles), intent(in) :: linkedbinlist
  integer, dimension(numparticles,3), intent(in) :: ParticleBin
  integer, dimension(numbinsx,numbinsy,numbinsz), intent(in) :: bins
  double precision, intent(in) :: t, tend, tstart, tstep, height, length, width, rho
  !double precision, dimension(3), intent(in) :: Vwall
  double precision, dimension(numparticles), intent(in) :: Radius, Mass
  double precision, dimension(numparticles,3), intent(in) ::  VOld, PositionOld
  double precision, dimension(numparticles,3), intent(out) :: PsiTot
  double precision, dimension(numparticles,3) :: CenterCoords, PsiCont, PsiFric, PsiAdh, PsiEnv, PsiGrav
  double precision, dimension(3) :: normal, tangent, psicontwall, psifricwall, psicontpart, psifricpart, psiadhwall, psiadhpart, psienvwall, psienvpart, VTan, Vtan2
  double precision, dimension(numparticles) :: deltawallx1, deltawallx2, deltawally1, deltawally2, deltawallz1, deltawallz2, deltawalloldx1, deltawalloldx2, deltawalloldy1, deltawalloldy2, deltawalloldz1, deltawalloldz2
  double precision, dimension(numparticles,numparticles) :: deltaparticleold, deltaparticle
  double precision :: Bpart, Bwall, MuD, MuS, zdist, xdist, ydist, distance, epsilonwall, epsilonpart, speed, DampCoef, g, NuPart, NuWall, EPart, EWall, EstarWall, EstarPart, Rstar, deltawall, deltapart, rhoN2, viscosityN2, Reynolds, deltawalldot, deltapartdot, gamma1, gamma2, mstar, zeta
  integer :: i, j, k, m, n, h1, h2, h3
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462



  !! Parameters
  Bwall = 5.0E-6 !Calculated analytically (see 10/29/13) Deposition
  Bpart = 0.333*Bwall !0.2*Bwall
  DampCoef = 5.0E-8
  MuD = 0.1
  MuS = 0.2
  g = 9.81 ! downward gravitational acceleration [m/s^2]
  NuPart = 0.260D0  ! 316L SS at 366 K
  NuWall = 0.260D0 ! 316L SS at 366 K
  EPart = 193.0e9  ! [Pa], 316L SS at 366 K
  Ewall = 193.0e9 ! [Pa] 316L SS
  EstarWall = (EPart*Ewall)/(Ewall*(1-NuPart**2.0)+Epart*(1-NuWall**2.0))  ! Effective Modulus for wall
  EstarPart = (EPart*EPart)/(EPart*(1-NuPart**2.0)+EPart*(1-NuPart**2.0))  ! Effective Modulus for particles all made of same material, HX
  rhoN2 = 1.165 ! [kg/m^3] at NTP (20 C, 1 atm) assuming a Nitrogen surrounding atmosphere
  viscosityN2 = 1.747e-5 ! [kg/m-s or N-s/m^2] at NTP (20 C, 1 atm)
  gamma1 = 1.0e-5 ! damping coefficient for Hertzian contact collisions
  gamma2 = 1.0e1  ! damping coefficient for Hertzian contact collisions with delta^(0.5) term
  zeta = 0.1 ! damping coefficienct for Hertzian contact damping, Eduardo/ Wriggers method (z = 1 -- critically damped, z > 1 -- overdamped, z < 1 -- underdamped)

  normal(:) = 0.0D0
  tangent(:) = 0.0D0
  psicontwall(:) = 0.0D0
  psicontpart(:) = 0.0D0
  psifricwall(:) = 0.0D0
  psifricpart(:) = 0.0D0
  psiadhwall(:) = 0.0D0
  psiadhpart(:) = 0.0D0
  psienvwall(:) = 0.0D0
  psienvpart(:) = 0.0D0
  psicont(:,:) = 0.0D0
  psifric(:,:) = 0.0D0
  psiadh(:,:) = 0.0D0
  psienv(:,:) = 0.0D0
  psigrav(:,:) = 0.0D0
  psitot(:,:) = 0.0D0
  CenterCoords(:,:) = PositionOld(:,:)
 ! deltawallx1(:) = 0.0D0
 ! deltawallx2(:) = 0.0D0
 ! deltawally1(:) = 0.0D0
 ! deltawally2(:) = 0.0D0
 ! deltawallz1(:) = 0.0D0
 ! deltawallz2(:) = 0.0D0
 ! deltawalloldx1(:) = 0.0D0
 ! deltawalloldx2(:) = 0.0D0
 ! deltawalloldy1(:) = 0.0D0
 ! deltawalloldy2(:) = 0.0D0
 ! deltawalloldz1(:) = 0.0D0
 ! deltawalloldz2(:) = 0.0D0
 ! deltaparticle(:,:) = 0.0D0
 ! deltaparticleold(:,:) = 0.0D0

!!$!!!!!WALL MOVES DOWN ---EDIT THIS LATER---  
!!$if (t < 0.8*tend) then
!!$   height = oldheight+Vwall(3)*tstep
!!$else
!!$   height = oldheight
!!$end if

  ! Calculate Gravitational Force and environmental damping (drag) force
  do i = 1,numparticles
     PsiGrav(i,3) = -Mass(i)*g
     Reynolds = rhoN2*norm2(Vold(i,:))*Radius(i)*2.0/viscosityN2
     PsiEnv(i,:) = -6.0*PI*viscosityN2*Radius(i)*Vold(i,:) ! Stoke's Law for drag force of Low Reynold's Number (small) spheres
 !    PsiEnv(i,:) = -DampCoef*(Vold(i,:)-0.0) 

  end do


!!!!!CHECK FOR COLLISIONS, 1st particle-wall, 2nd particle-particle
!!!! Particle-top/ bottom (z) wall collisions
!$OMP PARALLEL DO
  do i = 1,numparticles
     zdist = abs(CenterCoords(i,3)-height)
     if (zdist < Radius(i)) then
        !      epsilonwall = abs(zdist-Radius(i))/Radius(i)
        normal(1) = 0.0D0
        normal(2) = 0.0D0
        normal(3) = 1.0D0
        !      psicontwall(3) = Bwall*epsilonwall*normal(3)
        deltawall = abs(zdist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - gamma1*deltawalldot*normal(3)
!!$        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - gamma2*deltawall**(0.5)*deltawalldot*normal(3)
        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(3)
        PsiCont(i,3) = PsiCont(i,3) + psicontwall(3)  ! CONTACT FORCE

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(1) = (0.0 - Vtan(1)) / norm2(-Vtan)
!!$        tangent(2) = (0.0 - Vtan(2)) / norm2(-Vtan)
!!$        psifricwall(1) = MuD*abs(psicontwall(3))*tangent(1)
!!$        psifricwall(2) = Mud*abs(psicontwall(3))*tangent(2)
!!$        PsiFric(i,1) = PsiFric(i,1) + psifricwall(1)
!!$        PsiFric(i,2) = PsiFric(i,2) + psifricwall(2)  !  FRICTION FORCE
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall
   
    !    PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) ! DAMPING FORCE

     end if
  end do
!$OMP END PARALLEL DO


!$OMP PARALLEL DO
  do i = 1,numparticles
     zdist = abs(CenterCoords(i,3)-0.0)
     if (zdist < Radius(i)) then
        !     epsilonwall = abs(zdist-Radius(i))/Radius(i)
        normal(1) = 0.0D0
        normal(2) = 0.0D0
        normal(3) = -1.0D0
        !     psicontwall(3) = Bwall*epsilonwall*normal(3)
        deltawall = abs(zdist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - gamma1*deltawalldot*normal(3)
!!$        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - gamma2*deltawall**(0.5)*deltawalldot*normal(3)
        psicontwall(3) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(3) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(3)
        PsiCont(i,3) = PsiCont(i,3) + psicontwall(3)  ! CONTACT FORCE

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(1) = (0.0 - Vtan(1)) / norm2(-Vtan)
!!$        tangent(2) = (0.0 - Vtan(2)) / norm2(-Vtan)
!!$        psifricwall(1) = MuD*abs(psicontwall(3))*tangent(1)
!!$        psifricwall(2) = Mud*abs(psicontwall(3))*tangent(2)
!!$        PsiFric(i,1) = PsiFric(i,1) + psifricwall(1)
!!$        PsiFric(i,2) = PsiFric(i,2) + psifricwall(2)  !  FRICTION FORCE
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall

  !      PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) 

     end if
  end do
!$OMP END PARALLEL DO


!!!! Particle-x wall collisions
!$OMP PARALLEL DO
  do i = 1,numparticles
     xdist = abs(CenterCoords(i,1)-length)
     if (xdist < Radius(i)) then
        !     epsilonwall = abs(xdist-Radius(i))/Radius(i)
        normal(1) = 1.0D0
        normal(2) = 0.0D0
        normal(3) = 0.0D0
        !     psicontwall(1) = Bwall*epsilonwall*normal(1)
        deltawall = abs(xdist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - gamma1*deltawalldot*normal(1)
!!$        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - gamma2*deltawall**(0.5)*deltawalldot*normal(1)
        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(1)
        PsiCont(i,1) = PsiCont(i,1) + psicontwall(1)

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(2) = (0.0 - Vtan(2)) / norm2(-Vtan)
!!$        tangent(3) = (0.0 - Vtan(3)) / norm2(-Vtan)
!!$        psifricwall(2) = MuD*abs(psicontwall(1))*tangent(2)
!!$        psifricwall(3) = Mud*abs(psicontwall(1))*tangent(3) 
!!$        PsiFric(i,2) = PsiFric(i,2) + psifricwall(2)
!!$        PsiFric(i,3) = PsiFric(i,3) + psifricwall(3)
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall

    !    PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) 

     end if
  end do
!$OMP END PARALLEL DO

!$OMP  PARALLEL DO
  do i = 1,numparticles
     xdist = abs(CenterCoords(i,1)-0.0)
     if (xdist < Radius(i)) then
        !      epsilonwall = abs(xdist-Radius(i))/Radius(i)
        normal(1) = -1.0D0
        normal(2) = 0.0D0
        normal(3) = 0.0D0
        !      psicontwall(1) = Bwall*epsilonwall*normal(1)
        deltawall = abs(xdist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - gamma1*deltawalldot*normal(1)
!!$        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - gamma2*deltawall**(0.5)*deltawalldot*normal(1)
        psicontwall(1) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(1) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(1)
        PsiCont(i,1) = PsiCont(i,1) + psicontwall(1)

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(2) = (0.0 - Vtan(2)) / norm2(-Vtan)
!!$        tangent(3) = (0.0 - Vtan(3)) / norm2(-Vtan)
!!$        psifricwall(2) = MuD*abs(psicontwall(1))*tangent(2)
!!$        psifricwall(3) = Mud*abs(psicontwall(1))*tangent(3) 
!!$        PsiFric(i,2) = PsiFric(i,2) + psifricwall(2)
!!$        PsiFric(i,3) = PsiFric(i,3) + psifricwall(3)
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall

    !    PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) 

     end if
  end do
!$OMP END PARALLEL DO


!!!! Particle-y wall collisions
!$OMP PARALLEL DO
  do i = 1,numparticles
     ydist = abs(CenterCoords(i,2)-width)
     if (ydist < Radius(i)) then
        !     epsilonwall = abs(ydist-Radius(i))/Radius(i)
        normal(1) = 0.0D0
        normal(2) = 1.0D0
        normal(3) = 0.0D0
        !      psicontwall(2) = Bwall*epsilonwall*normal(2)
        deltawall = abs(ydist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - gamma1*deltawalldot*normal(2)
!!$        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - gamma2*deltawall**(0.5)*deltawalldot*normal(2)
        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(2)
        PsiCont(i,2) = PsiCont(i,2) + psicontwall(2)

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(1) = (0.0 - Vtan(1)) / norm2(-Vtan)
!!$        tangent(3) = (0.0 - Vtan(3)) / norm2(-Vtan)
!!$        psifricwall(1) = MuD*abs(psicontwall(2))*tangent(1)
!!$        psifricwall(3) = Mud*abs(psicontwall(2))*tangent(3) 
!!$        PsiFric(i,1) = PsiFric(i,1) + psifricwall(1)
!!$        PsiFric(i,3) = PsiFric(i,3) + psifricwall(3)
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall

  !      PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) 

     end if
  end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO
  do i = 1,numparticles
     ydist = abs(CenterCoords(i,2)-0.0)
     if (ydist < Radius(i)) then
        !     epsilonwall = abs(ydist-Radius(i))/Radius(i)
        normal(1) = 0.0D0
        normal(2) = -1.0D0
        normal(3) = 0.0D0
        !     psicontwall(2) = Bwall*epsilonwall*normal(2)
        deltawall = abs(ydist-Radius(i))
        deltawalldot = dot_product(Vold(i,:),normal)
        RStar = Radius(i)
        Mstar = Mass(i)
!!$        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - gamma1*deltawalldot*normal(2)
!!$        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - gamma2*deltawall**(0.5)*deltawalldot*normal(2)
        psicontwall(2) = -4.0/3.0*Rstar**(0.5)*EstarWall*deltawall**(1.5)*normal(2) - 2.0*zeta*(2.0*EstarWall*Mstar)**(0.5)*(Rstar*deltawall)**(0.25)*deltawalldot*normal(2)
        PsiCont(i,2) = PsiCont(i,2) + psicontwall(2)

        Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$        tangent(1) = (0.0 - Vtan(1)) / norm2(-Vtan)
!!$        tangent(3) = (0.0 - Vtan(3)) / norm2(-Vtan)
!!$        psifricwall(1) = MuD*abs(psicontwall(2))*tangent(1)
!!$        psifricwall(3) = Mud*abs(psicontwall(2))*tangent(3) 
!!$        PsiFric(i,1) = PsiFric(i,1) + psifricwall(1)
!!$        PsiFric(i,3) = PsiFric(i,3) + psifricwall(3)
        if (norm2(vtan) > 0.0) then
           tangent = (-Vtan)/norm2(Vtan)
           psifricwall = MuD*norm2(psicontwall)*tangent
        else
           psifricwall(:) = 0.0D0
        end if
        PsiFric(i,:) = PsiFric(i,:) + Psifricwall

    !    PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-0.0) 
     end if
  end do
!$OMP END PARALLEL DO



!!!!PARTICLE TO PARTICLE COLLISIONS - NO BINNING
!!$i = 1
!!$j = 1
!!$do i = 1,numparticles
!!$   do j = 1,numparticles
!!$      distance = norm2(CenterCoords(i,:)-CenterCoords(j,:))
!!$      speed = norm2(Vold(i,:)-Vold(j,:))
!!$
!!$      if ((distance < (Radius(i)+Radius(j))) .AND. i/=j) then
!!$         epsilonpart = abs(distance-(Radius(i)+Radius(j))) / (Radius(i)+Radius(j))
!!$         normal(1) = (CenterCoords(i,1)-CenterCoords(j,1)) / distance
!!$         normal(2) = (CenterCoords(i,2)-CenterCoords(j,2)) / distance
!!$         normal(3) = (CenterCoords(i,3)-CenterCoords(j,3)) / distance
!!$         psicontpart(1) = Bpart*epsilonpart*normal(1)
!!$         psicontpart(2) = Bpart*epsilonpart*normal(2)
!!$         psicontpart(3) = Bpart*epsilonpart*normal(3)
!!$         PsiCont(i,1) = Psicont(i,1) + Psicontpart(1)
!!$         PsiCont(i,2) = Psicont(i,2) + Psicontpart(2)
!!$         PsiCont(i,3) = Psicont(i,3) + Psicontpart(3)
!!$
!!$         Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
!!$         Vtan2 = Vold(j,:) - (dot_product(Vold(j,:),normal))*normal
!!$
!!$         tangent(1) = (Vtan2(1))-Vtan(1)) / norm2(Vtan-Vtan2)
!!$         tangent(2) = (Vtan2(2))-Vtan(2)) / norm2(Vtan-Vtan2)
!!$         tangent(3) = (Vtan2(3))-Vtan(3)) / norm2(Vtan-Vtan2)
!!$         psifricpart(1) = MuD*norm2(psicontpart)*tangent(1)
!!$         psifricpart(2) = MuD*norm2(psicontpart)*tangent(2)
!!$         psifricpart(3) = MuD*norm2(psicontpart)*tangent(3)
!!$         Psifric(i,1) = Psifric(i,1) + psifricpart(1)
!!$         Psifric(i,2) = Psifric(i,2) + psifricpart(2)
!!$         Psifric(i,3) = Psifric(i,3) + psifricpart(3)
!!$
!!$         PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-Vold(j,:)) 
!!$      end if
!!$   end do
!!$end do

!!!!PARTICLE TO PARTICLE COLLISIONS - WITH BINNING AND OPENMP PARALLELIZATION

!$OMP PARALLEL DO
  do i = 1,numparticles

     ! Search through neighboring bins
     h1 = particlebin(i,1)
     h2 = particlebin(i,2)
     h3 = particlebin(i,3)

     do k = -1,1,1
        do m = -1,1,1
           do n = -1,1,1
              if ((h1+k > 0) .AND. (h1+k <= numbinsx) .AND. (h2+m > 0) .AND. (h2+m <= numbinsy) .AND. (h3+n > 0) .AND. (h3+n <= numbinsz)) then

                 j = bins(h1+k,h2+m,h3+n)

                 do while ((j/=0) .AND. (j <= numparticles))

                    distance = norm2(CenterCoords(i,:)-CenterCoords(j,:))
                    speed = norm2(Vold(i,:)-Vold(j,:))

                    if ((distance < (Radius(i)+Radius(j))) .AND. i/=j) then
                       !                     epsilonpart = abs(distance-(Radius(i)+Radius(j))) / (Radius(i)+Radius(j))
                       normal(1) = (CenterCoords(j,1)-CenterCoords(i,1)) / distance
                       normal(2) = (CenterCoords(j,2)-CenterCoords(i,2)) / distance
                       normal(3) = (CenterCoords(j,3)-CenterCoords(i,3)) / distance

!!!!!!!!!!!!!!!!! HERTZ MODIFICATIONS
                       deltaPart = abs(distance-(Radius(i)+Radius(j)))
                       deltapartdot = dot_product((Vold(i,:)-Vold(j,:)),normal)

                       RStar = Radius(i)*Radius(j)/(Radius(i)+Radius(j))
                       Mstar = Mass(i)*Mass(j)/(Mass(i)+Mass(j))
!!$                       psicontpart(1) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(1) - gamma1*deltapartdot*normal(1)
!!$                       psicontpart(2) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(2) - gamma1*deltapartdot*normal(2)
!!$                       psicontpart(3) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(3) - gamma1*deltapartdot*normal(3)
!!$                       psicontpart(1) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(1) - gamma2*deltapart**(0.5)*deltapartdot*normal(1)
!!$                       psicontpart(2) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(2) - gamma2*deltapart**(0.5)*deltapartdot*normal(2)
!!$                       psicontpart(3) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(3) - gamma2*deltapart**(0.5)*deltapartdot*normal(3)
                       psicontpart(1) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(1) - 2.0*zeta*(2.0*EstarPart*Mstar)**(0.5)*(Rstar*deltaPart)**(0.25)*deltapartdot*normal(1)
                       psicontpart(2) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(2) - 2.0*zeta*(2.0*EstarPart*Mstar)**(0.5)*(Rstar*deltaPart)**(0.25)*deltapartdot*normal(2)
                       psicontpart(3) = -4.0/3.0*Rstar**(0.5)*EstarPart*deltaPart**(1.5)*normal(3) - 2.0*zeta*(2.0*EstarPart*Mstar)**(0.5)*(Rstar*deltaPart)**(0.25)*deltapartdot*normal(3)     
                       !                    psicontpart(1) = Bpart*epsilonpart*normal(1)
                       !                    psicontpart(2) = Bpart*epsilonpart*normal(2)
                       !                    psicontpart(3) = Bpart*epsilonpart*normal(3)
                       PsiCont(i,1) = Psicont(i,1) + Psicontpart(1)
                       PsiCont(i,2) = Psicont(i,2) + Psicontpart(2)
                       PsiCont(i,3) = Psicont(i,3) + Psicontpart(3)

                       Vtan = Vold(i,:) - (dot_product(Vold(i,:),normal))*normal
                       Vtan2 = Vold(j,:) - (dot_product(Vold(j,:),(-normal)))*(-normal)
                       if (norm2(Vtan-Vtan2) > 0.0) then
                          tangent = (Vtan2-Vtan) / norm2(Vtan2-Vtan)
                          psifricpart = MuD*norm2(psicontpart)*tangent
                       else
                          psifricpart(:) = 0.0D0
                       end if
                       PsiFric(i,:) = Psifric(i,:) + Psifricpart

!!$                       tangent(1) = (Vtan2(1)-Vtan(1)) / norm2(Vtan-Vtan2)
!!$                       tangent(2) = (Vtan2(2)-Vtan(2)) / norm2(Vtan-Vtan2)
!!$                       tangent(3) = (Vtan2(3)-Vtan(3)) / norm2(Vtan-Vtan2)
!!$                       psifricpart(1) = MuD*norm2(psicontpart)*tangent(1)
!!$                       psifricpart(2) = MuD*norm2(psicontpart)*tangent(2)
!!$                       psifricpart(3) = MuD*norm2(psicontpart)*tangent(3)
!!$                       Psifric(i,1) = Psifric(i,1) + psifricpart(1)
!!$                       Psifric(i,2) = Psifric(i,2) + psifricpart(2)
!!$                       Psifric(i,3) = Psifric(i,3) + psifricpart(3)



       !                PsiEnv(i,:) = PsiEnv(i,:) + -DampCoef*(Vold(i,:)-Vold(j,:)) 
                    end if

                    j = linkedbinlist(j)

                 end do

              end if

           end do
        end do
     end do

  end do
!$OMP END PARALLEL DO

!$OMP BARRIER  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!SUM UP TOTAL FORCES

  Psitot(:,1) = PsiCont(:,1) + Psifric(:,1) + PsiEnv(:,1) + PsiGrav(:,1)
  Psitot(:,2) = PsiCont(:,2) + Psifric(:,2) + PsiEnv(:,2) + PsiGrav(:,2)
  Psitot(:,3) = PsiCont(:,3) + Psifric(:,3) + PsiEnv(:,3) + PsiGrav(:,3)

  ! Isolate Problem Forces - no friction to begin with
!!$  Psitot(:,1) = PsiCont(:,1)  + PsiEnv(:,1) + PsiGrav(:,1)
!!$  Psitot(:,2) = PsiCont(:,2) + PsiEnv(:,2) + PsiGrav(:,2)
!!$  Psitot(:,3) = PsiCont(:,3) + PsiEnv(:,3) + PsiGrav(:,3)

!!$  Psitot(:,1) = PsiCont(:,1)  + PsiGrav(:,1)
!!$  Psitot(:,2) = PsiCont(:,2) + PsiGrav(:,2)
!!$  Psitot(:,3) = PsiCont(:,3) + PsiGrav(:,3)


END SUBROUTINE Force




!!RANDOM SEED SUBROUTINE
SUBROUTINE init_random_seed()
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          


            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)
END SUBROUTINE




!!!!CONSTRAIN PARTICLES FROM LEAVING BOX USING ELASTIC COLLISION AT WALLS
!!$   do i=1,numparticles
!!$
!!$      if (PositionOld(i,1) < 0 + Radius(i)) then
!!$         PositionOld(i,1) = Radius(i)
!!$         Vold(i,1) = - 0.7*Vold(i,1)
!!$         Psitot(i,1) = - Psitot(i,1)
!!$      elseif (PositionOld(i,1) > length - Radius(i)) then
!!$         PositionOld(i,1) = length-Radius(i)
!!$         Vold(i,1) = - 0.7*Vold(i,1)
!!$         Psitot(i,1) = - Psitot(i,1)
!!$      end if
!!$
!!$      if (PositionOld(i,2) < 0 + Radius(i)) then
!!$         PositionOld(i,2) = Radius(i)
!!$         Vold(i,2) = - 0.7*Vold(i,2)
!!$         Psitot(i,2) = - Psitot(i,2)
!!$      elseif (PositionOld(i,2) > width - Radius(i)) then
!!$         PositionOld(i,2) = width-Radius(i)
!!$         Vold(i,2) = -0.7*Vold(i,2)
!!$         Psitot(i,2) = - Psitot(i,2)
!!$      end if
!!$
!!$      if (PositionOld(i,3) < 0 + Radius(i)) then
!!$         PositionOld(i,3) = Radius(i)
!!$         Vold(i,3) = - 0.7*Vold(i,3)
!!$         Psitot(i,3) = - Psitot(i,3)
!!$      elseif (PositionOld(i,3) > height - Radius(i)) then
!!$         PositionOld(i,3) = height - Radius(i)
!!$         Vold(i,3) = - 0.7*Vold(i,3)
!!$         Psitot(i,3) = - Psitot(i,3)
!!$      end if
!!$   
!!$   end do
!!$         

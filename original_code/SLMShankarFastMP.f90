program SLMShankarFastMP
  use omp_lib
  implicit none

  !! coupled DEM code + FD Code + Particles soften as they cool - includes CombinedThermoMech

!!!MultiSLM + CombinedThermoMech + softening of particles


  !! DEM VARIABLES
  double precision, dimension(:,:), allocatable :: Particles, origparticles, PositionOld, VOld, PsiTot1, PsiTot2, PsiTot3, PsiTot4, YF1, YF2, YF3, YF4, PositionNew, VNew
  integer :: numparticles, i, numcol, counter, j, k, n_z, num_y, layers, totalnumparticles, n, orignumparticles, numbinsx, numbinsy, numbinsz, totalnumbins, laserpass, count
  integer, dimension(:), allocatable :: inside, state, linkedbinlist
  integer, dimension(:,:), allocatable :: particlebin
  double precision :: tstart, tend, t, tstep, avgthermalcond, Acon, Alaser, A0, Amax, xdist, ydist, zdist, distance, thermalcondsol, thermalcondliq, hatchspace
  double precision :: Csolid, Cliquid, LatentHeatMelt, LatentHeatVap, LatentHeatSolid, HeatCapacity, Tmelt, Tvapor, length, width, height, dist, alpha, rho, ThermalConductivity, absorptivity, I0, deltaTemp, CTE, porosity, avgRadius, avgHeight
  double precision :: LaserPower, LaserDiameter, LaserSpeed, ThetaStart, difftol, rho1, rho2, rho3, binlength, binwidth, binheight, LaserHeight, ThermalCondSubstrate, ThetaSubstrate, EYoungSubstrate, EYoungStart, EYoungLiquid, thermalcondgas, Cgas
  double precision :: Tmelt1, Tmelt2, Tmelt3, Tsint1, Tsint2, Tsint3, ThermalCond1, ThermalCond2, ThermalCond3, HeatCap1, HeatCap2, HeatCap3, temp, maxheight
  double precision, dimension(:), allocatable :: thermalcond, Q1, Q2, Q3, Q4, H, Vol, Radius, ThetaNew, ThetaOld, mass, density, volume, HeatCap, EYoungPart
  double precision, dimension(:), allocatable :: Y1, Y2, Y3, Y4, RandNumbers, OrigThermalCond, OrigHeatCap, AverageTemp, Melted, PercentMelted, Particle1, PercentVaporized, Vaporized
  integer, dimension(:,:,:), allocatable :: bins
  double precision, dimension(2) :: LaserPosition
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462

  !! FD VARIABLES
  double precision :: lengthFD, widthFD, heightFD, hx, hy, hz, DivGradX, DivGradY, DivGradZ, ExtinctionCoefFD, Theta0, U0, densitysolid, ycoord, xcoord, zcoord, densityliquid
  double precision, dimension(:,:,:), allocatable :: DensityFD, HeatCapFD, ThermalCondFD, U, ThetaFD, ThetaOldFD, ThetaDotFD, DivGradX1, DivGradY1, DivGradZ1 
  integer, dimension(:,:,:), allocatable :: stateFD
  double precision, dimension(:,:), allocatable :: CoordsFD, GhostNodesTop, FluxBC
  double precision, dimension(:), allocatable :: AvgT, AvgT1, AvgT2, AvgT3, AvgT4, AvgT1Laser, AvgT2Laser, AvgT3Laser, AvgT4Laser, AvgT0, AvgT5, time
  integer :: iterationDEM, iterationFD, numtimesteps, nodesx, nodesy, nodesz, nodestotal


  !!OPEN TECPLOT OUTPUT FILE FOR DEM
  Open(unit=100, file='SLMShankarResults2800MP.dat')
  write(unit=100,fmt='(A)') 'Title="SLM Particles"'
  write(unit=100,fmt='(A)') 'Variables="X","Y","Z","R","T","In"'

!!$  Open(unit=150, file='SLMShankarResultsSideView2800MP.dat')
!!$  write(unit=150,fmt='(A)') 'Title="SLM Particles"'
!!$  write(unit=150,fmt='(A)') 'Variables="X","Y","Z","R","T","In"'

  Open(unit=200, file='SLMShankarResultsFinal2800MP.dat')
!!$  open(unit=300, file='MultiSLMParticleAvgTemp.dat')
!!$  open(unit=400, file='MultiSLMParticleMeltVap.dat')

  !!OPEN TECPLOT OUTPUT FILE FOR FDM
!!$  Open(unit=10, file='MultiSLMCoolFDMatlabResults.dat')
!!$  write(unit=10,fmt='(A)') 'Title="MultiSLM Results"'
!!$  write(unit=10,fmt='(A)') 'Variables="time","AvgT","AvgT0","AvgT1"'

!!$  Open(unit=20, file='MultiSLMCoolFDResults.dat')
!!$  write(unit=20,fmt='(A)') 'Title="FD Temperature Results"'
!!$80 format (F13.6,2x,F13.6,2x,F13.6,2x,F13.6)
!!$
!!$  Open(unit=30, file='MultiSLMCoolFDResultsFinal.dat')
!!$  write(unit=30,fmt='(A)') 'Title="Final FD Temperature Results"'
!!$

  ! Read in particle coordinates 
  open(unit=110,file='PowderDepositionShankarMPFinal2800.dat')
  read(unit=110,fmt='(I10)') totalnumparticles


  !! INITIALIZE DEM PARAMETERS
  layers = 1
  numparticles = totalnumparticles/layers ! number of particles per layer assuming the same in each layer
  orignumparticles = numparticles

  length = 0.00040D0 ! dimensions of box in m
  width = 0.00080D0
  height = 0.0020D0

  !!Allocate and define variables
  allocate(particles(totalnumparticles,4))
  allocate(inside(totalnumparticles))
  allocate(state(totalnumparticles))
  allocate(Radius(totalnumparticles))
  allocate(thermalcond(totalnumparticles))
  allocate(Q1(totalnumparticles))
  allocate(Q2(totalnumparticles))
  allocate(Q3(totalnumparticles))
  allocate(Q4(totalnumparticles))
  allocate(H(totalnumparticles))
  allocate(Volume(totalnumparticles))
  allocate(ThetaNew(totalnumparticles))
  allocate(ThetaOld(totalnumparticles))
  allocate(mass(totalnumparticles))
  allocate(density(totalnumparticles))
  allocate(HeatCap(totalnumparticles))
  allocate(Y1(totalnumparticles))
  allocate(Y2(totalnumparticles))
  allocate(Y3(totalnumparticles))
  allocate(Y4(totalnumparticles))
  allocate(RandNumbers(totalnumparticles))
  allocate(OrigThermalCond(totalnumparticles))
  allocate(OrigHeatCap(totalnumparticles))
  allocate(Melted(totalnumparticles))
  allocate(Vaporized(totalnumparticles))
  allocate(AverageTemp(500000))
  allocate(PercentMelted(500000))
  allocate(Particle1(500000))
  allocate(PercentVaporized(500000))

  allocate(PositionOld(numparticles,3))
  allocate(PositionNew(numparticles,3))
  allocate(VOld(numparticles,3))
  allocate(VNew(numparticles,3))
  allocate(PsiTot1(numparticles,3))
  allocate(PsiTot2(numparticles,3))
  allocate(PsiTot3(numparticles,3))
  allocate(PsiTot4(numparticles,3))
  allocate(YF1(numparticles,6))
  allocate(YF2(numparticles,6))
  allocate(YF3(numparticles,6))
  allocate(YF4(numparticles,6))
  allocate(EYoungPart(numparticles))



  ! Reads in x, y, z, radius of each particle
  do i=1,totalnumparticles
     read(unit=110,fmt='(F13.8,2x,F13.8,2x,F13.8,2x,F13.8)') Particles(i,1), Particles(i,2), Particles(i,3), Particles(i,4)
  end do

  close(unit=110)

  !Convert coordinates from mm to m
  do i = 1,numparticles
     Particles(i,1) = Particles(i,1)/1.0e3
     Particles(i,2) = Particles(i,2)/1.0e3
     Particles(i,3) = Particles(i,3)/1.0e3
     Particles(i,4) = Particles(i,4)/1.0e3
  end do

  Radius(:) = Particles(:,4)
  PositionNew(:,:) = Particles(:,1:3)
  PositionOld(:,:) = Particles(:,1:3)
  avgradius = sum(Radius) / numparticles
  avgheight = sum(Particles(:,3)) / numparticles
  maxheight = maxval(Particles(:,3))

  tstart = 0.0D0
  tend = 1.50D0

  t = tstart
  tstep = 5.0e-8
  counter = 1
  k = 1

  AverageTemp(:) = 0.0D0
  ThetaStart = 298.0D0 ! Start temp in K (room temp.)
  ThetaOld(:) = ThetaStart
  ThetaNew(:) = ThetaStart
  absorptivity = 0.44  ! approximate value for 316LSS from Gusarov(2005) and Khairallah(2014) - range is from 0.33 - 0.44
  LaserPower = absorptivity*300.0 ! Watts 
  LaserSpeed = 0.015D0 ! m/s
  LaserDiameter = 1.0e-3 ! m
  porosity = 0.43  ! powder bed porosity of PowderDepositionShankarMPFinal2800.dat
  alpha = 3.0*(1-porosity)/(4.0*porosity*avgRadius) ! [m^-1] extinction coefficient for Beer-Lambert Law
!!$  I0 = 2.0*LaserPower/(PI*LaserDiameter**2.0)
  I0 = LaserPower/(PI*LaserDiameter**2.0) ! W/m^2 for top hat laser
  HatchSpace = LaserDiameter
  LaserHeight = maxheight - 1.2*avgRadius ! height below which laser beam intensity begins diminishing [m]
!!$  LaserHeight = avgHeight 
  print *, LaserHeight

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
  do i = 1,totalnumparticles
     particlebin(i,1) = ceiling(Particles(i,1)/binlength) ! assign bin/ grid coordinates to each particle
     particlebin(i,2) = ceiling(Particles(i,2)/binwidth) 
     particlebin(i,3) = ceiling(Particles(i,3)/binheight)
  end do

  !create bin grid and linked list according to Alejandro's method
  bins(:,:,:) = 0
  linkedbinlist(:) = 0
  do i = totalnumparticles,1,-1
     linkedbinlist(i) = bins(particlebin(i,1),particlebin(i,2),particlebin(i,3))
     bins(particlebin(i,1),particlebin(i,2),particlebin(i,3)) = i
  end do




  ! Particle Material properties for 316L SS at 366 K
  rho = 7919.0 ! kg/m^3
  EYoungStart = 1.93e11 ! Pa
  EYoungLiquid = 1700.0 ! Pa, based off cooling calculations on 11/13/14 for cooling rate of 10^5 K/s
  Tmelt = 1700.0 ! K 
  Tvapor = 3130.0 ! K (guess based on iron)
  thermalcondsol = 14.28 ! W/m-k at 366 K
  thermalcondliq = 32.4 ! at 1700 K according to Hodge et al (Implementation of a Thermomechanical model for the simulation of SLM, 2013)
  thermalcondgas = thermalcondliq*0.01 ! approx thermal conductivity of air at high temps (0.3 W/m-K)
  Csolid = 485.34 ! J/kg*K at 373 K
  Cliquid = 815.0 ! from Gusarov 2007
  Cgas = Cliquid*1.5 ! generally higher than for liquids - need more justification for this
  LatentHeatMelt = 2.99e5 ! J/kg from Gusarov 2007
  LatentHeatVap = 6.09e6 ! J/kg (guess based on iron)
  LatentHeatSolid = 2.99e5 ! J/kg 
  !CTE = 13.82E-6
  EYoungPart(:) = EYoungstart
  deltaTemp = 1.0D0

  ! Substrate Material Properties (currently assumed to be 316L at 366 K)
  ThermalCondSubstrate = 14.28 ! W/m-K at 366 K
  ThetaSubstrate = ThetaStart ! 363 K or 90 C
  EYoungSubstrate = EYoungstart 



  ! Assign particle properties, with some randomly varying around the mean values listed above
  Call init_Random_Seed()
  Call Random_Number(RandNumbers)

  do i = 1,totalnumparticles
     Volume(i) = 4.0/3.0*PI*Radius(i)**3.0 ! m^3
     Mass(i) = Volume(i)*rho ! kg
     HeatCap(i) = Csolid !- 100 + RandNumbers(i)*200
     ThermalCond(i) = thermalcondsol !-10 + RandNumbers(i)*20 
  end do

  OrigThermalCond(:) = ThermalCond(:)
  OrigHeatCap(:) = HeatCap(:)

  Q1(:) = 0.0D0
  H(:) = 0.0D0


!!$!!! INITIALIZE FDM PARAMETERS
!!$  lengthFD = length ! [m]
!!$  widthFD = width
!!$  heightFD = 0.00050D0
!!$
!!$  nodesx = 201
!!$  nodesy = 201
!!$  nodesz = 101
!!$  nodestotal = nodesx*nodesy*nodesz
!!$
!!$  hx = lengthFD / (nodesx-1.0)
!!$  hy = widthFD / (nodesy-1.0)
!!$  hz = heightFD / (nodesz-1.0)
!!$
!!$
!!$  allocate(DensityFD(nodesx,nodesy,nodesz))
!!$  allocate(HeatCapFD(nodesx,nodesy,nodesz))
!!$  allocate(ThermalCondFD(nodesx,nodesy,nodesz))
!!$  allocate(U(nodesx,nodesy,nodesz))
!!$  allocate(ThetaFD(nodesx,nodesy,nodesz))
!!$  allocate(ThetaOldFD(nodesx,nodesy,nodesz))
!!$  allocate(ThetaDotFD(nodesx,nodesy,nodesz))
!!$  allocate(StateFD(nodesx,nodesy,nodesz))
!!$  allocate(GhostNodesTop(nodesx,nodesy))
!!$  allocate(FluxBC(nodesx,nodesy))
!!$  allocate(DivGradX1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(DivGradY1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(DivGradZ1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(CoordsFD(nodestotal,3))
!!$  allocate(AvgT(500000))
!!$  allocate(AvgT1(500000))
!!$  allocate(AvgT2(200000))
!!$  allocate(AvgT3(200000))
!!$  allocate(AvgT4(200000))
!!$  allocate(AvgT0(500000))
!!$  allocate(AvgT5(200000))
!!$  allocate(AvgT1Laser(200000))
!!$  allocate(AvgT2Laser(200000))
!!$  allocate(AvgT3Laser(200000))
!!$  allocate(AvgT4Laser(200000))
!!$  allocate(time(500000))


  ! Parameters for 316L at 366 K, (90 C or 363 K is a typical preheat temp.)
  densitysolid = 7919.0 ! kg/m^3
  densityliquid = 7300.0 ! kg/m^3


!!$  StateFD(:,:,:) = 1 ! 1 - solid, 2 - liquid, 3 - gas
!!$  DensityFD(:,:,:) = densitysolid ! [kg/m^3] initial value at all nodes
!!$  HeatCapFD(:,:,:) = Csolid ! [J/kg-K]
!!$  ThermalCondFD(:,:,:) = thermalcondsol ! [W/m-K]
!!$  ExtinctionCoefFD = 4.0*pi*10.0/1.08e-6  ! extinction coef for Aluminum (approx) under fiber laser - (p.50 of J.R. Master Thesis), very large number that will ensure all laser energy is absorbed in top nodes of FD mesh
!!$  U0 = I0*exp(-alpha*(LaserHeight-0.0)) / hz ! Volumetric heat source from laser at bottom of particles or top of substrate (W/m^3)
!!$
!!$  Theta0 = ThetaStart ! [K] Initial temp and B.C. on faces of cube
!!$
!!$  do i=1,numtimesteps+1
!!$     time(i) = tstart + (i-1)*tstep
!!$  end do
!!$
!!$  ThetaOldFD(:,:,:) = Theta0
!!$  ThetaFD(:,:,:) = Theta0
!!$  ThetaDotFD(:,:,:) = 0.0
!!$  time(:) = 0.0
!!$  GhostNodesTop(:,:) = Theta0
!!$  FluxBC(:,:) = 0.0D0



  !!Let us not seek to satisfy our thirst for freedom by drinking from the cup of bitterness and hatred

  !!BEGIN TIME LOOP
  do n = 1,layers  !! n represents the layers as they get stacked on top of each other

     numparticles = n*orignumparticles  ! assumes the same number of particles are in each layer

!!$   LaserPosition(1) = -LaserDiameter/2.0
!!$   LaserPosition(2) = LaserDiameter/2.0
     LaserPosition(1) = -LaserDiameter/2.0
     LaserPosition(2) = width/2.0
     LaserPosition(2) = 0.0
     t = tstart
     k = 1
     counter = 1
     laserpass = 1
     iterationDEM = 1
     iterationFD = 1

!!$   do while (LaserPosition(2) < width)
!!$     do while (LaserPosition(1) <= (length + LaserDiameter/2.0))

     do while (t <= tend)

!!! Define laser path to move in a straight path once through the middle of the domain
        LaserPosition(1) = LaserPosition(1) + LaserSpeed*tstep

!!$      if (LaserPosition(1) > (length + LaserDiameter/2.0)) then
!!$         LaserPosition(1) = length + LaserDiameter/2.0
!!$         LaserPosition(2) = LaserPosition(2) + HatchSpace
!!$         LaserSpeed = -LaserSpeed
!!$      elseif (LaserPosition(1) < (-LaserDiameter/2.0)) then
!!$         LaserPosition(1) = -LaserDiameter/2.0
!!$         LaserPosition(2) = LaserPosition(2) + HatchSpace
!!$         LaserSpeed = - LaserSpeed
!!$      end if


!!!!!!!!!!!!!!!!!!!!!!!!! DEM PORTION OF CODE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Determine which particles lie within laser zone and receive heat from it (Top Hat Laser)
!$OMP PARALLEL DO
        do i = 1,numparticles
           dist = ((Particles(i,1)-LaserPosition(1))**2.0 + (Particles(i,2)-LaserPosition(2))**2.0)**(0.5)

           if (dist <= (LaserDiameter/2.0)) then
              Inside(i) = 1

!!!! USE BEER LAMBERT LAW TO DETERMINE HEATING AS A FUNCTION OF DEPTH
              if (Particles(i,3) > LaserHeight) then
                 H(i) = I0*PI*Radius(i)**2.0
              else
                 H(i) = I0*exp(-alpha*(LaserHeight-Particles(i,3)))*PI*Radius(i)**2.0
              end if

           else
              Inside(i) = 0
              H(i) = 0.0D0
           end if

        end do
!$OMP END PARALLEL DO

!$OMP BARRIER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$!!!! USE BEER LAMBERT LAW TO DETERMINE HEATING AS A FUNCTION OF DEPTH AND TOP HAT LASER
!!$           if (Particles(i,3) > LaserHeight) then
!!$              H(i) = I0
!!$           else
!!$              H(i) = I0*exp(-alpha*(LaserHeight-Particles(i,3)))
!!$           end if
!!$
!!$        end do


!!!!!!!!!!!!!!!!!!!!!
!!$      do i = 1,numparticles
!!$         do j=1,numbinsx
!!$            do k=1,numbinsy
!!$               NEED TO SEARCH IN particlebin(j,k,1) FOR maxval(Particles(:,3))
!!$               LaserHeight = maxval(Particles(:,3)
!!$               LaserHeight(particlebin(i,1),particlebin(i,2))
!!!!!!!!!!!!!!!!!!!!!!!



!!!!! UPDATE POSITION AND VELOCITY USING EXPLICIT RK-4 METHOD
        call Force(numparticles,PositionOld,VOld,Radius,Mass,height,length,width,PsiTot1,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),EYoungPart)

        YF1(:,1) = Vold(:,1)
        YF1(:,2) = Vold(:,2)
        YF1(:,3) = Vold(:,3)

        YF2(:,1) = Vold(:,1) + tstep/2.0*1.0/mass(:)*Psitot1(:,1)
        YF2(:,2) = Vold(:,2) + tstep/2.0*1.0/mass(:)*Psitot1(:,2)
        YF2(:,3) = Vold(:,3) + tstep/2.0*1.0/mass(:)*Psitot1(:,3)


        call Force(numparticles,PositionOld,YF2(:,1:3),Radius,Mass,height,length,width,PsiTot2,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),EYoungPart) !Psi(Vold=YF2(:,1:3),t=t+tstep/2)

        YF3(:,1) = Vold(:,1) + tstep/2.0*1.0/mass(:)*Psitot2(:,1)  
        YF3(:,2) = Vold(:,2) + tstep/2.0*1.0/mass(:)*Psitot2(:,2)  
        YF3(:,3) = Vold(:,3) + tstep/2.0*1.0/mass(:)*Psitot2(:,3) 


        call Force(numparticles,PositionOld,YF3(:,1:3),Radius,Mass,height,length,width,PsiTot3,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),EYoungPart) !Psi(Vold=YF3(:,1:3),t=t+tstep/2)

        YF4(:,1) = Vold(:,1) + tstep*1.0/mass(:)*Psitot3(:,1)  
        YF4(:,2) = Vold(:,2) + tstep*1.0/mass(:)*Psitot3(:,2)  
        YF4(:,3) = Vold(:,3) + tstep*1.0/mass(:)*Psitot3(:,3)  


        call Force(numparticles,PositionOld,YF4(:,1:3),Radius,Mass,height,length,width,PsiTot4,numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),EYoungPart) !Psi(Vold=YF4(:,1:3),t=t+tstep)

        Vnew(:,1) = Vold(:,1) + tstep/6.0*(1.0/mass(:)*Psitot1(:,1) + 2.0/mass(:)*Psitot2(:,1) + 2.0/mass(:)*Psitot3(:,1) + 1.0/mass(:)*Psitot4(:,1)) 
        Vnew(:,2) = Vold(:,2) + tstep/6.0*(1.0/mass(:)*Psitot1(:,2) + 2.0/mass(:)*Psitot2(:,2) + 2.0/mass(:)*Psitot3(:,2) + 1.0/mass(:)*Psitot4(:,2)) 
        Vnew(:,3) = Vold(:,3) + tstep/6.0*(1.0/mass(:)*Psitot1(:,3) + 2.0/mass(:)*Psitot2(:,3) + 2.0/mass(:)*Psitot3(:,3) + 1.0/mass(:)*Psitot4(:,3)) 

        YF1(:,4) = PositionOld(:,1)
        YF2(:,4) = PositionOld(:,1) + tstep/2.0*YF1(:,1)
        YF3(:,4) = PositionOld(:,1) + tstep/2.0*YF2(:,1)
        YF4(:,4) = PositionOld(:,1) + tstep*YF3(:,1)

        YF1(:,5) = PositionOld(:,2)
        YF2(:,5) = PositionOld(:,2) + tstep/2.0*YF1(:,2)
        YF3(:,5) = PositionOld(:,2) + tstep/2.0*YF2(:,2)
        YF4(:,5) = PositionOld(:,2) + tstep*YF3(:,2)

        YF1(:,6) = PositionOld(:,3)
        YF2(:,6) = PositionOld(:,3) + tstep/2.0*YF1(:,3)
        YF3(:,6) = PositionOld(:,3) + tstep/2.0*YF2(:,3)
        YF4(:,6) = PositionOld(:,3) + tstep*YF3(:,3)

        PositionNew(:,1) = PositionOld(:,1) + tstep/6.0*(YF1(:,1) + 2.0*YF2(:,1) + 2.0*YF3(:,1) + YF4(:,1))
        PositionNew(:,2) = PositionOld(:,2) + tstep/6.0*(YF1(:,2) + 2.0*YF2(:,2) + 2.0*YF3(:,2) + YF4(:,2))
        PositionNew(:,3) = PositionOld(:,3) + tstep/6.0*(YF1(:,3) + 2.0*YF2(:,3) + 2.0*YF3(:,3) + YF4(:,3))



!!!CALCULATE NEW PARTICLE TEMPERATURE FORWARD EULER
        call HeatFlux(numparticles,Particles(1:numparticles,:),ThetaOld(1:numparticles),ThermalCond(1:numparticles),Q1(1:numparticles),numbinsx,numbinsy,numbinsz,bins,linkedbinlist(1:numparticles),ParticleBin(1:numparticles,:),LaserHeight)

        ThetaNew(:) = ThetaOld(:) + tstep/(mass(:)*HeatCap(:))*(Q1(:)+H(:))




!!!!UPDATE THERMALLY DEPENDENT MATERIAL PROPERTIES BEFORE MOVING TO NEXT TIME STEP 
!$OMP PARALLEL DO
        do i = 1,numparticles
           if (ThetaNew(i) < 1400.0) then  ! data up til 1255 K, extrapolate up to 1400 K then hold properties constant until Tmelt at 1658 K
              HeatCap(i) = 0.213256*ThetaNew(i) + 407.130411
              ThermalCond(i) = 0.013290*ThetaNew(i) + 9.461666
         !     EYoungPart(i) = (-0.082478*ThetaNew(i) + 223.898672)*1.0e9
              !            CTE = (0.003292*ThetaNew(i) + 12.526519)*1e-6
              !            Radius(i) = (Radius(i)**3/(1-3*CTE*(ThetaNew(i)-ThetaStart)))**(1.0/3.0)
           elseif (ThetaNew(i) >= 1400.0 .AND. ThetaNew(i) < (Tmelt-deltaTemp)) then
              HeatCap(i) = 0.213256*1400.0 + 407.130411
              ThermalCond(i) = 0.013290*1400.0 + 9.461666
            !  EYoungPart(i) = (-0.082478*1400.0 + 223.898672)*1.0e9
              !            CTE = (0.003292*1200.0 + 12.526519)*1e-6
           end if
        end do
!$OMP END PARALLEL DO


!!!ACCOUNT FOR PHASE CHANGES
!$OMP PARALLEL DO
        do i = 1,numparticles
           if (ThetaOld(i) < (Tmelt-deltaTemp) .AND. ThetaNew(i) < (Tmelt-deltaTemp)) then  ! SOLID
              HeatCap(i) = HeatCap(i)
              ThermalCond(i) = ThermalCond(i)  ! Defined to be solid properties
           elseif (ThetaOld(i) < Tmelt .AND. ThetaNew(i) >= Tmelt) then  ! MELTING
              ThetaNew(i) = Tmelt + deltaTemp/100.0
              HeatCap(i) = HeatCap(i) + LatentHeatMelt/deltaTemp
              EYoungPart(i) = minval((/EYoungStart**0.5,EYoungPart(i)/)) ! becoming softer by 1/2 each step
              ThermalCond(i) = ThermalCond(i)
           elseif (ThetaOld(i) >= (Tmelt+deltaTemp) .AND. ThetaNew(i) >= (Tmelt+deltaTemp) .AND. ThetaOld(i) < Tvapor .AND. ThetaNew(i) < Tvapor) then  ! LIQUID
              HeatCap(i) = Cliquid
              EYoungPart(i) = EYoungLiquid
              ThermalCond(i) = thermalcondliq ! Liquid properties
           elseif (ThetaOld(i) >= Tmelt .AND. ThetaOld(i) < Tvapor .AND. ThetaNew(i) >= Tvapor) then  ! BOILING
              ThetaNew(i) = Tvapor + deltaTemp/100.0
              HeatCap(i) = Cliquid + LatentHeatVap/deltaTemp
              ThermalCond(i) = thermalcondliq
           elseif (ThetaOld(i) >=(Tvapor+deltaTemp) .AND. ThetaNew(i) >= (Tvapor+deltaTemp)) then  ! GAS
!!$              ThetaNew(i) = Tvapor + deltaTemp + 1  !!!Constrain Particle Temp
              HeatCap(i) = Cgas
              ThermalCond(i) = thermalcondgas ! Gas properties
           elseif (ThetaOld(i) >= Tmelt .AND. ThetaNew(i) < Tmelt) then  ! SOLIDIFYING
              ThetaNew(i) = Tmelt - deltaTemp/100.0
              HeatCap(i) = Cliquid + LatentHeatSolid/deltaTemp
           elseif (ThetaOld(i) >= Tvapor .AND. ThetaNew(i) < Tvapor) then  ! CONDENSING
              ThetaNew(i) = Tvapor - deltaTemp/100.0
              HeatCap(i) = Cgas + LatentHeatVap/deltaTemp
           end if
        end do
!$OMP END PARALLEL DO

!$OMP BARRIER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$        ! CALCULATE AVERAGE TEMP
!!$        temp = 0.0
!!$        do i = 1,numparticles
!!$           temp = temp + ThetaNew(i)
!!$        end do
!!$!!!$        AverageTemp(counter) = temp/numparticles
!!$!!!$        Particle1(counter) = ThetaNew(1)
!!$
!!$        ! Determine percentage of particles melted and vaporized
!!$        do i = 1,numparticles
!!$
!!$           if (ThetaNew(i) > Tmelt .AND. ThetaNew(i) < (Tvapor+deltaTemp)) then
!!$              Melted(i) = 1.0
!!$           else
!!$              Melted(i) = 0.0
!!$           end if
!!$
!!$           if (ThetaNew(i) > (Tvapor+deltaTemp)) then
!!$              Vaporized(i) = 1.0
!!$           else
!!$              Vaporized(i) = 0.0
!!$           end if
!!$
!!$        end do
!!$!!!$        PercentMelted(counter) = sum(Melted)/numparticles
!!$!!!$        PercentVaporized(counter) = sum(Vaporized)/numparticles



!!!OUTPUT DEM DATA TO TECPLOT every 50000 time steps, 600 total - total simulation time of 1.5 seconds
        if ((counter==50000*iterationDEM) .OR. (counter==1)) then
           write(unit=100,fmt='(A)') 'Zone'
           write(unit=150,fmt='(A)') 'Zone'
           do i = 1,numParticles 
              write(unit=100,fmt=70) 1e3*Particles(i,1), 1e3*Particles(i,2), 1e3*Particles(i,3), 1e3*Particles(i,4), ThetaNew(i), Inside(i)
            
!!$              if (Particles(i,2) > width/2.0) then
!!$                 write(unit=150,fmt=70) 1e3*Particles(i,1), 1e3*Particles(i,2), 1e3*Particles(i,3), 1e3*Particles(i,4), ThetaNew(i), Inside(i)
!!$              end if
           end do
70         format (F13.6,2x,F13.6,2x,F13.6,2x,F13.6,2x,F13.6,2x,I4)
           print *, iterationDEM
           iterationDEM = iterationDEM + 1
        end if





!!!!!!!!!!!!!!!!!!!! FDM PORTION OF CODE  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        AvgT(counter) = sum(sum(sum(ThetaFD,1),1)) / (nodestotal*1.0)
!!$        AvgT0(counter) = sum(sum(ThetaFD(:,:,1),1))/(nodesx*nodesy*1.0);
!!$        AvgT1(counter) = sum(sum(ThetaFD(:,:,11),1))/(nodesx*nodesy*1.0);
!!$     AvgT2(iteration) = sum(sum(ThetaFD(:,:,9),1))/(nodesx*nodesy*1.0);
!!$     AvgT3(iteration) = sum(sum(ThetaFD(:,:,13),1))/(nodesx*nodesy*1.0);
!!$     AvgT4(iteration) = sum(sum(ThetaFD(:,:,17),1))/(nodesx*nodesy*1.0);
!!$     AvgT5(iteration) = sum(sum(ThetaFD(:,:,21),1))/(nodesx*nodesy*1.0);
!!$     AvgT1Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,5),1))/25.0;
!!$     AvgT2Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,9),1))/25.0;    
!!$     AvgT3Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,13),1))/25.0;    
!!$     AvgT4Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,17),1))/25.0; 


      !  FluxBC(:,:) = 0.0

!!$        do k = 1,nodesz-1
!!$           do j = 2,nodesy-1
!!$              do i = 2,nodesx-1
!!$
!!$                 !! Determine laser heat input at each node
!!$                 xcoord = (i-1.0)*hx
!!$                 ycoord = (j-1.0)*hy
!!$                 zcoord = -(k-1.0)*hz
!!$                 dist = ((xcoord-LaserPosition(1))**2.0 + (ycoord-LaserPosition(2))**2.0)**(0.5)
!!$
!!$                 U(i,j,k) = U0*exp(-2.0*dist**2.0/LaserDiameter**2.0)*exp(-ExtinctionCoefFD*(0.0-zcoord))
!!$
!!$                 if (k==1) then
!!$                    ! 2nd order, 3 pt stencil for Neumann BC
!!$                    GhostNodesTop(i,j) = (2*hz*FluxBC(i,j)/ThermalCondFD(i,j,k) + 4.0*ThetaFD(i,j,k) - ThetaFD(i,j,k+1)) / 3.0

                 ! 1st order Neumann BC
        !         GhostNodesTop(i,j) = hz*FluxBC(i,j)/ThermalCondFD(i,j,k) + ThetaFD(i,j,k)
!!$
!!$                    DivGradZ = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k+1))/2.0 * (ThetaFD(i,j,k+1)-ThetaFD(i,j,k))/hz**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k))/2.0 * (ThetaFD(i,j,k)-GhostNodesTop(i,j))/hz**2.0
!!$
!!$                 else
!!$
!!$                    DivGradZ = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k+1))/2.0 * (ThetaFD(i,j,k+1)-ThetaFD(i,j,k))/hz**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k-1))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i,j,k-1))/hz**2.0
!!$
!!$                 end if
!!$
!!$
!!$                 !! Calculate heat conduction using FD method
!!$                 DivGradX = (ThermalCondFD(i,j,k)+ThermalCondFD(i+1,j,k))/2.0 * (ThetaFD(i+1,j,k)-ThetaFD(i,j,k))/hx**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i-1,j,k))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i-1,j,k))/hx**2.0
!!$                 DivGradY = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j+1,k))/2.0 * (ThetaFD(i,j+1,k)-ThetaFD(i,j,k))/hy**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j-1,k))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i,j-1,k))/hy**2.0
!!$
!!$
!!$                 ThetaDotFD(i,j,k) = 1/(DensityFD(i,j,k)*HeatCapFD(i,j,k)) * (DivGradX+DivGradY+DivGradZ + U(i,j,k))
!!$
!!$              end do
!!$           end do
!!$        end do
!!$
!!$
!!$        ! Dirichlet BCs everywhere except for Neumann BC at top face
!!$        ThetaOldFD = ThetaFD
!!$        ThetaFD = ThetaFD + tstep * ThetaDotFD ! FWD EULER 
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!! Update Material Properties
!!$        do i = 1,nodesx
!!$           do j = 1,nodesy
!!$              do k = 1,nodesz
!!$
!!$!!!! UPDATE PROPERTIES, data up til 1255 K, extrapolate up to 1400 K then hold properties constant until Tmelt at 1700 K
!!$                 if (ThetaFD(i,j,k) < 1400.0) then   
!!$                    HeatCapFD(i,j,k) = 0.213256*ThetaFD(i,j,k) + 407.130411
!!$                    ThermalCondFD(i,j,k) = 0.013290*ThetaFD(i,j,k) + 9.461666
!!$                    DensityFD(i,j,k) = -0.428874*ThetaFD(i,j,k) + 8081.792237
!!$                    !            CTE = (0.003292*ThetaNew(i) + 12.526519)*1e-6
!!$                    !            Radius(i) = (Radius(i)**3/(1-3*CTE*(ThetaNew(i)-ThetaStart)))**(1.0/3.0)
!!$                 elseif (ThetaFD(i,j,k) >= 1400.0 .AND. ThetaFD(i,j,k) < (Tmelt-deltaTemp)) then
!!$                    HeatCapFD(i,j,k) = 0.213256*1400.0 + 407.130411
!!$                    ThermalCondFD(i,j,k) = 0.013290*1400.0 + 9.461666
!!$                    DensityFD(i,j,k) = -0.428874*1400.0 + 8081.792237
!!$                    !            CTE = (0.003292*1200.0 + 12.526519)*1e-6
!!$                 end if
!!$
!!$!!!! PHASE CHANGES
!!$                 if (ThetaOldFD(i,j,k) < Tmelt .AND. ThetaFD(i,j,k) < Tmelt) then ! SOLID
!!$                    HeatCapFD(i,j,k) = HeatCapFD(i,j,k)
!!$                    StateFD(i,j,k) = 1
!!$                    ThermalCondFD(i,j,k) = ThermalCondFD(i,j,k) 
!!$                 elseif (ThetaOldFd(i,j,k) < Tmelt .AND. ThetaFD(i,j,k) >= Tmelt) then ! MELTING
!!$                    ThetaFD(i,j,k) = Tmelt + deltaTemp/100.0
!!$                    HeatCapFD(i,j,k) = HeatCapFD(i,j,k) + LatentHeatMelt/deltaTemp
!!$                    !    ThermalCond(i) = OrigThermalCond(i)
!!$                 elseif (ThetaOldFD(i,j,k) >= (Tmelt+deltaTemp) .AND. ThetaFD(i,j,k) >= (Tmelt+deltaTemp) .AND. ThetaOldFD(i,j,k) < Tvapor .AND. ThetaFD(i,j,k) < Tvapor) then ! LIQUID
!!$                    HeatCapFD(i,j,k) = Cliquid
!!$                    DensityFD(i,j,k) = densityliquid
!!$                    StateFD(i,j,k) = 2
!!$                    ThermalCondFD(i,j,k) = thermalcondliq ! Liquid properties
!!$                 elseif (ThetaOldFD(i,j,k) >= Tmelt .AND. ThetaOldFD(i,j,k) < Tvapor .AND. ThetaFD(i,j,k) >= Tvapor) then ! BOILING
!!$                    ThetaFD(i,j,k) = Tvapor + deltaTemp/100.0
!!$                    HeatCapFD(i,j,k) = Cliquid + LatentHeatVap/deltaTemp
!!$                    ThermalCondFD(i,j,k) = 0.5*thermalcondsol
!!$                 elseif (ThetaOldFD(i,j,k) >=(Tvapor+deltaTemp) .AND. ThetaFD(i,j,k) >= (Tvapor+deltaTemp)) then ! GAS
!!$                    ThetaFD(i,j,k) = Tvapor + deltaTemp + 1  !!!Constrain Particle Temp
!!$                    ThermalCondFD(i,j,k) = 0.01*thermalcondsol
!!$                    StateFD(i,j,k) = 3
!!$                 elseif (ThetaOldFD(i,j,k) >= Tmelt .AND. ThetaFD(i,j,k) < Tmelt) then ! SOLIDIFYING
!!$                    ThetaFD(i,j,k) = Tmelt - deltaTemp/100.0
!!$                    HeatCapFD(i,j,k) = Cliquid + LatentHeatSolid/deltaTemp
!!$                    !    ThermalCond(i) = 1.5*OrigThermalCond(i)
!!$                    !         else
!!$                    !            print *, 'heat capacity case missing'
!!$                 end if
!!$
!!$              end do
!!$           end do
!!$        end do
!!$

!!!!!!!!!!!!!!!!!!!!!!!!!!




!!$!!!OUTPUT data TO TECPLOT every 12000 time steps: x, y, z [mm], Temp [K]
!!$        if (counter==1) then
!!$           write(20,*) 'Variables="X","Y","Z","R","T"'
!!$           write(20,*) 'Zone I=',nodesx, 'J=',nodesy, 'K=',nodesz, 'F=POINT'
!!$           do k = 1,nodesz
!!$              do j = 1,nodesy
!!$                 do i = 1,nodesx
!!$                    write(unit=20,fmt=70) (i-1.0)*hx*1.0e3, (j-1.0)*hy*1.0e3, -(k-1.0)*hz*1.0e3, 1.0, ThetaFD(i,j,k)
!!$                 end do
!!$              end do
!!$           end do
!!$           print *, iterationFD
!!$           iterationFD = iterationFD + 1
!!$        elseif (counter==12000*iterationFD) then
!!$            write(20,*) 'Variables="X","Y","Z","R","T"'
!!$           write(20,*) 'Zone I=',nodesx, 'J=',nodesy, 'K=',nodesz, 'F=POINT', 'VARSHARELIST=([1,2,3,4]=1)'
!!$           do k = 1,nodesz
!!$              do j = 1,nodesy
!!$                 do i = 1,nodesx
!!$                    write(unit=20,fmt='(F13.6)') ThetaFD(i,j,k)
!!$                 end do
!!$              end do
!!$           end do
!!$           print *, iterationFD
!!$           iterationFD = iterationFD + 1
!!$        end if



!!! GIVE VARIABLES NEW VALUES
        ThetaOld(:) = ThetaNew(:)
        Particles(:,1:3) = PositionNew(:,:)
        PositionOld(:,:) = PositionNew(:,:)
        Vold(:,:) = Vnew(:,:)

        counter = counter + 1
        t = t + tstep
!!$        time(counter) = time(counter-1) + tstep

     end do



  end do



!!!! Output Final FD Results
!!$  AvgT(counter) = sum(sum(sum(ThetaFD,1),1)) / (nodestotal*1.0)
!!$  AvgT0(counter) = sum(sum(ThetaFD(:,:,1),1))/(nodesx*nodesy*1.0);
!!$  AvgT1(counter) = sum(sum(ThetaFD(:,:,11),1))/(nodesx*nodesy*1.0);

!!$  write(30,*) 'Variables="X","Y","Z","R","T"'
!!$  write(30,*) 'Zone I=',nodesx, 'J=',nodesy, 'K=',nodesz, 'F=POINT'
!!$  do k = 1,nodesz
!!$     do j = 1,nodesy
!!$        do i = 1,nodesx
!!$           write(unit=30,fmt=70) (i-1.0)*hx*1.0e3, (j-1.0)*hy*1.0e3, -(k-1.0)*hz*1.0e3, 1.0, ThetaFD(i,j,k)
!!$        end do
!!$     end do
!!$  end do



!!$  ! Output data
!!$  do i = 1,counter
!!$     write(unit=10,fmt=80) time(i), AvgT(i), AvgT0(i), AvgT1(i)
!!$  end do

!!$  close(unit=10)
!!$  close(unit=20)
!!$  close(unit=30)


!!$  print *, time(counter)
  !print *, minval(minval(LaserHeight,1))
!!$  print *, maxval(maxval(maxval(ThetaFD,1),1))



  ! Output final DEM results to file
  write(unit=200,fmt='(I10,2x,F13.6,2x,F13.6)') numparticles, absorptivity, LaserPower/absorptivity

  do i = 1,numParticles 
     write(unit=200,fmt=70) 1e3*Particles(i,1), 1e3*Particles(i,2), 1e3*Particles(i,3), 1e3*Particles(i,4), ThetaNew(i), Inside(i)
  end do

!!$  do i = 1,counter-1
!!$     write(unit=300,fmt='(F13.6,2x,F13.6)') AverageTemp(i), Particle1(i)
!!$     write(unit=400,fmt='(F13.6,2x,F13.6)') PercentMelted(i), PercentVaporized(i)
!!$  end do

  print *, t

  close(unit=100)
!!$  close(unit=150)
  close(unit=200)
!!$  close(unit=300)
!!$  close(unit=400)

end program SLMShankarFastMP







!!!!! FORCE CALCULATION SUBROUTINE !!!!!!
SUBROUTINE Force(numparticles,PositionOld,VOld,Radius,Mass,height,length,width,PsiTot,numbinsx,numbinsy,numbinsz,bins,linkedbinlist,ParticleBin,EPart)
  implicit none

  integer, intent(in) :: numparticles, numbinsx, numbinsy, numbinsz
  integer, dimension(numparticles), intent(in) :: linkedbinlist
  integer, dimension(numparticles,3), intent(in) :: ParticleBin
  integer, dimension(numbinsx,numbinsy,numbinsz), intent(in) :: bins
  double precision, intent(in) :: height, length, width
  !double precision, dimension(3), intent(in) :: Vwall
  double precision, dimension(numparticles), intent(in) :: Radius, Mass, EPart
  double precision, dimension(numparticles,3), intent(in) ::  VOld, PositionOld
  double precision, dimension(numparticles,3), intent(out) :: PsiTot
  double precision, dimension(numparticles,3) :: CenterCoords, PsiCont, PsiFric, PsiAdh, PsiEnv, PsiGrav
  double precision, dimension(3) :: normal, tangent, psicontwall, psifricwall, psicontpart, psifricpart, psiadhwall, psiadhpart, psienvwall, psienvpart, VTan, Vtan2
  double precision, dimension(numparticles) :: deltawallx1, deltawallx2, deltawally1, deltawally2, deltawallz1, deltawallz2, deltawalloldx1, deltawalloldx2, deltawalloldy1, deltawalloldy2, deltawalloldz1, deltawalloldz2
  double precision, dimension(numparticles,numparticles) :: deltaparticleold, deltaparticle
  double precision :: Bpart, Bwall, MuD, MuS, zdist, xdist, ydist, distance, epsilonwall, epsilonpart, speed, DampCoef, g, NuPart, NuWall, EWall, EstarWall, EstarPart, Rstar, deltawall, deltapart, rhoN2, viscosityN2, Reynolds, deltawalldot, deltapartdot, gamma1, gamma2, mstar, zeta
  integer :: i, j, k, m, n, h1, h2, h3
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462



  !! Parameters
  Bwall = 5.0E-6 !Calculated analytically (see 10/29/13) Deposition
  Bpart = 0.333*Bwall !0.2*Bwall
  DampCoef = 5.0E-8
  MuD = 0.1
  MuS = 0.2
  g = 9.81 ! downward gravitational acceleration [m/s^2]
  NuPart = 0.26D0  ! Poisson's ratio of 316L SS at 422 K
  NuWall = 0.26D0 ! 
 ! EPart = 187.0e9  ! [Pa], Young's Modulus of HX at 373 K
  Ewall = 1.93e11 ! [Pa] 316L SS at 366 K
 ! EstarWall = (EPart*Ewall)/(Ewall*(1-NuPart**2.0)+Epart*(1-NuWall**2.0))  ! Effective Modulus for wall
 ! EstarPart = (EPart*EPart)/(EPart*(1-NuPart**2.0)+EPart*(1-NuPart**2.0))  ! Effective Modulus for particles all made of same material, HX
  rhoN2 = 1.165 ! [kg/m^3] at NTP (20 C, 1 atm) assuming a Nitrogen surrounding atmosphere
  viscosityN2 = 1.75e-5 ! [kg/m-s or Pa*s] Nitrogen viscosity at NTP (20 C, 1 atm), Argon is 2.1e-5
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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
        EstarWall = (EPart(i)*Ewall)/(Ewall*(1-NuPart**2.0)+Epart(i)*(1-NuWall**2.0))  ! Effective Modulus for wall
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


!!!!PARTICLE TO PARTICLE COLLISIONS - WITH BINNING AND OPENMP

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
                       EstarPart = (EPart(i)*EPart(j))/(EPart(j)*(1-NuPart**2.0)+EPart(i)*(1-NuPart**2.0))
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







SUBROUTINE HeatFlux(numparticles,Particles,ThetaOld,ThermalCond,Q,numbinsx,numbinsy,numbinsz,bins,linkedbinlist,ParticleBin,LaserHeight)
  implicit none

  integer, intent(in) :: numparticles, numbinsx, numbinsy, numbinsz !, nodesx, nodesy
  double precision, intent(in) :: LaserHeight !, hx, hy
  double precision, dimension(numparticles,4), intent(in) :: Particles
  double precision, dimension(numparticles), intent(in) :: ThetaOld, ThermalCond
!  double precision, dimension(nodesx,nodesy), intent(in) :: ThetaFDTop, ThermalCondFDTop
  integer, dimension(numparticles), intent(in) :: linkedbinlist
  integer, dimension(numparticles,3), intent(in) :: ParticleBin
  integer, dimension(numbinsx,numbinsy,numbinsz), intent(in) :: bins
  double precision, dimension(numparticles), intent(out) :: Q
 ! double precision, dimension(nodesx,nodesy), intent(out) :: FluxBC
  double precision :: xdist, ydist, zdist, distance, A0, Amax, Acon, avgthermalcond, theta, r, ThetaEnv, emissivity, SB, HTC, avgradius, xcoord, ycoord, ThermalCondSubstrate, ThetaSubstrate, LocalConductionHT
  integer :: i, j, k, m, n, h1, h2, h3, point, count, Xnode, Ynode
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462


!!!NEW VARIABLES TO PASS IN: ParticleBin, numbinsx, numbinsy, numbinsz, bins, linkedbinlist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!HEAT EQUATION FOR LASER!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ThetaEnv = 363 ! [K]
  emissivity = 0.30 ! emissivity at avg temp. (1020 K) for polished 316 SS from Radiative Heat Transfer by Modest
  SB = 5.6704e-8 ! [W/m^2-K^4]
  HTC = 40.0 ! [W/m^2-K] (range of HTC from 20 - 50 W/m^2-K, see 'Nusselt.m' for calculations), using value close to high end since spheres will cause more turbulence than a flat plate
  Q(:) = 0.0D0
!!$  FluxBC(:,:) = 0.0D0
  count = 0
  avgradius = sum(Particles(:,4)) / numparticles
  ThermalcondSubstrate = 14.28 ! W/m-K of 316L SS at 366 K
  ThetaSubstrate = 298 ! substrate fixed at room temp
  

!$OMP PARALLEL DO
  do i = 1,numparticles

!!! PARTICLE TO PARTICLE CONDUCTION
     ! Search through neighboring bins and update heat transfer as necessary
     h1 = particlebin(i,1)
     h2 = particlebin(i,2)
     h3 = particlebin(i,3)

     do k = -1,1,1
        do m = -1,1,1
           do n = -1,1,1
              if ((h1+k > 0) .AND. (h1+k <= numbinsx) .AND. (h2+m > 0) .AND. (h2+m <= numbinsy) .AND. (h3+n > 0) .AND. (h3+n <= numbinsz)) then

                 j = bins(h1+k,h2+m,h3+n)

                 do while ((j/=0) .AND. (j <= numparticles))

                    xdist = Particles(i,1)-Particles(j,1)
                    ydist = Particles(i,2)-Particles(j,2)
                    zdist = Particles(i,3)-Particles(j,3)
                    distance = (xdist**2 + ydist**2 + zdist**2)**(0.5)

                    if ((distance <= (Particles(i,4)+Particles(j,4))) .AND. i/=j) then
                       avgthermalcond = (thermalcond(i)+thermalcond(j)) / 2.0

                       Theta = acos((Particles(j,4)**2.0-Particles(i,4)**2.0-distance**2.0)/(-2.0*distance*Particles(i,4)))
                       r = Particles(i,4)*sin(Theta)
                       Acon = PI*r**2.0

!!$                     if (count < 10) then
!!$                        print *, Acon
!!$                     end if
!!$
!!$                     count = count + 1

!!$                     A0 = PI/8.0*(Particles(i,4)+Particles(j,4))**2.0
!!$                     Amax = PI/4.0*(Particles(i,4)+Particles(j,4))**2.0
!!$                     Acon = A0*(Particles(i,4)+Particles(j,4)) / distance
!!$                     if (Acon > Amax) then
!!$                        Acon = Amax
!!$                     else 
!!$                        Acon = Acon
!!$                     end if

                       Q(i) = Q(i) + Acon*avgthermalcond*(ThetaOld(j)-ThetaOld(i)) / distance

                    end if

                    j = linkedbinlist(j)

                 end do

              end if

           end do
        end do
     end do


!!!PARTICLE TO BOTTOM SUBSTRATE CONDUCTION, SUBSTRATE MODELED VIA FDM, HT USED AS FLUX BC FOR FD MESH
     zdist = abs(Particles(i,3)-0.0)   
     if (zdist < Particles(i,4)) then
!!$        xcoord = Particles(i,1)
!!$        ycoord = Particles(i,2)
!!$        Xnode = floor(xcoord/hx) + 1
!!$        Ynode = floor(ycoord/hy) + 1 
!!$        ThermalcondSubstrate = (ThermalCondFDTop(Xnode,Ynode) + ThermalCondFDTop(Xnode+1,Ynode) + ThermalCondFDTop(Xnode,Ynode+1) + ThermalCondFDTop(Xnode+1,Ynode+1)) / 4.0
!!$        ThetaSubstrate = (ThetaFDTop(Xnode,Ynode) + ThetaFDTop(Xnode+1,Ynode) + ThetaFDTop(Xnode,Ynode+1) + ThetaFDTop(Xnode+1,Ynode+1)) / 4.0

        avgthermalcond = (thermalcond(i)+ThermalCondSubstrate) / 2.0
        Acon = PI*(Particles(i,4)**2.0-zdist**2.0)

!!$      if(Particles(i,4) >= 9.5e-6 .AND. Particles(i,4) <= 10.5e-6) then
!!$         print *, Particles(i,4)
!!$         print *, Acon
!!$      end if
!!$
        LocalConductionHT = Acon*avgthermalcond*(ThetaSubstrate-ThetaOld(i)) / zdist ! [W]
        Q(i) = Q(i) + LocalConductionHT

!!$        ! Average Conductive HT from particle into flux B.C. distributed equally over 4 corresponding FD nodes
!!$        FluxBC(Xnode,Ynode) = FluxBC(Xnode,Ynode) - LocalConductionHT / (hx*hy*4.0)    ! [W/m^2]
!!$        FluxBC(Xnode+1,Ynode) = FluxBC(Xnode+1,Ynode) - LocalConductionHT / (hx*hy*4.0) 
!!$        FluxBC(Xnode,Ynode+1) = FluxBC(Xnode,Ynode+1) - LocalConductionHT / (hx*hy*4.0) 
!!$        FluxBC(Xnode+1,Ynode+1) = FluxBC(Xnode+1,Ynode+1) - LocalConductionHT / (hx*hy*4.0) 

     end if


!!!CONVECTION + RADIATION FROM TOP LAYER OF PARTICLES 
     if (Particles(i,3) >= (LaserHeight-avgRadius)) then
        Q(i) = Q(i) + HTC*(ThetaEnv - ThetaOld(i))*PI*Particles(i,4)**2.0  ! Convection
        Q(i) = Q(i) + emissivity*SB*(ThetaEnv**4.0 - ThetaOld(i)**4.0)*PI*Particles(i,4)**2.0  ! Radiation
     end if


  end do
!$OMP END PARALLEL DO

!$OMP BARRIER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




END SUBROUTINE HeatFlux





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
END SUBROUTINE init_random_seed






!!$SUBROUTINE FDSubstrate
!!$
!!$
!!$  double precision :: tstep, tstart, tend, lengthFD, widthFD, heightFD, hx, hy, hz, DivGradX, DivGradY, DivGradZ, absorptivity, ExtinctionCoefFD, I0, Theta0, U0, LaserPower, LaserSpeed, LaserDiameter, densitysolid, Tmelt, Tvapor, thermalcondsol, thermalcondliq, Csolid, Cliquid, LatentHeatMelt, LatentHeatVap, LatentHeatSolid, deltaTemp, CTE, dist, ycoord, xcoord, zcoord, densityliquid, LaserHeight
!!$  double precision, dimension(:,:,:), allocatable :: DensityFD, HeatCapFD, ThermalCondFD, U, ThetaFD, ThetaOldFD, ThetaDotFD, DivGradX1, DivGradY1, DivGradZ1 
!!$  integer, dimension(:,:,:), allocatable :: stateFD
!!$  double precision, dimension(:,:), allocatable :: CoordsFD, GhostNodesTop, FluxBC
!!$  double precision, dimension(:), allocatable :: AvgT, AvgT1, AvgT2, AvgT3, AvgT4, AvgT1Laser, AvgT2Laser, AvgT3Laser, AvgT4Laser, AvgT0, AvgT5, time, LaserPosition
!!$  integer :: i, j, k, iteration, counter, numtimesteps, nodesx, nodesy, nodesz, nodestotal
!!$  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793238462
!!$
!!$
!!$
!!$  !!OPEN TECPLOT OUTPUT FILE
!!$  Open(unit=100, file='SubstrateMatlabResults.dat')
!!$  write(unit=100,fmt='(A)') 'Title="ME201 Project 5 Results"'
!!$  write(unit=100,fmt='(A)') 'Variables="time","AvgT","AvgT1","AvgT2","AvgT3","AvgT4","AvgT1Laser","AvgT2Laser","AvgT3Laser","AvgT4Laser","AvgT0","AvgT5"'
!!$
!!$  Open(unit=200, file='SubstrateFDResults.dat')
!!$  write(unit=200,fmt='(A)') 'Title="FD Temperature Results"'
!!$80 format (F13.6,2x,F13.6,2x,F13.6,2x,F13.6)
!!$
!!$  Open(unit=300, file='SubstrateFDResultsFinal.dat')
!!$  write(unit=300,fmt='(A)') 'Title="Final FD Temperature Results"'
!!$
!!$
!!$  tstart = 0.0D0
!!$!!!$  tend = 10.0D0
!!$  tstep = 1.0e-6*0.25
!!$!!!$  numtimesteps = ceiling((tend-tstart) / tstep)
!!$
!!$  lengthFD = 0.0010D0 ! [m]
!!$  widthFD = 0.0010D0
!!$  heightFD = 0.00050D0
!!$
!!$  nodesx = 301
!!$  nodesy = 301
!!$  nodesz = 151
!!$  nodestotal = nodesx*nodesy*nodesz
!!$
!!$  hx = lengthFD / (nodesx-1.0)
!!$  hy = widthFD / (nodesy-1.0)
!!$  hz = heightFD / (nodesz-1.0)
!!$
!!$
!!$  allocate(DensityFD(nodesx,nodesy,nodesz))
!!$  allocate(HeatCapFD(nodesx,nodesy,nodesz))
!!$  allocate(ThermalCondFD(nodesx,nodesy,nodesz))
!!$  allocate(U(nodesx,nodesy,nodesz))
!!$  allocate(ThetaFD(nodesx,nodesy,nodesz))
!!$  allocate(ThetaOldFD(nodesx,nodesy,nodesz))
!!$  allocate(ThetaDotFD(nodesx,nodesy,nodesz))
!!$  allocate(StateFD(nodesx,nodesy,nodesz))
!!$  allocate(GhostNodesTop(nodesx,nodesy))
!!$  allocate(FluxBC(nodesx,nodesy))
!!$  allocate(DivGradX1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(DivGradY1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(DivGradZ1(nodesx-2,nodesy-2,nodesz-2))
!!$  allocate(CoordsFD(nodestotal,3))
!!$  allocate(LaserPosition(2))
!!$  allocate(AvgT(200000))
!!$  allocate(AvgT1(200000))
!!$!!!$  allocate(AvgT2(200000))
!!$!!!$  allocate(AvgT3(200000))
!!$!!!$  allocate(AvgT4(200000))
!!$  allocate(AvgT0(200000))
!!$!!!$  allocate(AvgT5(200000))
!!$!!!$  allocate(AvgT1Laser(200000))
!!$!!!$  allocate(AvgT2Laser(200000))
!!$!!!$  allocate(AvgT3Laser(200000))
!!$!!!$  allocate(AvgT4Laser(200000))
!!$  allocate(time(200000))
!!$
!!$
!!$  ! Parameters for 316L at 366 K, (90 C or 363 K is a typical preheat temp.)
!!$  densitysolid = 7919.0 ! kg/m^3
!!$  densityliquid = 7300.0 ! kg/m^3
!!$  Tmelt = 1700.0 ! K (commonly used value in literature, i.e. Gusarov 2009) 
!!$  Tvapor = 3130.0 ! K (guess based on iron)
!!$  thermalcondsol = 14.28 ! W/m-k at 366 K
!!$  ! thermalcondliq = 1.5*thermalcondsol
!!$  Csolid = 485.34 ! J/kg*K at 373 K
!!$  Cliquid = 815.0 ! from Gusarov 2007
!!$  LatentHeatMelt = 2.99e5 ! J/kg from Gusarov 2007
!!$  LatentHeatVap = 6.09e6 ! J/kg (guess based on iron)
!!$  LatentHeatSolid = 2.99e5 ! J/kg 
!!$  !CTE = 13.82E-6
!!$  deltaTemp = 1.0D0
!!$
!!$
!!$  StateFD(:,:,:) = 1 ! 1 - solid, 2 - liquid, 3 - gas
!!$  DensityFD(:,:,:) = densitysolid ! [kg/m^3] initial value at all nodes
!!$  HeatCapFD(:,:,:) = Csolid ! [J/kg-K]
!!$  ThermalCondFD(:,:,:) = thermalcondsol ! [W/m-K]
!!$  ExtinctionCoefFD = 4.0*pi*10.0/1.08e-6  ! extinction coef for Aluminum (approx) under fiber laser - (p.50 of J.R. Master Thesis), very large number that will ensure all laser energy is absorbed in top nodes of FD mesh
!!$  absorptivity = 0.44  ! approximate value for 316LSS from Gusarov(2005) and Khairallah(2014) - range is from 0.33 - 0.44
!!$  LaserPower = absorptivity*45.0 ! Watts 
!!$  LaserSpeed = 0.30D0 ! m/s
!!$  LaserDiameter = 51.0e-6 ! m
!!$  LaserHeight = 0.0D0
!!$  I0 = 2.0*LaserPower/(PI*LaserDiameter**2.0) ! Guassian Laser Beam (W/m^2)
!!$  U0 = I0 / hz ! Volumetric heat source from laser (W/m^3)
!!$
!!$  Theta0 = 363.0 ! [K] Initial temp and B.C. on faces of cube
!!$
!  do i=1,numtimesteps+1
!    time(i) = tstart + (i-1)*tstep
!end do
!!$
!!$  ThetaOldFD(:,:,:) = Theta0
!!$  ThetaFD(:,:,:) = Theta0
!!$  ThetaDotFD(:,:,:) = 0.0
!!$  time(:) = 0.0
!!$  GhostNodesTop(:,:) = Theta0
!!$  FluxBC(:,:) = 0.0D0
!!$
!!$  LaserPosition(1) = -LaserDiameter*3.0
!!$  LaserPosition(2) = widthFD/2.0
!!$  counter = 1
!!$  iteration = 1
!!$
!!$
!!$  ! Begin time loop
!!$  do while (LaserPosition(1) <= (lengthFD + LaserDiameter*3.0))
!!$
!!$!!! Define laser path to move in a straight path once through the middle of the domain
!!$     LaserPosition(1) = LaserPosition(1) + LaserSpeed*tstep
!!$
!!$
!!$     AvgT(iteration) = sum(sum(sum(ThetaFD,1),1)) / (nodestotal*1.0)
!!$     AvgT0(iteration) = sum(sum(ThetaFD(:,:,1),1))/(nodesx*nodesy*1.0);
!!$     AvgT1(iteration) = sum(sum(ThetaFD(:,:,11),1))/(nodesx*nodesy*1.0);
!    AvgT2(iteration) = sum(sum(ThetaFD(:,:,9),1))/(nodesx*nodesy*1.0);
!   AvgT3(iteration) = sum(sum(ThetaFD(:,:,13),1))/(nodesx*nodesy*1.0);
!  AvgT4(iteration) = sum(sum(ThetaFD(:,:,17),1))/(nodesx*nodesy*1.0);
! AvgT5(iteration) = sum(sum(ThetaFD(:,:,21),1))/(nodesx*nodesy*1.0);
!AvgT1Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,5),1))/25.0;
!     AvgT2Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,9),1))/25.0;    
!    AvgT3Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,13),1))/25.0;    
!   AvgT4Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,17),1))/25.0; 
!!$
!!$
!!$
!!$     do k = 1,nodesz-1
!!$        do j = 2,nodesy-1
!!$           do i = 2,nodesx-1
!!$
!!$              !! Determine laser heat input at each node
!!$              xcoord = (i-1.0)*hx
!!$              ycoord = (j-1.0)*hy
!!$              zcoord = -(k-1.0)*hz
!!$              dist = ((xcoord-LaserPosition(1))**2.0 + (ycoord-LaserPosition(2))**2.0)**(0.5)
!!$
!!$              U(i,j,k) = U0*exp(-2.0*dist**2.0/LaserDiameter**2.0)*exp(-ExtinctionCoefFD*(LaserHeight(i,j)-zcoord))
!!$
!!$              if (k==1) then
!!$                 ! 2nd order, 3 pt stencil for Neumann BC
!!$                 GhostNodesTop(i,j) = (2*hz*FluxBC(i,j)/ThermalCondFD(i,j,k) + 4.0*ThetaFD(i,j,k) - ThetaFD(i,j,k+1)) / 3.0

! 1st order Neumann BC
!               GhostNodesTop(i,j) = hz*FluxBC(i,j)/ThermalCondFD(i,j,k) + ThetaFD(i,j,k)
!!$
!!$                 DivGradZ = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k+1))/2.0 * (ThetaFD(i,j,k+1)-ThetaFD(i,j,k))/hz**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k))/2.0 * (ThetaFD(i,j,k)-GhostNodesTop(i,j))/hz**2.0
!!$
!!$              else
!!$
!!$                 DivGradZ = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k+1))/2.0 * (ThetaFD(i,j,k+1)-ThetaFD(i,j,k))/hz**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j,k-1))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i,j,k-1))/hz**2.0
!!$
!!$              end if
!!$
!!$
!!$
!!$              !! Calculate heat conduction using FD method
!!$              DivGradX = (ThermalCondFD(i,j,k)+ThermalCondFD(i+1,j,k))/2.0 * (ThetaFD(i+1,j,k)-ThetaFD(i,j,k))/hx**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i-1,j,k))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i-1,j,k))/hx**2.0
!!$              DivGradY = (ThermalCondFD(i,j,k)+ThermalCondFD(i,j+1,k))/2.0 * (ThetaFD(i,j+1,k)-ThetaFD(i,j,k))/hy**2.0 - (ThermalCondFD(i,j,k)+ThermalCondFD(i,j-1,k))/2.0 * (ThetaFD(i,j,k)-ThetaFD(i,j-1,k))/hy**2.0
!!$
!!$
!!$              ThetaDotFD(i,j,k) = 1/(DensityFD(i,j,k)*HeatCapFD(i,j,k)) * (DivGradX+DivGradY+DivGradZ + U(i,j,k))
!!$
!!$
!!$
!!$           end do
!!$        end do
!!$     end do
!!$
!!$
!!$     ! Dirichlet BCs everywhere except for Neumann BC at top face
!!$     ThetaOldFD = ThetaFD
!!$     ThetaFD = ThetaFD + tstep * ThetaDotFD ! FWD EULER FOR NOW - LATER CHANGE TO RK4
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!
!!$!!!! Update Material Properties
!!$     do i = 1,nodesx
!!$        do j = 1,nodesy
!!$           do k = 1,nodesz
!!$
!!$!!!! UPDATE PROPERTIES, data up til 1255 K, extrapolate up to 1400 K then hold properties constant until Tmelt at 1700 K
!!$              if (ThetaFD(i,j,k) < 1400.0) then   
!!$                 HeatCapFD(i,j,k) = 0.213256*ThetaFD(i,j,k) + 407.130411
!!$                 ThermalCondFD(i,j,k) = 0.013290*ThetaFD(i,j,k) + 9.461666
!!$                 DensityFD(i,j,k) = -0.428874*ThetaFD(i,j,k) + 8081.792237
!!$                 !            CTE = (0.003292*ThetaNew(i) + 12.526519)*1e-6
!!$                 !            Radius(i) = (Radius(i)**3/(1-3*CTE*(ThetaNew(i)-ThetaStart)))**(1.0/3.0)
!!$              elseif (ThetaFD(i,j,k) >= 1400.0 .AND. ThetaFD(i,j,k) < (Tmelt-deltaTemp)) then
!!$                 HeatCapFD(i,j,k) = 0.213256*1400.0 + 407.130411
!!$                 ThermalCondFD(i,j,k) = 0.013290*1400.0 + 9.461666
!!$                 DensityFD(i,j,k) = -0.428874*1400.0 + 8081.792237
!!$                 !            CTE = (0.003292*1200.0 + 12.526519)*1e-6
!!$              end if
!!$
!!$!!!! PHASE CHANGES
!!$              if (ThetaOldFD(i,j,k) < Tmelt .AND. ThetaFD(i,j,k) < Tmelt) then
!!$                 HeatCapFD(i,j,k) = HeatCapFD(i,j,k)
!!$                 StateFD(i,j,k) = 1
!!$                 !  ThermalCond(i) = OrigThermalCond(i)  ! Defined to be solid properties
!!$              elseif (ThetaOldFd(i,j,k) < Tmelt .AND. ThetaFD(i,j,k) >= Tmelt) then
!!$                 ThetaFD(i,j,k) = Tmelt + deltaTemp/100.0
!!$                 HeatCapFD(i,j,k) = HeatCapFD(i,j,k) + LatentHeatMelt/deltaTemp
!!$                 !    ThermalCond(i) = OrigThermalCond(i)
!!$              elseif (ThetaOldFD(i,j,k) >= (Tmelt+deltaTemp) .AND. ThetaFD(i,j,k) >= (Tmelt+deltaTemp) .AND. ThetaOldFD(i,j,k) < Tvapor .AND. ThetaFD(i,j,k) < Tvapor) then
!!$                 HeatCapFD(i,j,k) = Cliquid
!!$                 DensityFD(i,j,k) = densityliquid
!!$                 StateFD(i,j,k) = 2
!!$                 !   ThermalCond(i) = 1.5*OrigThermalCond(i) ! Liquid properties
!!$              elseif (ThetaOldFD(i,j,k) >= Tmelt .AND. ThetaOldFD(i,j,k) < Tvapor .AND. ThetaFD(i,j,k) >= Tvapor) then
!!$                 ThetaFD(i,j,k) = Tvapor + deltaTemp/100.0
!!$                 HeatCapFD(i,j,k) = Cliquid + LatentHeatVap/deltaTemp
!!$                 !           ThermalCond(i) = 0.5*OrigThermalCond(i)
!!$              elseif (ThetaOldFD(i,j,k) >=(Tvapor+deltaTemp) .AND. ThetaFD(i,j,k) >= (Tvapor+deltaTemp)) then
!!$                 ThetaFD(i,j,k) = Tvapor + deltaTemp + 1  !!!Constrain Particle Temp
!!$                 ThermalCondFD(i,j,k) = 0.01*thermalcondsol
!!$                 StateFD(i,j,k) = 3
!!$              elseif (ThetaOldFD(i,j,k) >= Tmelt .AND. ThetaFD(i,j,k) < Tmelt) then
!!$                 ThetaFD(i,j,k) = Tmelt - deltaTemp/100.0
!!$                 HeatCapFD(i,j,k) = Cliquid + LatentHeatSolid/deltaTemp
!!$                 !    ThermalCond(i) = 1.5*OrigThermalCond(i)
!!$                 !         else
!!$                 !            print *, 'heat capacity case missing'
!!$              end if
!!$
!!$           end do
!!$        end do
!!$     end do
!!$
!!$
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$
!!$
!!$
!!$
!!$!!!OUTPUT data TO TECPLOT every 100 time steps: x, y, z [mm], Temp [K]
!!$     if ((iteration==100*counter+1) .OR. (iteration==1)) then
!!$!!!$        write(200,*) 'Variables="X","Y","Z","T"'
!!$!!!$        write(200,*) 'Zone I=',nodesx, 'J=',nodesy, 'K=',nodesz, 'F=POINT'
!!$!!!$        do k = 1,nodesz
!!$!!!$           do j = 1,nodesy
!!$!!!$              do i = 1,nodesx
!!$!!!$                 write(unit=200,fmt=80) (i-1.0)*hx*1.0e3, (j-1.0)*hy*1.0e3, -(k-1.0)*hz*1.0e3, ThetaFD(i,j,k)
!!$!!!$              end do
!!$!!!$           end do
!!$!!!$        end do
!!$        print *, counter
!!$        counter = counter + 1
!!$     end if
!!$
!!$     iteration = iteration + 1
!!$     time(iteration) = time(iteration-1) + tstep
!!$
!!$
!!$  end do
!!$
!!$  AvgT(iteration) = sum(sum(sum(ThetaFD,1),1)) / (nodestotal*1.0)
!!$  AvgT0(iteration) = sum(sum(ThetaFD(:,:,1),1))/(nodesx*nodesy*1.0);
!!$  AvgT1(iteration) = sum(sum(ThetaFD(:,:,11),1))/(nodesx*nodesy*1.0);
! AvgT2(iteration) = sum(sum(ThetaFD(:,:,9),1))/(nodesx*nodesy*1.0);
!AvgT3(iteration) = sum(sum(ThetaFD(:,:,13),1))/(nodesx*nodesy*1.0);
! AvgT4(iteration) = sum(sum(ThetaFD(:,:,17),1))/(nodesx*nodesy*1.0);
!AvgT5(iteration) = sum(sum(ThetaFD(:,:,21),1))/(nodesx*nodesy*1.0);
! AvgT1Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,5),1))/25.0;
! AvgT2Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,9),1))/25.0;    
!AvgT3Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,13),1))/25.0;    
!AvgT4Laser(iteration) = sum(sum(ThetaFD(9:13,9:13,17),1))/25.0; 
!!$
!!$
!!$  write(300,*) 'Variables="X","Y","Z","T"'
!!$  write(300,*) 'Zone I=',nodesx, 'J=',nodesy, 'K=',nodesz, 'F=POINT'
!!$  do k = 1,nodesz
!!$     do j = 1,nodesy
!!$        do i = 1,nodesx
!!$           write(unit=300,fmt=80) (i-1.0)*hx*1.0e3, (j-1.0)*hy*1.0e3, -(k-1.0)*hz*1.0e3, ThetaFD(i,j,k)
!!$        end do
!!$     end do
!!$  end do
!!$
!!$
!!$
!!$  ! Output data
!!$  do i = 1,iteration
!!$     write(unit=100,fmt=80) time(i), AvgT(i), AvgT0(i), AvgT1(i)
!!$  end do
!!$
!!$  close(unit=100)
!!$  close(unit=200)
!!$  close(unit=300)
!!$
!!$
!!$  print *, time(iteration)
!!$  print *, minval(minval(LaserHeight,1))
!!$  print *, maxval(maxval(maxval(ThetaFD,1),1))
!!$
!!$
!!$END SUBROUTINE FDSubstrate
!!$
!!$

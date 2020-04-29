program SpheretoCircle
  implicit none

!!! PROGRAM THAT READS IN PARTICLE COORDS FROM SLM SIMULATIONS AND OUTPUTS A PLANE OF THOSE AS CIRCLES FOR COOLING SIMULATIONS

  double precision, dimension(:,:), allocatable :: Particles
  integer, dimension(:,:), allocatable :: Inside
  integer :: numparticles, i, counter, j, k,  Nx, Ny, Nz, NodesTotal
  double precision :: length, width, height, powderdensity, distance, gridlength, avgheight, avgradius, refinements, middleY, hx, hy, hz, stopheight
  double precision, dimension(:), allocatable :: Radius, Xcoord, Ycoord, Zcoord

  ! Create output file
  open(unit=100,file='InsidePoints.dat')
  write(unit=100,fmt='(A)') 'Title = Points inside vs. outside of the particles'

  ! Read in particle coordinates 
  open(unit=110,file='PowderDepositionShankarMPFinal2800.dat')
  read(unit=110,fmt='(I10)') numparticles

  length = 0.00040D0*1.0e3 ! [m]
  width = 0.00080D0*1.0e3
  height = 0.00350D0*1.0e3
  stopheight = 0.5D0
  middleY = width / 2.0

  !!Allocate and define variables
  allocate(Particles(numparticles,4))
  allocate(Radius(numparticles))

  ! Reads in x, y, z, radius of each particle (in [mm])
  do i=1,numparticles
     read(unit=110,fmt='(F13.8,2x,F13.8,2x,F13.8,2x,F13.8,2x,I10,2x,F13.6)') Particles(i,1), Particles(i,2), Particles(i,3), Particles(i,4)
  end do

  close(unit=110)

  Radius(:) = Particles(:,4)



  ! Create Nodes and generate grid mesh - all coords in [mm]
  refinements = 5.0
  gridlength = minval(Radius) / refinements
  avgheight = sum(Particles(:,3)) / (1.0*numparticles)
  avgradius = sum(Radius) / (1.0*numparticles)
  Nx = ceiling(length/gridlength)+1
  Ny = ceiling(width/gridlength)+1
  !Nz = ceiling(2.0*avgradius/gridlength)
  Nz = ceiling(stopheight/gridlength)+1 ! goal is to get particles fully stacked up to 0.5 mm in height

  hx = length / (Nx-1.0)
  hy = width / (Ny-1.0)
  hz = stopheight / (Nz-1.0)

  print *, avgRadius, hx, hz


  ! 2D planar determination of whether particles reside within a sphere
  NodesTotal = (Nx)*(Nz)
  allocate(Xcoord(Nx))
  !allocate(Ycoord(NodesTotal))
  allocate(Zcoord(Nz))
  allocate(Inside(Nx,Nz))
  Inside(:,:) = 0
  counter = 1

  ! assign vector of coordinates
  do i = 1,Nx
     Xcoord(i) = (i-1.0)*hx
  end do

  do i = 1,Nz
     Zcoord(i) = (i-1.0)*hz 
  end do

!!$! 2D planar porosity calc
!!$NodesTotal = (Nx+1)*(Ny+1)
!!$allocate(Xcoord(NodesTotal))
!!$allocate(Ycoord(NodesTotal))
!!$allocate(Zcoord(NodesTotal))
!!$allocate(Points(NodesTotal))
!!$Points(:) = 0.0D0
!!$counter = 1
!!$do i = 1,Nx+1
!!$   do j = 1,Ny+1
!!$      Xcoord(counter) = (i-1)*gridlength
!!$      Ycoord(counter) = (j-1)*gridlength
!!$      Zcoord(counter) = avgRadius
!!$      counter = counter + 1 
!!$   end do
!!$end do

  ! Check to see which grid points lie within a sphere
  do i = 1,Nx
     do k = 1,Nz
        j = 1
        do while ((j <= numparticles) .AND. (Inside(i,k) < 0.5))
           distance = ((Xcoord(i)-Particles(j,1))**2.0 + (middleY-Particles(j,2))**2.0 + (Zcoord(k)-Particles(j,3))**2.0)**0.5
           if (distance <= Radius(j)) then
              Inside(i,k) = 1
           end if
           j = j+1
        end do
     end do
 !    print *, i
  end do


  ! Output data to Tecplot file
  write(100,*) 'Variables="X","Y","Inside"'
  write(100,*) 'Zone I=',Nx, 'J=',Nz, 'F=POINT'

  do i = 1,Nx
     do k = 1,Nz
        write(unit=100,fmt='(F13.6,2x,F13.6,2x,I6)') Xcoord(i), Zcoord(k), Inside(i,k)
     end do
  end do

  close(unit=100)


  ! Determine approximate density of particles
  powderdensity = sum(sum(Inside,1))*1.0 / (1.0*NodesTotal)

  print *, 'Density is ', powderdensity
  print *, 'Nx, Nz, numparticles', Nx, Nz, numparticles




end program SpheretoCircle

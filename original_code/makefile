all: 
	ifort LaserSinter.f90 -o LaserSinter -fast

compile2:
	ifort RDistributions.f90 LaserSinterParameterSiemens.f90 -o LaserSinterParameterSiemens -check all -warn all

compile3:
	ifort RDistributions.f90 PowderDeposition.f90 -o PowderDeposition -fast

compile4:
	ifort RDistributions.f90 PowderDepositionFast.f90 -o PowderDepositionFast -fast

compile5:
	ifort Porosity.f90 -o Porosity -fast

compile6:
	ifort LaserFDM.f90 -o LaserFDM -fast

compile7:
	ifort LaserFDSubstrate.f90 -o LaserFDSubstrate -fast

compile8:
	ifort RDistributions.f90 PowderDepositionLaserFlash.f90 -o PowderDepositionLaserFlash -fast

compile9:
	ifort MultiSLM.f90 -o MultiSLM -fast

compile10: 
	ifort SLMFast.f90 -o SLMSFast -fast

compile11:
	gfortran SLMFastOMP.f90 -o SLMFast0MP -fopenmp -ofast

compile12:
	ifort RDistributions.f90 PowderDepositionOMP.f90 -o PowderDepositionOMP -fast -openmp





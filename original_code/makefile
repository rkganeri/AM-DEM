all: 
	ifort LaserSinter.f90 -o LaserSinter -fast

compile2:
	ifort RDistributions.f90 LaserSinterParameterSiemens.f90 -o LaserSinterParameterSiemens -check all -warn all

compile3:
	ifort RDistributions.f90 PowderDepositionKhairallah.f90 -o PowderDepositionKhairallah -fast

compile4:
	ifort RDistributions.f90 PowderDepositionShankarFast.f90 -o PowderDepositionShankarFast -fast

compile5:
	ifort Porosity.f90 -o Porosity -fast

compile6:
	ifort LaserFDM.f90 -o LaserFDM -fast

compile7:
	ifort LaserFDSubstrate.f90 -o LaserFDSubstrate -fast

compile8:
	ifort RDistributions.f90 PowderDepositionLaserFlash.f90 -o PowderDepositionLaserFlash -fast

compile9:
	ifort MultiSLMShankar.f90 -o MultiSLMShankar -fast

compile10: 
	ifort SLMShankarFast.f90 -o SLMShankarFast -fast

compile11:
	gfortran SLMShankarFastMP2.f90 -o SLMShankarFast40MP2 -fopenmp -ofast

compile12:
	ifort RDistributions.f90 PowderDepositionShankarMP.f90 -o PowderDepositionShankarMP -fast -openmp





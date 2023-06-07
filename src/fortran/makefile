COMPILER=gfortran
mainflag=-fno-backslash -fbounds-check
mainlibs=-L`pwd` -lpgplot -L/usr/X11R6/lib -lX11 -lcfitsio -lnsl -lgcc -lm -lc -L/usr/include/libpng12 -lpng -lz -lFiltersv01 
ext_libs=L_FADO/64bits/OpenSUSE12.3
main_dir=/opt/BACKUP_S/FADOv01

FADOcom: *.f90
	#
	#
	# CREATE LINK for external libraries
	ln -s ${main_dir}/${ext_libs}/libcfitsio.a  
	ln -s ${main_dir}/${ext_libs}/libcpgplot.a  
	ln -s ${main_dir}/${ext_libs}/libpgplot.a

	#
	#
	# COMPILE individual routines
	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/M_FADO/ModDataType.f90 

	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/R_FADO/BASIC/*.f90 
	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/R_FADO/INTERPOLATION/*.f90 
	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/R_FADO/INTEGRATION/*.f90 

	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/F_FADO/ReadFilters.f90
	gfortran -fno-backslash -fbounds-check -O3 -c ${main_dir}/F_FADO/EvalFilters.f90

	#
	#
	# CREATE static library: libFiltersv01
	ar -x libcfitsio.a
	ar -x libpgplot.a
	ar rcs libFiltersv01.a *.mod *.o

	#
	#
	# Compile Main FILTERs code
	cp EvalFilters.CODE EvalFilters.CODE.f90
	gfortran -fno-backslash -fbounds-check EvalFilters.CODE.f90 -o EvalFilters.exe -O3 ${mainlibs}

	#	
	#
	# CLEANING
	rm *.o *.mod
	rm EvalFilters.CODE.f90

	#
	#
	# UNLINK
	unlink libcfitsio.a
	unlink libcpgplot.a
	unlink libpgplot.a

	# ENDING makefile of FILTERs

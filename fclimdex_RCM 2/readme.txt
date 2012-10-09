How to compile and run the code on Linux:

Assume the F90 compiler is ifort, and
the Netcdf lib is installed in /usr/local/netcdf/4.1.1-intel/.
We're using Fortran90 interface of netcdf now.

First load the netcdf modules:
   $ module load netcdf/4.1.1-intel

 Edit the Makefile if necessary, then
    $ make fclimdex.exe
 or
    $ ifort all.f90 -o fclimdex.exe -I/usr/local/netcdf/4.1.1-intel/include -L/usr/local/netcdf/4.1.1-intel/lib -lnetcdf

to run the code, edit the input.nml if necessary, then
    $ ./fclimdex.exe

If OpenMP is used, you can change processor core number # by setting
    $ ifort all.f90 -o fclimdex.exe -openmp -I/usr/local/netcdf/4.1.1-intel/include -L/usr/local/netcdf-intel/4.1.1/lib -lnetcdf
    $ export OMP_NUM_THREADS=#
    $ ./fclimedx.exe
by default, it uses all the cores on one machine.

To clean temporary files:
    $ make clean


-------------------------------------------------------------------
Note:
In case theres error, use -lnetcdff instead of -lnetcdf.

And make sure the lib environment is correct: LD_LIBRARY_PATH


---------------------------------------------------------------------
Features:
This version input and output netcdf files using Fortran90 interface

In this version, the input data names must contain "TX","TN" and "PR",
Otherwise please modify the subroutine "read_file" and re-compile the code.
You can separate the names (in "infilename_temp.txt") by blanks, tabulars or ","s.
you can add any comment by using "!" first, and information will be ignored.
The names can be in any order, and you can write names even without the ".nc".

Data (Tmax,Tmin,prcp) must have three dimensions in order of lon,lat,time.
And variables must be saved in order of lon, lat, time, data.
This also applies to the output of the index.
The code will be smarter later.

The output index can be put in one folder, or in different sub folders, controlled in inp.nml.

The code is not fast yet, will improve performance soon.

By Hongang Yang @ 2011.9.5:  hongang.yang@unsw.edu.au
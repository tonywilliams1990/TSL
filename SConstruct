import os.path
import os
import glob


# compiler
c_comp = 'g++'
f_comp = 'gfortran'

lapack = blas = slepc = petsc = 1

# library names
TSL_lib = 'TSL'
blas_lib = 'blas'
lapack_lib = 'lapack'
# the PETSC_DIR, SLEPC_DIR and PETSC_ARCH variables must be set first!
slepc_lib = 'slepc'
petsc_lib = 'petsc'
petsc_inc = ''
mpi_lib = 'mpi'
mpi_inc = ''

# default output format for messaging
red = "\033[1;31m"
yellow = "\033[1;33m"
green = "\033[1;32m"
blue = "\033[1;34m"
off = "\033[0m"

def message( col, text ):
    print col + " * " + text + off

# Output message to user
message( blue, " -----  Building -----")

# set the build dir
topdir     = os.getcwd()
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
libs_str   = TSL_lib + ' '
preproc    = ''
#opts = ' -O2 -std=c++14 -Wall -Wextra '
opts = ' -O2 '
link_flags = ' '

# Set the rpath for the linker to find petsc/slepc-3
rpath = []

# Get variables
petsc_dir  = os.environ.get('PETSC_DIR','default')
petsc_arch = os.environ.get('PETSC_ARCH','default')
slepc_dir  = os.environ.get('SLEPC_DIR','default')

if int(lapack):
    print( " *" )
    message( green, "LAPACK support is enabled.")
    message( blue, " I'm assuming a GNU compiler")
    libs_str += blas_lib + ' ' + lapack_lib + ' gfortran '

if int(petsc * slepc):
    print( " *" )
    message( green, "PETSc & SLEPc library support is enabled. ")
    message( blue, " $PETSC_DIR/x86_64-linux-gnu-complex/lib.")
    message( blue, " and $SLEPC_DIR/x86_64-linux-gnu-complex/lib")

    petsc_lib_dir = petsc_dir + "/" + petsc_arch + "/lib"
    petsc_inc     = petsc_dir + "/include " + petsc_dir + "/" + petsc_arch + "/include "
    rpath.append( petsc_lib_dir )

    slepc_lib_dir = slepc_dir + "/" + petsc_arch + "/lib"
    slepc_inc     = slepc_dir + "/include " + slepc_dir + "/" + petsc_arch + "/include "
    rpath.append( slepc_lib_dir )

    preproc += ' -DPETSC_Z '
    preproc += ' -DSLEPC '
    preproc += ' -DINC_MPI '

libs_str   += slepc_lib + ' ' + petsc_lib + ' ' + mpi_lib + ' '
libdir_str += slepc_lib_dir + ' ' + petsc_lib_dir + ' '
incdir_str += slepc_inc + ' ' + petsc_inc + ' ' + mpi_inc + ' '

# source is all cpp files, so let's glob them
src = glob.glob('src/*.cpp')
blas_src = glob.glob('BLAS/SRC/*.f')
lapack_src = glob.glob('LAPACK/SRC/*.f') + glob.glob('LAPACK/SRC/*.F') + glob.glob('LAPACK/INSTALL/*.f')

# Split the strings into list elements for lib/inc directoris
incdir = incdir_str.split()
libdir = libdir_str.split()

# Initialise the environment
env = Environment( FORTRAN = f_comp, CXX = c_comp, CPPPATH = incdir, CCFLAGS = opts + preproc, LINKFLAGS = link_flags, LIBPATH = libdir)

# Now check the environment is complete
conf = Configure( env )

'''
# Check for external libraries
if not conf.CheckLib( blas_lib ):
    # check blas
    message( red, "No blas library! It is required to enable LAPACK")
    blas = 0
if not conf.CheckLib( lapack_lib ):
    # check lapack
    message( red, "No lapack library!")
    lapack = 0
if ( blas * lapack  == 0 ):
    message( red, "Check your library path ... without the above")
    message( red, "listed libraries, I can't compile with LAPACK support.")
    Exit(1)
else:
    message( green, "Found BLAS & LAPACK support.")
    message( green, "LAPACK solvers will be used in preference to the native ones.")

if not conf.CheckLib( petsc_lib ):
    message( red, "No libpetsc!")
else:
	message( green, "Found PETSC.")

if not conf.CheckLib( slepc_lib ):
    message( red, "No libslepc!")
else:
	message( green, "Found SLEPC.")

slepc = blas = lapack = mpi = 1
if not conf.CheckLib('mpi'):
    message( red, "No libmpi!")
    slepc = 0
if not conf.CheckLib('slepc'):
    message( red, "No libslepc!")
    slepc = 0
if not conf.CheckLib('petsc'):
    message( red, "No libpetsc!")
    petsc = 0
if not conf.CheckLib('lapack'):
    message( red, "No liblapack!")
    lapack = 0
if not conf.CheckLib('blas'):
    message( red, "No libblas!")
    blas = 0
if ( slepc * petsc * lapack * blas == 0 ):
    message( red, "SLEPC support has failed.")
    Exit(1)
else:
    message( green, "Found SLEPc and PETSc, including support for sparse matrix eigensolvers.")
'''

env = conf.Finish()

#libs_str   += blas_lib + ' ' + lapack_lib + ' gfortran'
libs = libs_str.split()

# Build the libraries in ./lib
env.StaticLibrary('lib/' + blas_lib, blas_src)
env.StaticLibrary('lib/' + lapack_lib, [blas_src, lapack_src]) # LAPACK requires BLAS
env.StaticLibrary('lib/' + TSL_lib, [blas_src, lapack_src, src] ) # TSL requires LAPACK & BLAS

SConscript('Examples/SConscript', exports='env opts preproc link_flags rpath incdir libdir topdir libs' )

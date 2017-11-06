import os.path
import glob

# compiler
c_comp = 'g++'
f_comp = 'gfortran'

# library names
TSL_lib = 'TSL'
blas_lib = 'blas'
lapack_lib = 'lapack'

# source is all cpp files, so let's glob them
src = glob.glob('src/*.cpp')
blas_src = glob.glob('BLAS/SRC/*.f')

# set the build dir
topdir = os.getcwd()
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
#opts = ' -O2 -std=c++14 -Wall -Wextra '
opts = ' -O2 -Wall '

# Split the strings into list elements for lib/inc directoris
incdir = incdir_str.split()
libdir = libdir_str.split()

# Initialise the environment
env = Environment( FORTRAN = f_comp, CXX = c_comp, CPPPATH = incdir, CCFLAGS = opts, LIBPATH = libdir)

# Now check the environment is complete
conf = Configure( env )
env = conf.Finish()

libs_str   = TSL_lib + ' ' + blas_lib
libs = libs_str.split()

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

# Build the libraries in ./lib
env.StaticLibrary('lib/blas', blas_src)
env.StaticLibrary('lib/TSL', src )

SConscript('Examples/SConscript', exports='env opts incdir libdir topdir libs' )

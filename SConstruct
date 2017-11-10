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
lapack_src = glob.glob('LAPACK/SRC/*.f') + glob.glob('LAPACK/SRC/*.F') + glob.glob('LAPACK/INSTALL/*.f')

# set the build dir
topdir = os.getcwd()
incdir_str = topdir + '/include '
libdir_str = topdir + '/lib '
#opts = ' -O2 -std=c++14 -Wall -Wextra '
opts = ' -O2 -Wall '
link_flags = ''

# Split the strings into list elements for lib/inc directoris
incdir = incdir_str.split()
libdir = libdir_str.split()

# Initialise the environment
env = Environment( FORTRAN = f_comp, CXX = c_comp, CPPPATH = incdir, CCFLAGS = opts, LINKFLAGS = link_flags, LIBPATH = libdir)

# Now check the environment is complete
conf = Configure( env )
env = conf.Finish()

libs_str   = blas_lib + ' ' + lapack_lib  + ' ' + TSL_lib + ' gfortran'
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
env.StaticLibrary('lib/' + blas_lib, blas_src)
env.StaticLibrary('lib/' + lapack_lib, [blas_src, lapack_src]) # LAPACK requires BLAS
env.StaticLibrary('lib/' + TSL_lib, [blas_src, lapack_src, src] ) # TSL requires LAPACK & BLAS



SConscript('Examples/SConscript', exports='env opts link_flags incdir libdir topdir libs' )

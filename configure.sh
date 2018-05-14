#!/bin/bash
home=/home/tony/Desktop/TSL
#echo $home
tar -xvzf petsc-3.9.1.tar.gz
tar -xvzf slepc-3.9.1.tar.gz

export PETSC_DIR=$home/petsc-3.9.1
export SLEPC_DIR=$home/slepc-3.9.1
export PETSC_ARCH=x86_64-linux-gnu-complex

cd $home/petsc-3.9.1
./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --download-fblaslapack --download-superlu_dist --download-mpich --download-metis --download-parmetis --download-cmake --download-scalapack --download-mumps --with-scalar-type=complex
make PETSC_DIR=$home/petsc-3.9.1 PETSC_ARCH=x86_64-linux-gnu-complex all
cd $home/slepc-3.9.1
./configure
make SLEPC_DIR=$home/slepc-3.9.1 PETSC_DIR=$home/petsc-3.9.1 PETSC_ARCH=x86_64-linux-gnu-complex
cd $home
export LD_LIBRARY_PATH=$home/lib:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib

# put new variables in bash_profile
echo "export PETSC_DIR=$home/petsc-3.9.1" >> ~/.bash_profile
echo "export SLEPC_DIR=$home/slepc-3.9.1" >> ~/.bash_profile
echo "export PETSC_ARCH=x86_64-linux-gnu-complex" >> ~/.bash_profile
echo "export LD_LIBRARY_PATH=$home/lib:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib" >> ~/.bash_profile

# update shell
source ~/.bash_profile

# uncomment stuff below if bashrc profile is different

# put new variables in bashrc
#echo "export PETSC_DIR=$home/petsc-3.8.2" >> ~/.bashrc
#echo "export SLEPC_DIR=$home/slepc-3.8.2" >> ~/.bashrc
#echo "export PETSC_ARCH=x86_64-linux-gnu-complex" >> ~/.bashrc
#echo "export LD_LIBRARY_PATH=$home/lib:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib" >> ~/.bashrc

# update shell
#source ~/.bashrc

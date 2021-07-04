#!/bin/zsh

rm -rf *.o
rm -rf *.mexa*

/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_mex.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_normalize.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_eval.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_coloc.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_coloc_deriv.cpp
/media/sdb1/Matlab_Linux/bin/mex -g -c -largeArrayDims-v GCC=/usr/bin/gcc-6 COMPFLAGS="$COMPFLAGS /openmp /O2 /favor:EM64T -DNTHREADS=4" bbs_bending_ur.cpp


/media/sdb1/Matlab_Linux/bin/mex -v GCC='/usr/bin/gcc-6' bbs_normalize.o bbs.o bbs_mex.o
/media/sdb1/Matlab_Linux/bin/mex -v GCC='/usr/bin/gcc-6' bbs_eval.o bbs.o bbs_mex.o
/media/sdb1/Matlab_Linux/bin/mex -v GCC='/usr/bin/gcc-6' bbs_coloc.o bbs.o bbs_mex.o
/media/sdb1/Matlab_Linux/bin/mex -v GCC='/usr/bin/gcc-6' bbs_coloc_deriv.o bbs.o bbs_mex.o
/media/sdb1/Matlab_Linux/bin/mex -v GCC='/usr/bin/gcc-6' bbs_bending_ur.o bbs.o bbs_mex.o

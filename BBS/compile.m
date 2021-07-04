function compile

% BBS is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
% 
% BBS is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

debug = 0;
nthreads = 4;
cpp_flags = '/openmp /O2 /favor:EM64T';
mex_opts = '-largeArrayDims';

delete *.o
delete *.obj
delete *.mex*
delete *.pdb
delete *.ilk

files = {'bbs.cpp', 'bbs_mex.cpp', 'bbs_normalize.cpp', 'bbs_eval.cpp', 'bbs_coloc.cpp', 'bbs_coloc_deriv.cpp', 'bbs_bending_ur.cpp'};

%mex_opts = [mex_opts '-v GCC=/usr/bin/gcc-6'];

for ii = 1:numel(files)
    fprintf('Compiling "%s"\n', files{ii});
    fprintf('mex -g -c %s COMPFLAGS="$COMPFLAGS %s -DNTHREADS=%d" %s\n', mex_opts, cpp_flags, nthreads, files{ii});
    eval(sprintf('mex -g -c %s COMPFLAGS="$COMPFLAGS %s -DNTHREADS=%d" %s', mex_opts, cpp_flags, nthreads, files{ii}))
end

mex -v GCC='/usr/bin/gcc-6' bbs_normalize.o bbs.o bbs_mex.o
mex -v GCC='/usr/bin/gcc-6' bbs_eval.o bbs.o bbs_mex.o
mex -v GCC='/usr/bin/gcc-6' bbs_coloc.o bbs.o bbs_mex.o
mex -v GCC='/usr/bin/gcc-6' bbs_coloc_deriv.o bbs.o bbs_mex.o
mex -v GCC='/usr/bin/gcc-6' bbs_bending_ur.o bbs.o bbs_mex.o

delete *.o

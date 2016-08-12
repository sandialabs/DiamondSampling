This directory contains both Matlab and C/C++ code and is more 
efficient than the pure Matlab code.  The C/C++ code uses primitives
from the CSparse library, which is included in the CSparse directory.

1. Open MATLAB
2. Type "mex -setup" if mex is not already setup on the machine
3. Navigate to this directory
4. Type "setup" to build CSparse functions
  -- compiling CSparse will take several seconds 
5. Type "compile" to build mads code
6. Type "script_check_exact" and "script_check_sampling" to test the code
  -- script_check_exact should return 3 successes
  -- script_check_sampling should generate two loglog plots showing
     error decreasing as the number of samples increases
7. (Optional) Type "clean" to delete mex-compiled files


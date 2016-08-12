% clean MADS object files and executables
delete *.o
eval(['delete *.' mexext])

% clean CSparse object files and executables
cd CSparse/MATLAB/CSparse
delete *.o
eval(['delete *.' mexext])
cd ../../..
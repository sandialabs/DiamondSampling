% start building mex command, use 64-bit integers and optimizations
mexcmd = 'mex -largeArrayDims -O ';

% define path to CSparse library (assume it has been mexed)
CS_path = 'CSparse/';
  
% add include directives for cs.h and cs_mex.h
mexcmd = [mexcmd '-I' CS_path 'MATLAB/CSparse -I' CS_path '/Include'];

% compile sampling functions
eval([mexcmd ' -c sampling.c'])

try
    % ispc does not appear in MATLAB 5.3
    pc = ispc ;
catch
    % if ispc fails, assume we are on a Windows PC if it's not unix
    pc = ~isunix ;
end

if (pc)
    obj = '.obj' ;
else
    obj = '.o' ;
end

% add CSparse object files to link against
cs = { 'cs_mex', 'cs_util', 'cs_malloc', 'cs_multiply', 'cs_scatter', ...
    'cs_transpose', 'cs_cumsum', 'cs_fkeep', 'cs_dropzeros', 'cs_dupl', ...
    'cs_compress', 'cs_print', 'cs_norm'} ;
CS = [];
for i = 1:length (cs)
    CS = [CS ' ' CS_path 'MATLAB/CSparse/' cs{i} obj] ;                                              
end    
    
% compile mex functions
eval([mexcmd ' exact_topk_mex.cpp -lmwblas' CS])
eval([mexcmd ' sample_ab_mex.c sampling' obj ' ' CS])


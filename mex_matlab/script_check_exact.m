A = sprandn(10,9,.7);
B = sprandn(9,11,.4);
k = 5;

% compute A*B in full and sort for top k entries
C = A*B;
[ii,jj,vv] = find(C);
[~,idx] = sort(abs(vv),'descend');
Ctop = sparse(ii(idx(1:k)),jj(idx(1:k)),vv(idx(1:k)),size(C,1),size(C,2));

% compute top k entries of A*B using priority queue
D = exact_topk_mex(A,B,k,'sparse');

assert(nnz(D-Ctop)==0);
fprintf('success for sparse case\n');

% get strict upper triangle of A*B and sort for top k entries
Cupper = triu(C,1);
[ii,jj,vv] = find(Cupper);
[~,idx] = sort(abs(vv),'descend');
Ctopupper = sparse(ii(idx(1:k)),jj(idx(1:k)),vv(idx(1:k)),size(C,1),size(C,2));

% compute top k entries of strict upper triangle of A*B using priority queue
Dupper = exact_topk_mex(A,B,k,'sparse','upper');

assert(nnz(Dupper-Ctopupper)==0);
fprintf('success for sparse upper case\n');




A = randn(10,3);
B = randn(3,11);
k = 5;

% compute A*B in full and sort for top k entries
C = A*B;
[ii,jj,vv] = find(C);
[~,idx] = sort(abs(vv),'descend');
Ctop = sparse(ii(idx(1:k)),jj(idx(1:k)),vv(idx(1:k)),size(C,1),size(C,2));

% compute top k entries of A*B using priority queue
D = exact_topk_mex(A,B,k,'dense');

assert(nnz(D-Ctop)==0);
fprintf('success for dense case\n');
function [C,info,adata,bdata] = atb(A,B,nsamples,varargin)
%ATB Sampling to find largest entries of A'*B.
%
%   C = ATB(A,B,N) computes the top 10 largest magnitude entries of the
%   matrix product A'*B via a statistical sampling approach. The result C
%   has only the top 10 entries, as estimated by the procedure.
%
%   C = ATB(A,B,N,param,value,...) accepts parameter-value pairs as
%   follows...
%
%   'ntop'        - Number of top entries to compute. Default: 10.
%   'nbudget'     - Budget of dot products for computing top entries.
%                   Default: 10*ntop.
%   'type'        - Type of sampling, either 'wedge' or 'diamond'. 
%                   Default: 'diamond'.
%   'nodiag'      - Do not consider diagonal entries of A'*B in computation
%                   of the largest entries. This only makes sense in the
%                   case that B == A. Default: false. 
%   'binary'      - Treat the inputs as binary. Default: autodetected.
%   'nonnegative' - Treat the inputs as nonnegative. Default: autodetected.
%   'dense'       - Treat the inputs as dense. Default: autodetected.
%   'ata'         - B == A. Default: autodetected.
%   'asqr'        - B == A and A == A'. Default: autodetected.
%   'skipchecks'  - Skip all autodetect and error checks. If true,
%                   'binary', 'nonnegative', 'ata', and 'asqr' are set to
%                   be false unless explicitly specified. Default: false. 
%   'adata'       - Reuse certain data about A matrix from previous run.
%                   This is useful if the same A matrix is used repeatedly
%                   with different B matrices. This data is returned as the
%                   third argument. Default: [].
%   'bdata'       - Reuse certain data about B matrix from previous run.
%                   This is useful if the same B matrix is used repeatedly
%                   with different A matrices. This data is returned as the
%                   fourth argument. Default: [].
%
%   

%% Process inputs
starttime = tic;

params = inputParser;
params.addParameter('binary',[]);
params.addParameter('nonnegative',[]);
params.addParameter('dense',[]);
params.addParameter('nodiag',false);
params.addParameter('ata',[]);
params.addParameter('asqr',[]);
params.addParameter('ntop',10);
params.addParameter('nbudget',[]);
params.addParameter('type','diamond',@(x) ismember(x,{'diamond','wedge'}));
params.addParameter('skipchecks',false);
params.addParameter('adata',[]);
params.addParameter('bdata',[]);
params.parse(varargin{:});

skipchecks = params.Results.skipchecks;
binary = params.Results.binary;
nonnegative = params.Results.nonnegative;
nodiag = params.Results.nodiag;
ata = params.Results.ata;
asqr = params.Results.asqr;
ntop = params.Results.ntop;
nbudget = params.Results.nbudget;
sampletype = params.Results.type;
dense = params.Results.dense;
adata = params.Results.adata;
bdata = params.Results.bdata;

wedges = isequal(sampletype,'wedge');

% Budget
if isempty(nbudget)
    nbudget = 10*ntop;
end

% Binary
if skipchecks
    if isempty(binary)
        binary = false;
    end
else
    binary_actual = all(nonzeros(A)==1) & all(nonzeros(B)==1);
    if isempty(binary)
        binary = binary_actual;
    elseif binary && ~binary_actual
        warning('Converting non-binary matrices to binary');
        A = spones(A);
        B = spones(B);
    elseif ~binary && binary_actual
        fprintf('*** WARNING: Treating binary matrix as weighted ***\n');
    end
end

% Nonnegative
if skipchecks
    if isempty(nonnegative)
        nonnegative = false;
    end
else
    nonnegative_actual = all(nonzeros(A) > 0) & all(nonzeros(B) > 0);
    if isempty(nonnegative)
        nonnegative = nonnegative_actual;
    elseif nonnegative && ~nonnegative_actual
        error('Matrices have negative entries but nonnegative = true');
    elseif ~nonnegative && nonnegative_actual
        warning('Since nonnegative = false, treating nonnegative matrices as general weighted');
    end
end

if skipchecks
    if isempty(ata)
        ata = false;
    end
else    
    ata_actual = isequal(A,B);
    if isempty(ata)
        ata = ata_actual;
    elseif ata && ~ata_actual
        error('Do not treat asymmetric result as symmetric when A ~= B');
    elseif ~ata && ata_actual
        warning('Since ata = false, ignoring the fact that B=A');
    end
end

if skipchecks
    if isempty(asqr)
        asqr = false;
    end
else
    asqr_actual = ata && isequal(A,A');
    if isempty(asqr)
        asqr = asqr_actual;
    elseif asqr && ~asqr_actual
        error('A is not symmetric but asqr = true.');
    elseif ~asqr && asqr_actual
        warning('Since asqr = false, ignoring the fact that A=A');
    end
end

if skipchecks
    if isempty(dense)
        dense = false;
    end    
else
    density_A = nnz(A)/numel(A);
    density_B = nnz(B)/numel(B);
    dense_actual = (density_A > 0.25) || (density_B > 0.25);
    if isempty(dense)
        dense = dense_actual;
    elseif strcmp(dense,'forcefalse')
        dense = false;
    elseif dense && ~dense_actual
        warning('Treating problem as dense even though it is sparse. Set option ''dense'' to ''forcefalse'' to supress this warning.');
    elseif ~dense && dense_actual
        warning('Treating problem as sparse even though it is dense');
    end
end

if nodiag && ~ata
    error('Do not omit diagonal when A~= B');
end

if binary || nonnegative
    absA = A;
    absB = B;
else
    if isempty(adata)
        absA = abs(A);
    else
        absA = adata{1};
    end
    if isempty(bdata)     
        absB = abs(B);
    else
        absB = bdata{1};
    end
end

if isempty(adata)
    absAt = absA';
else
    absAt = adata{2};
end
if isempty(bdata)
    absBt = absB';
else
    absBt = bdata{2};
end

nrejects_nodiag = 0;
nrejects_noselfedges = 0;
ndiamonds = 0;

n1 = size(A,2);
n2 = size(B,2);
m = size(A,1);
if m ~= size(B,1)
    error('B must have the same number of rows as A');
end

time_preprocess = toc(starttime);

%% Sampling
starttime = tic;

Arw = absA*ones(n1,1);
if ata
    Brw = Arw;
else
    Brw = absB*ones(n2,1);
end

if wedges 
    if binary    
        ctrweight = Arw .* Brw;
        K = random_sample(ctrweight,nsamples);
        I = sample_nghbrs(K,absAt,Arw);
        J = sample_nghbrs(K,absBt,Brw);
    else
        ctrweight = Arw .* Brw;
        K = random_sample(ctrweight,nsamples);
        I = sample_weighted_nghbrs(K,absAt);
        J = sample_weighted_nghbrs(K,absBt);
    end
    nghbr1 = I;
    nghbr2 = J;
    sumw = sum(ctrweight);
    if nonnegative == false
        idxa = tt_sub2ind(size(A),[K,I]);
        sgna = full(sign(A(idxa)));
        idxb = tt_sub2ind(size(B),[K,J]);
        sgnb = full(sign(B(idxb)));
        vals = sgna .* sgnb;
    else
        vals = ones(size(I));
    end
else
       
    Acw = absAt*ones(m,1);
    
    if dense
        if isempty(adata)
            W1 = bsxfun(@times, Acw, absAt);
        else
            W1 = adata{3};
        end
        W = bsxfun(@times, W1, Brw');
        sumw = sum(W(:));
        edgesample = random_sample(W(:),nsamples);
        subs = tt_ind2sub(size(W),edgesample);
        I = subs(:,1);
        K = subs(:,2);
    else
        Dbr = spdiags(Brw,0,m,m);
        Dac = spdiags(Acw,0,n1,n1);
        W = Dac * absAt * Dbr;
        sumw = sum(nonzeros(W));
        [WI,WK,WV] = find(W);
        edgesample = random_sample(WV,nsamples);
        I = WI(edgesample);
        K = WK(edgesample);
    end
    
    if binary   
        J = sample_nghbrs(K,absBt,Brw);
        L = sample_nghbrs(I,absA,Acw);        
    else
        J = sample_weighted_nghbrs(K,absBt);
        L = sample_weighted_nghbrs(I,absA);
    end
    
    
    tmp = tt_sub2ind(size(B),[L,J]); % This is a Tensor Toolbox function
    vals = full(B(tmp));
    tf = ~(vals == 0);
    ndiamonds = sum(tf);
    I = I(tf);
    J = J(tf);
    K = K(tf);
    L = L(tf);
    vals = vals(tf);

    if nonnegative == false
        idxaki = tt_sub2ind(size(A),[K I]);
        sgnaki = full(sign(A(idxaki)));
        idxbkj = tt_sub2ind(size(B),[K J]);
        sgnbkj = full(sign(B(idxbkj)));
        idxali = tt_sub2ind(size(A),[L I]);
        sgnali = full(sign(A(idxali)));
        vals = sgnaki .* sgnbkj .* sgnali .* vals;
    end
    
    if asqr   
        nghbr1 = [I; K];
        nghbr2 = [J; L];
        vals = [vals; vals];
    else
        nghbr1 = I;
        nghbr2 = J;
        vals = vals;
    end

end

time_sample = toc(starttime);

%% Postprocess samples
starttime = tic;

% Remove any diagonal values
if nodiag
    tf = (nghbr1 == nghbr2);
    nrejects_nodiag = sum(tf);
    nghbr1 = nghbr1(~tf);
    nghbr2 = nghbr2(~tf);
    vals = vals(~tf);
end

% Make everything upper triangular (may be repeats - that's okay)
if ata
    tf = (nghbr2 < nghbr1);
    tmp = nghbr2;
    nghbr2(tf) = nghbr1(tf);
    nghbr1(tf) = tmp(tf);
    
    tf = (nghbr1 == nghbr2);
    vals(tf) = vals(tf)*2;
    
    vals = vals ./ 2;
end

if asqr && ~wedges,
    vals = vals ./ 2;
end

% Construct X
X = sparse(nghbr1,nghbr2,vals,n1,n2);

% Adjust sizes, if needed
nbudget = min(nbudget,nnz(X));
ntop = min(ntop,nnz(X));

% Construct C, using the biggest entries of X
[ii,jj,vv] = find(X);
[~,idx] = sort(abs(vv),'descend');
idx = idx(1:nbudget);
ii = ii(idx);
jj = jj(idx);

% Compute the dot products for all the entries in the budget
dp = zeros(nbudget,1);
for k = 1:nbudget
    dp(k) = A(:,ii(k))'*B(:,jj(k));
end

% Use the very top entries to construct C
[~,idx] = sort(abs(dp),'descend');
idx = idx(1:ntop);
C = sparse(ii(idx), jj(idx), dp(idx), n1, n2);

time_postprocess = toc(starttime);

%% Save data
if nargout >= 2
    info.X = X;
    info.sumw = sumw;
    info.binary = binary;
    info.nonnegative = nonnegative;
    info.nodiag = nodiag;
    info.ata = ata;
    info.asqr = asqr;
    info.sampletype = sampletype;
    info.ntop = ntop;
    info.nbudget = nbudget;
    info.nrejects_noselfedges = nrejects_noselfedges;
    info.nrejects_nodiag = nrejects_nodiag;
    info.ndiamonds = ndiamonds;
    info.time_preprocess = time_preprocess;
    info.time_sample = time_sample;
    info.time_postprocess = time_postprocess;
end

if nargout >= 3
    if ~exist('W1','var')
        W1 = [];
    end
    adata = {absA, absAt, W1};
    bdata = {absB, absBt};
end



function samples = sample_nghbrs(root,A,deg)
%SAMPLE_COLUMNS Unweighted samples of columns of A

nsamples = size(root,1);

% Set-up: the index computation exploits the fact that the entries of A are
% sorted by column index. 
nbr_start = [1; 1+cumsum(deg(1:end-1))];
[nbr,~,~] = find(A);

% Sample 1st neighbor
idx1 = floor(deg(root) .* rand(nsamples,1));
samples = nbr(nbr_start(root) + idx1);

function samples = sample_weighted_nghbrs(root,A)
%SAMPLE_WEIGHTED NBGHRS Weighted sampling of columns of A

n = size(A,2);
[rootsrt,rootidx] = sort(root,'ascend');
rootcnt = accumarray(rootsrt,1,[n 1]);
idx = [1; cumsum(rootcnt)+1];
samples = zeros(size(root));
for k = 1:n
    if rootcnt(k) == 0
        continue;
    end
    [i,~,v] = find(A(:,k));
    tmpsample = random_sample(v,rootcnt(k));
    nghbrsrt1(idx(k):idx(k+1)-1) = i(tmpsample);
end
samples(rootidx,1) = nghbrsrt1;
function Ctop = topk(X,A,B,k,nbudget,varargin)
% recompute top k entries of C=A'*B based on X estimates,
% using re-computation of nbudget sparse dot products

    % Check for AtA flag
    params = inputParser;
    params.addParameter('AtA',0);
    params.parse(varargin{:});
    if params.Results.AtA ~= 0,
        X = X + X'; % symmetrize
        X = triu(X,1); % ignore diagonal
    end
   
    % Adjust sizes, if needed
    nbudget = min(nbudget,nnz(X));
    k = min(k,nnz(X));

    
    % Construct Ctop, using the biggest entries of X
    [ii,jj,vv] = find(X);
    if nbudget < nnz(X),    
        [~,idx] = sort(abs(vv),'descend');
        idx = idx(1:nbudget);
        ii = ii(idx);
        jj = jj(idx);
    end

    % Compute the dot products for all the entries in the budget
    dp = zeros(nbudget,1);
    for kk = 1:nbudget
        dp(kk) = A(:,ii(kk))'*B(:,jj(kk));
    end

    % Use the very top entries to construct Ctop
    [~,idx] = sort(abs(dp),'descend');
    idx = idx(1:k);
    Ctop = sparse(ii(idx), jj(idx), dp(idx), size(A,2), size(B,2));
    
end
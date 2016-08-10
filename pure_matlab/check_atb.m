function [relerr,A,B,X] = check_atb(K,M,N,p,nsamples,producttype,matrixtype,type)

switch producttype
    case {'asqr'}
        A = sprandn(K,K,p);
        A = (A+A')/2;
        B = A;
        P = triu(ones(K),0); % Pattern
    case {'ata'}
        A = sprandn(K,M,p);
        B = A;
        P = triu(ones(M),0); % Pattern
    case {'atb'}
        A = sprandn(K,M,p);
        B = sprandn(K,N,p);
        P = ones(M,N); % Pattern
end

switch matrixtype
    case {'binary'}
        A = spones(A);
        B = spones(B);
    case {'nonnegative'}
        A = abs(A);
        B = abs(B);
    case {'pmones'}
        A = sign(A).*spones(A);
        B = sign(B).*spones(B);            
end

switch type
    case{'wedge'}
        E = A'*B;
        Eplus = abs(A)'*abs(B);
    case{'diamond'}
        E = (A'*B).^2;
        Eplus = abs(A)'*abs(B);
end

P = (Eplus~=0) & P; % Pattern of valid nonzeros in X
Evals = full(E(P==1));
Epvals = full(Eplus(P==1));

[~,info] = atb(A,B,nsamples,'type',type,'dense','forcefalse');
X = info.X;
Xbadvals = X(P==0);
if ~isempty(find(Xbadvals,1))
    error('Invalid count in X!');
end

ns = nsamples;
Xvals = full(X(P==1));
Yvals = (info.sumw/ns) * Xvals;
relerr = norm(Yvals-Evals)./norm(Evals);
fprintf('Relative error = %g\n', relerr)

end

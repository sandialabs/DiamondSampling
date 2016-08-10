function s = random_sample(p, nsamples, method)
%RANDOM_SAMPLE creates a random sample proportional to the given counts.
%
%   S = RANDOM_SAMPLE(P) choose N = round(sum(P)) samples (with
%   replacement) from {1,...,length(P)} proportional to the values in P.
%   So, if P = [2 1 1], then we might expect S (sorted) to be [ 1 1 2 3 ].
%   However, we also allow for P to be non-integral.
%
%   S = RANDOM_SAMPLE(P,N) creates N random samples with replacement
%   proportional to the values in P.
%
%   S = RANDOM_SAMPLE(P,N,'alias') used the alias method to compute the
%   random samples. This requires O(length(P)) setup time, but then only
%   O(N) time for the samples.
%
%   The alias method is adapted from the following blog post:
%   http://possiblywrong.wordpress.com/2012/02/05/
%   the-alias-method-and-double-precision/ 
%   


if ~exist('nsamples','var')
    nsamples = round(sum(p));
end

if ~exist('method','var')
    method = 'standard';
end

if strcmpi(method,'standard')

%p = p .* nsamples / sum(p);
cumdist = [0; cumsum(p)];
bins = cumdist / cumdist(end);

testval = abs(bins(end) - 1);
if  testval > eps
    warning('Last entry of bins is not exactly 1. Diff = %e.', testval);
end

[~, s] = histc(rand(nsamples,1),bins);

elseif strcmpi(method(1:5),'alias')
    
    % Setup
    n = numel(p);
    prob = nan(1,n);
    alias = nan(1,n);
    p = p * (n / sum(p));
    small = nan(1,n);
    large = nan(1,n);
    ns = 0;
    nl = 0;
    for j = 1:n
        if p(j) > 1
            nl = nl + 1;
            large(nl) = j;
        else
            ns = ns + 1;
            small(ns) = j;
        end
    end
    while ns ~= 0 && nl ~= 0
        j = small(ns);
        ns = ns - 1;
        k = large(nl);
        nl = nl - 1;
        prob(j) = p(j);
        alias(j) = k;
        p(k) = p(k) + p(j) - 1;
        if p(k) > 1
            nl = nl + 1;
            large(nl) = k;
        else
            ns = ns + 1;
            small(ns) = k;
        end
    end
    while ns ~= 0
        prob(small(ns)) = 1;
        ns = ns - 1;
    end
    while nl ~= 0
        prob(large(nl)) = 1;
        nl = nl - 1;
    end
    
    % Sample
    if strcmpi(method(end-2:end),'vec')
        u = n * rand(1,nsamples);
        j = floor(u);
        s = alias(j + 1);
        idx = (u - j <= prob(j + 1));
        s(idx) = j(idx) + 1;
    else
        s = zeros(1,nsamples);
        u = n * rand(1,nsamples);
        for i = 1:nsamples
            j = floor(u(i));
            if u(i) - j <= prob(j + 1)
                s(i) = j + 1;
            else
                s(i) = alias(j + 1);
            end
        end
    end
    
else
    error('Invalid method');
end
    
    
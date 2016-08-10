%% Test correctness of A'*B code

% For 'asqr', A = K x K
% For 'ata', A = K x M
% For 'atb', A = K x M, B = K x N

% K = 10; M = 5; N = 5; p = .5;
% nsamples = 1e2;
% producttype = 'atb';
% matrixtype = 'binary';
% type = 'diamond';
% [relerr,A,B,X] = check_atb(K,M,N,p,nsamples,producttype,matrixtype,type);
% relerr
% full(A'*B)
% full(round(X))

%% Make some plots

producttype = 'asqr';
matrixtype = 'binary';
K = 10;
p = 0.3;

tvals = {'wedge','diamond'};
nsvals = 10.^(2:7);

clear relerr
for t = 1:2,
    for ns = 1:length(nsvals)
        fprintf('%s, %d samples\n',tvals{t},nsvals(ns));
        relerr(t,ns) = check_atb(K,[],[],p,nsvals(ns),producttype,matrixtype,tvals{t});        
    end    
end
  
%%
figure(1);
semilogx(nsvals,relerr(1,:),nsvals,relerr(2,:))
legend('Wedges','Diamond');
xlabel('nsamples');
ylabel('relative error');
tstr = sprintf('Binary A^2 for %d x %d matrix with %d%% nonzeros', K,K,round(p*100));
title(tstr);

%%
producttype = 'ata';
matrixtype = 'binary';
K = 10;
M = 10;
p = 0.3;

tvals = {'wedge','diamond'};
nsvals = 10.^(2:7);

clear relerr
for t = [1:2]
    for ns = 1:length(nsvals)
        fprintf('%s, %d samples\n',tvals{t},nsvals(ns));
        relerr(t,ns) = check_atb(K,M,[],p,nsvals(ns),producttype,matrixtype,tvals{t});        
    end    
end

figure(2);
semilogx(nsvals,relerr(1,:),nsvals,relerr(2,:))
legend('Wedges','Diamond');
xlabel('nsamples');
ylabel('relative error');
tstr = sprintf('Binary AtA for %d x %d matrix with %d%% nonzeros', K,M,round(p*100));
title(tstr);
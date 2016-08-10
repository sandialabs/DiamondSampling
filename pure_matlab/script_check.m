%% Test correctness of A'*B code

%% Run both wedge and diamond sampling, make some plots to show
%% error decreasing as nsamples increases

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
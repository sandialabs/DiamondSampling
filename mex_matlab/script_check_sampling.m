%% Test correctness of A'*B code
inputs = {'sparse','dense'};

for i = 1:length(inputs),
    sparsedense = inputs{i};
    
    % Set numbers of samples to test
    nsamples = 10.^(2:6);

    % Set up problem
    switch sparsedense
        case {'sparse'}
            M = 7; K = 5; N = 9; p = .8;
            A = sprandn(M,K,p);
            B = sprandn(K,N,p);
            tstr = sprintf('Sprandn data for %d x %d x %d multiplication with %d%% input nonzeros',M,K,N,round(p*100));
        case {'dense'}
            M = 5; K = 9; N = 7;
            A = randn(M,K);
            B = randn(K,N);
            tstr = sprintf('Dense randn data for %d x %d x %d multiplication',M,K,N);
%             K = 5; M = 3; N = 3;
%             [A ~] = qr(randn(K,M),0);
%             B = A';
%             tstr = sprintf('Dense orthogonal A^T*A data (%d x %d x %d multiplication)',M,K,N);

    end

    % Perform exact computations (P is E's structural pattern)
    E  = A*B;
    E2 = E.^2;
    Ep = abs(A)*abs(B);
    P  = (Ep~=0);
    wedge_total        = ones(1,M)*abs(A)*abs(B)*ones(N,1);
    three_path_B_total = ones(1,M)*abs(A)*abs(B)*abs(B)'*ones(K,1);
    three_path_A_total = ones(1,K)*abs(A)'*abs(A)*abs(B)*ones(N,1);   

    % Loop over numbers of samples, compute relative error
    re_wedge     = zeros(size(nsamples));
    re_wedge_t   = zeros(size(nsamples));
    re_diamond   = zeros(size(nsamples));
    re_diamond_t = zeros(size(nsamples));
    for j = 1:length(nsamples)

        % Perform sampling
        Xw  = sample_ab_mex(A,B,sparsedense,'wedge','optimized',nsamples(j)); 
        Xwt = sample_ab_mex(B',A',sparsedense,'wedge','naive',nsamples(j));
        Xwt = Xwt';
        Xd  = sample_ab_mex(A,B,sparsedense,'diamond','optimized',nsamples(j));
        Xdt = sample_ab_mex(B',A',sparsedense,'diamond','naive',nsamples(j));
        Xdt = Xdt';

        % Check for samples corresponding to structural zeros
        if ~isempty(find(Xw(P==0),1)) || ~isempty(find(Xwt(P==0),1))
            error('Structural zero found in X from wedge sampling!');
        end
        if ~isempty(find(Xd(P==0),1)) || ~isempty(find(Xdt(P==0),1))
            error('Structural zero found in X from diamond sampling!');
        end

        % Compute proportionality constants
        wc  = wedge_total        / nsamples(j);
        dc  = three_path_B_total / nsamples(j);
        dtc = three_path_A_total / nsamples(j);

        % Compute estimates (diamond sampling converges to entries squared)
        Yw  =  wc * Xw;
        Ywt =  wc * Xwt;
        Yd  =  dc * Xd;
        Ydt = dtc * Xdt;

        % Compute relative errors
        re_wedge(j)     = norm(Yw-E,'fro')  / norm(E,'fro');
        re_wedge_t(j)   = norm(Ywt-E,'fro') / norm(E,'fro');
        re_diamond(j)   = norm(Yd-E2,'fro')  / norm(E2,'fro');
        re_diamond_t(j) = norm(Ydt-E2,'fro') / norm(E2,'fro');

    end

    % Make plot
    figure(i); clf;
    loglog(nsamples,re_wedge,nsamples,re_wedge_t,nsamples,re_diamond,nsamples,re_diamond_t)
    legend('Wedge','Wedges (t)','Diamond','Diamond (t)');
    xlabel('nsamples');
    ylabel('relative error');
    title(tstr);
end


%----- Correlted Total Viariation based Tensor Robust Principle Component Analysis -----%
function [X, E, err, iter] = PTCTV_KMC(M, opts)
% Solve the p-order Tensor Robust Principle Component Analysis via Tensor Correlted Total Viariation(TCTV) norm minimization by ADMM
% the transform in high-order TSVD uses DFT (default)


%% default paremeters setting 
dim = size(M);
d   = ndims(M);

transform  = 'DFT';
for i = 3:d
transform_matrices{i-2} = dftmtx(dim(i)); 
end

% 超参数
lambda     = 0.85/sqrt(prod(dim)/max(dim(1),dim(2)));
% mu         = 0.11;
mu         = 0.15;
N = 1;

% lambda = 0.05;
directions = 1:2; 
% tol        = 0; 
tol        = 7e-3;
max_iter   = 500;
rho        = 1.05;
% rho        = 1.1;
% mu         = 7e-2;
% mu         = 9e-1;

max_mu     = 10;
detail     = 0;

if ~exist('opts', 'var')
    opts = [];
end   
if isfield(opts, 'tol');                tol                = opts.tol;                end
if isfield(opts, 'max_iter');           max_iter           = opts.max_iter;           end
if isfield(opts, 'rho');                rho                = opts.rho;                end
if isfield(opts, 'max_mu');             max_mu             = opts.max_mu;             end
if isfield(opts, 'detail');             detail             = opts.detail;             end

if isfield(opts, 'N');                 N                 = opts.N;                 end
if isfield(opts, 'mu');                 mu                 = opts.mu;                 end
if isfield(opts, 'lambda');             lambda             = opts.lambda/sqrt(prod(dim)/max(dim(1),dim(2)));             end
%  N  

%% variables initialization
n = length(directions);
X        = zeros(dim);
% E        = M-IT;
%________________
E        = zeros(dim);
Lambda   = zeros(dim);
 
for i = 1:n
    index        = directions(i);
    G{index}     = porder_diff(X,index); 
    Gamma{index} = zeros(dim); 
end

%% FFT setting
T = zeros(dim);
for i = 1:n
    Eny = diff_element(dim,directions(i));
    T   = T + Eny; 
end

%% main loop
iter = 0;
Xk = X;
Ek = E;

while iter<max_iter
    iter = iter + 1;  
    %% Update X -- solve TV by FFT 
	H = zeros(dim);
    for i = 1:n
       index = directions(i);
       H = H + porder_diff_T(mu*G{index}-Gamma{index},index); 
    end
    X = real( ifftn( fftn( mu*(M-Ek)+Lambda+H)./(mu*(1+T)) ) );
  
    %% Updata Gi -- proximal operator of TNN
    for i = 1:n
        index = directions(i);

        [G{index}] = prox_pstnn(porder_diff(X,index)+Gamma{index}/mu,N,1/(n*mu)); 

    end
    
    %% Update E 
	E = prox_l1(M-X+Lambda/mu, lambda/mu);
%     Weighten = 1./(abs(E)+0.001);
%     E(E<0)=0;
    
    %% Stop criterion
    dY   = M-X-E;
    err = norm(dY(:))/norm(M(:));
%     dY = Xk-X;
%     err = norm(dY(:))/norm(M(:));
%     dG = G{1} + G{2}- Gk{1}- Gk{2};
% 
%     GkM = Gk{1} + Gk{2};
%     err = norm(dG(:))/norm(GkM(:));
    errlist(iter) = err;
    if iter>1
        if err < tol
            break;
        end 
%         if nnz(E) == obj
%             break;
%         end 
    end
    obj = nnz(E);
    %% Update detail display
    if detail
        C = E;
        obj = nnz(E);
%         if iter == 1 || mod(iter, 10) == 0
        if iter > 1  
%             obj = sum(cell2mat(tnn_G))/n;

%             err = norm(dY(:),'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update mulipliers: Lambda, Gamma, and mu
    Lambda = Lambda+mu*dY;
    for i = 1:n
        index = directions(i); 
        Gamma{index} = Gamma{index}+mu*(porder_diff(X,index)-G{index});
    end
    mu = min(rho*mu,max_mu);    
    
    Ek = E;
    Xk = X;
end
save errlist.mat errlist
end


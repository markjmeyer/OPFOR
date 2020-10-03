function [ waic lppd DIC] = WAIC( Y, res )
%waic function to estimate WAIC, lppd, and DIC for OPFORs
%

%%
N   = size(res.model.X,1);
P   = size(res.model.X,2);
T   = size(res.bhat,2);
B   = size(res.MCMC_beta,1);
L   = length(unique(Y));
Bs  = NaN(B, T, P);
Ys  = NaN(N, T, L);

%%
for c = 1:P
    i           = c - 1;
    Bs(:,:,c)	= res.MCMC_beta(:,(i*T+1):(c*T));
    
    valy    = sprintf('Y%d', i);
    assign('caller', valy, Yc)
end

for l = 1:L
    i           = l - 1;
    Ys(:,:,l)	= 1*(Y == i);
end

LikMat      = zeros(B,N);





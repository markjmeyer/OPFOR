%%%%%%%%%%%%%%%%%%%%%%%% OPFOR: Sample Script %%%%%%%%%%%%%%%%%%%%%%%
% Demonstration of OPFOR                                            %
% Simulated dataset                                                 %
% Specifications: N = 40, peak curve, T = 365                       %
%                                                                   %
% Before running with O-Splines:                                    %
%           Place osullivan.R in ~/MATLAB folder                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add paths %%
% addpath('~/Code');
addpath('/Users/mjm556/Documents/Research/OPFM/Code')

%% paths to access R from matlab (only run once for O-Splines) %%
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1)

%% data parameters %%
N       = 40;
T       = 365;
t       = linspace(0,1,T);
pi_d    = 0.55;

%% simulate covariates and build beta(t) %%
rng(6915);
x1      = normrnd(0, 1, N, 1);
b1      = (2/5)*normpdf(t,0.6,0.15); % peak setting from manuscript

%% build model %%
model.X         = [ones(N,1) x1]; % set design matrix
model.alf       = 0.05;
model.delt      = norminv(pi_d);
model.B         = 1000; % number of retained samples per covariate

%% exponential %%
nu1     = 0.0049; % 0.0049
w1      = 0.33; % 0.33
a1      = 0.01; % 0.01
CovStr  = zeros(T);
for i = 1:T
    for j = 1:T
        CovStr(i,j) = nu1*exp(-0.5*w1*(t(i) - t(j))^2) + ((a1)^2)*t(i)*t(j);
    end
end

%% build errors %%
E       = NaN(N, T);
for n = 1:N
    E(n,:) = mvnrnd(zeros(1,T), CovStr);
end

%% generate probabilities %%
cuts    = [0 0.5 1];
Yst     = x1*b1 + E;
P0      = normcdf(cuts(1)-Yst);
P1      = normcdf(cuts(2)-Yst) - P0;
P2      = normcdf(cuts(3)-Yst) - normcdf(cuts(2)-Yst);
P3      = 1 - normcdf(cuts(3)-Yst);
        
%% simulate Y %%
Y   = NaN(N,T);
for n = 1:N
    for i = 1:T
        pni     = [P0(n,i) P1(n,i) P2(n,i) P3(n,i)];
        mult    = mnrnd(1, pni);
        Y(n,i)  = find(mult == 1)-1;
    end
end

%% run O-Spline model using defaults %%
model.basis = 'O-Spline'; % K = 4; Kp = 2;
res_OS4     = opfor(Y, model);

%% run B-Spline model using defaults %%
model.basis = 'B-Spline'; % K = 10; Kp = 2; degreee = 3; alps = 0.01;
res_BS10    = opfor(Y, model);

%% run B-Spline model using defaults %%
model.basis = 'wavelet'; % wavelet = sym8; wtmode = sym; nlevels = 6;
res_SYM6    = opfor(Y, model);

%% run O-Spline model setting specs %%
model.basis = 'O-Spline';
Bspecs.K    = 2;    % number of knots for main effects
Bspecs.Kp   = 2;    % number of FPCA basis functions to estimate
res_OS2     = opfor(Y, model, Bspecs);

%% run B-spline model setting specs %%
model.basis     = 'B-Spline';
Bspecs.K        = 5;    % number of knots for main effects
Bspecs.Kp       = 2;    % number of FPCA basis functions to estimate
Bspecs.degree   = 3;
Bspecs.alps     = 0.01;  % tuning parameter balancing 2nd-derivative and 0th-derivative penalties, Goldsmith & Kitago (2015)
res_BS5         = opfor(Y, model, Bspecs);

%% run wavelet model setting specs %%
model.basis     = 'wavelet';
Bspecs.wavelet  = 'sym8'; % wavelet basis function
Bspecs.wtmode   = 'sym'; % boundary padding
Bspecs.nlevels  = 8; % levels of decomposition
res_SYM8        = opfor(Y, model, Bspecs);

%% set res to model of interest %%
res         = res_OS4;
% res         = res_OS2;
% res         = res_BS10;
% res         = res_BS5;
% res         = res_SYM6;
% res         = res_SYM8;

%% extract covar-specific coefs %%
Bint        = res.MCMC_beta(:,1:T);
B1          = res.MCMC_beta(:,(T+1):(2*T));

%% lppd %%
B           = size(res.MCMC_beta,1);
LikMat      = zeros(B,N);

%% generate level-specific matrices %%
Y0  = 1*(Y == 0);
Y1  = 1*(Y == 1);
Y2  = 1*(Y == 2);
Y3  = 1*(Y == 3);

%%
for b = 1:B
    Betab   = [Bint(b,:)' B1(b,:)'];
    Yst     = model.X*Betab';

    %%
    cuts    = res.MCMCspecs.cuts;
    P0      = normcdf(cuts(1)-Yst);
    P1      = normcdf(cuts(2)-Yst) - P0;
    P2      = normcdf(cuts(3)-Yst) - normcdf(cuts(2)-Yst);
    P3      = 1 - normcdf(cuts(3)-Yst);
    
    %% WAIC %%
    L0      = Y0.*P0;
    L1      = Y1.*P1;
    L2      = Y2.*P2;
    L3      = Y3.*P3;

    %%
    L       = L0 + L1 + L2 + L3;
    
    %%
    LikMat(b,:)     = prod(L,2);

end

%%
LikMat(LikMat == 0) = min(LikMat(LikMat > 0));
lLikMat     = log(LikMat);

%%
lppd    = sum(log(mean(LikMat)));
fprintf('lppd = %f\n', lppd)
pwaic   = sum((1/(B-1))*sum((lLikMat - mean(lLikMat)).^2));
WAIC    = -2*lppd + 2*pwaic;
fprintf('WAIC = %f\n', WAIC)

%% for older MATLAB %%
% meanMat     = ones(size(lLikMat));
% for i = 1:B
%     meanMat(i,:) = mean(lLikMat);
% end

%%
% lppd    = sum(log(mean(LikMat)));
% fprintf('lppd = %f\n', lppd)
% pwaic   = sum((1/(B-1))*sum((lLikMat - meanMat).^2));
% WAIC    = -2*lppd + 2*pwaic;
% fprintf('WAIC = %f\n', WAIC)

%% DIC %%
Betab   = res.bhat;
Yst     = model.X*Betab;

%% probs %%
cuts    = res.MCMCspecs.cuts;
P0      = normcdf(cuts(1)-Yst);
P1      = normcdf(cuts(2)-Yst) - P0;
P2      = normcdf(cuts(3)-Yst) - normcdf(cuts(2)-Yst);
P3      = 1 - normcdf(cuts(3)-Yst);

%% likelihood %%
Lik    = zeros(size(Y));
for i = 1:N
    for j = 1:T
        %%
        Pi          = [P0(i,j)' P1(i,j)' P2(i,j)' P3(i,j)'];
        Yi          = ismember(min(min(Y)):max(max(Y)), Y(i,j));
        Liki        = mnpdf(Yi, Pi);
        Lik(i,j)	= Liki;
    end
end

%%
lLik    = log(Lik);
DIC     = -2*sum(sum(lLik)) + 2*var(sum(lLik,2));
fprintf('DIC = %f\n', DIC)

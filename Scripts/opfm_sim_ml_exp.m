%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OPFM SIM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPFM Simulation based on Generated data for Odyssey              %
% Monkey Based, burnin 500, sample 500                             %
% Wavelet Bases: Coiflets, m = 4, J = 6                            %
%                Symmlets, m = 8, J = 6                            %
% Spline Bases:  B-spline, K = 5, 7                                %
%                O-Spline, K = 4, 6                                %
% x ~ N(0,1), L = 3, N = 40, T = 365                               %
% Independence Latent Covariance Structure                         %
%                                                                  %
% Created:      10/22/2018                                         %
% Modified:     10/24/2018                                         %
%                                                                  %
% By:           MJ Meyer                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add paths %%
% addpath('~/Documents/Code/OPFM/Code');
addpath('/Users/mjm556/Documents/Research/OPFM/Code')

%% data parameters %%
N       = 40;
T       = 365;
t       = linspace(1,108,T);
pi_d    = 0.55;

%% simulate covariates and build beta(t) %%
rng(6915);
x1      = normrnd(0, 1, N, 1);
b1      = (1/7)*exp(2*(1-(t./127)));

%% build model %%
model.X         = x1;
model.alf       = 0.05;
model.delt      = norminv(pi_d);

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

%% loop over settings %%
for b = 5:6
    if b == 1
        model.basis     = 'wavelet';
        Bspecs.nlevels  = 6;
    elseif b == 2
        model.basis     = 'wavelet';
        Bspecs.nlevels  = 8;
    elseif b == 3
        model.basis     = 'B-Spline';
        Bspecs.K        = 5;    % number of knots for main effects
    elseif b == 4
        model.basis     = 'B-Spline';
        Bspecs.K        = 10;    % number of knots for main effects
    elseif b == 5
        model.basis     = 'O-Spline';
        Bspecs.K        = 2;    % number of knots for main effects
    else
        model.basis     = 'O-Spline';
        Bspecs.K        = 4;    % number of knots for main effects
    end
    
    %% loop over datasets %%
    if b < 3
        for seed = 1:200
            %% set seed %%
            rng(seed);

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

            %% run simulated dataset
            res             = opfm_sim(Y, model, Bspecs);
    
            %% save output
            if b == 1
                fname   = sprintf('/Users/mjm556/Documents/Research/OPFM/Simulation/ML/Exponential/Symmlet6/mlIndB%dS%d.mat',b,seed);
            elseif b == 2
                fname   = sprintf('/Users/mjm556/Documents/Research/OPFM/Simulation/ML/Exponential/Symmlet8/mlIndB%dS%d.mat',b,seed);
            end
            save(fname, 'res');
    
            %% clear current model
            clear res
        end
    else
        %%
        count   = 1;
        seed    = 0;
        while count < 201
            %% set seed %%
            seed    = seed + 1;
            rng(seed);
            
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

            %% run simulated dataset %%
            try
                res             = opfm_sim(Y, model, Bspecs);
            catch ME
                continue
            end
            count   = count + 1;
            
            %% save output %%
            if b == 3
                fname   = sprintf('/Volumes/G-DRIVE mobile/Research/Probit/Simulation/ML/Exponential/Bspline5/mlExpB%dS%d.mat',b,seed);
            elseif b == 4
                fname   = sprintf('/Volumes/G-DRIVE mobile/Research/Probit/Simulation/ML/Exponential/Bspline10/mlExpB%dS%d.mat',b,seed);
            elseif b == 5
                fname   = sprintf('/Volumes/G-DRIVE mobile/Research/Probit/Simulation/ML/Exponential/Ospline4/mlExpB%dS%d.mat',b,seed);
            else
                fname   = sprintf('/Volumes/G-DRIVE mobile/Research/Probit/Simulation/ML/Exponential/Ospline6/mlExpB%dS%d.mat',b,seed);
            end
            save(fname, 'res');
    
            %% clear current model
            clear res
        end
    end
end

%%


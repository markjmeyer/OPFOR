function [ res ] = opfm_sim( Y, model, Bspecs )
%opfm_sim Function to run simulation settings
%   for OPFM

%% General MCMCspecs
MCMCspecs.thin              = 1;
MCMCspecs.time_update       = 100;
MCMCspecs.a                 = 0;
MCMCspecs.saveSigma         = 1;
MCMCspecs.updateLVvar       = 0;    % has some issues
MCMCspecs.updateCuts        = 0;
MCMCspecs.cuts              = [MCMCspecs.a 0.5 1]; %[MCMCspecs.a 1];

%% General model add-ons %%
model.C         = ones(size(Y,1),1);    % for wavelets only

%% set up basis specifications %%
switch model.basis
    case 'wavelet'
        %% basis-specific model add-ons %%
        Z   = [];   % use fixed effects only, required for sampler components of wavelets?
        model.Z{1}      = Z;                    % for wavelets only
        model.H         = 0;                    % for wavelets only
        model.Hstar     = 0;                    % for wavelets only
        model.m         = [];                   % for wavelets only

        %% set chain length %%
        MCMCspecs.B                 = size(model.X,2)*1000;
        MCMCspecs.burnin            = MCMCspecs.B;
        
        %% basis-specific MCMCspecs %%
        MCMCspecs.nj_nosmooth       = 1;        % can be 0, if is 0, don't do any constraint for pi_{ij}, bigTau_{ij}.
        MCMCspecs.minp              = 1e-14;
        MCMCspecs.maxO              = 1e20;
        MCMCspecs.minVC             = 1e-20;
        MCMCspecs.VC0_thresh        = 1e-4;
        MCMCspecs.delta_theta       = 1e-4;
        MCMCspecs.thetaMLE_maxiter  = 10^6;
        MCMCspecs.EmpBayes_maxiter  = 10^6;
        MCMCspecs.propsdTheta       = 0.015;
        MCMCspecs.tau_prior_var     = 1e3;      % the variance of tau_{ijk} when finding prior parameters for tau_{ijk}.
        MCMCspecs.tau_prior_idx     = 1;        % 1 indicate that a_tau and b_tau depend on ij, 0 indicate that they depend on jk. 
        MCMCspecs.PI_prior_var      = 0.06;     % this range should be in [0.02 0.09].

        %% set Bspecs to default %%
        if nargin < 3
            Bspecs.wavelet          = 'sym8';   % coif4, sym8
            Bspecs.wtmode           = 'sym';    % sym
        end
    case 'B-Spline'
        %% set chain length %%
        MCMCspecs.B             = size(model.X,2)*2000;
        MCMCspecs.burnin        = MCMCspecs.B;
        
        %% basis-specific MCMCspecs %%
        MCMCspecs.Asig          = 1;
        MCMCspecs.Bsig          = 1;
        MCMCspecs.sig2me        = 0.01; % starting value for measurement error variance, Goldsmith & Kitago (2015)
        
        %% set Bspecs %%
        if nargin < 3
            Bspecs.Kp               = 2;    % number of FPCA basis functions to estimate
            Bspecs.degree           = 3;
            Bspecs.alps             = 0.01;  % tuning parameter balancing 2nd-derivative and 0th-derivative penalties, Goldsmith & Kitago (2015)
        end
    case 'O-Spline'
        %% set chain length %%
        MCMCspecs.B             = size(model.X,2)*2500;
        MCMCspecs.burnin        = MCMCspecs.B;
        
        %% basis-specific MCMCspecs %%
        MCMCspecs.Asig          = 1;
        MCMCspecs.Bsig          = 1;
        MCMCspecs.sig2me        = 0.01; % starting value for measurement error variance, Goldsmith & Kitago (2015)
        
        %% set Bspecs %%
        if nargin < 3
            Bspecs.Kp               = 2;    % number of FPCA basis functions to estimate
        end
end

%% run model %%
res     = opfm(Y, model, Bspecs, MCMCspecs);


end


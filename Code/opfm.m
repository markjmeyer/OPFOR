function res = opfm(Y, model, Bspecs, MCMCspecs) 
%OPFM Performs wavelet or spline based ordinal probit function-on-scalar 
%           regression for categorical outcome with any number of
%           ordered levels. Y may be binary, three level, four level, etc.
%           Allows for any number of scalar covariates.
%
%   Created:    09/26/2016
%   Modified:   10/22/2018
%
%   By:         MJ Meyer

%% get model parameters %%
[ n, p ]    = size(model.X);   %%% n=# functions, p=# covariates
model.c     = size(model.C,2);
model.p     = p;
model.n     = n;
O           = size(Y,2); %% dim of Y
X           = model.X;
model.O     = O;
Ylevels     = unique(Y);
if length(Ylevels) > 2
    isOrdinal   = 1;
else
    isOrdinal               = 0;
    MCMCspecs.updateCuts    = null(1);
end
basis       = model.basis;

%% set initial values and priors for latent variables, cuts, and beta
beta        = (X'*X)^(-1)*X'*Y; % starting values for beta, naive OLS
if isOrdinal == 1
    cuts        = MCMCspecs.cuts;
else
    cuts        = 0;
end
if MCMCspecs.updateLVvar
    Sigma       = ones(O,1);
end

%% extract MCMC specs
B           = MCMCspecs.B;
burnin      = MCMCspecs.burnin;
thin        = MCMCspecs.thin;
if isempty(MCMCspecs.updateCuts)
    updateCuts  = 0;
else
    updateCuts  = MCMCspecs.updateCuts;
end

%% matrix to save beta samples
MCMC_beta	= NaN(B,p*O);
if updateCuts == 1
    MCMC_cuts           = NaN(B,length(cuts)-1);
end
if MCMCspecs.saveSigma
    MCMC_Sigma          = NaN(B,O);
end

%%
ii          = 0; % counter for MCMC samples to output
tic

%% MCMC Loop: B is desired # samples
for i = 1:(B*thin+burnin)
    
    %% (1) update latent variable
    if MCMCspecs.updateLVvar == 1
        [ yStar ]   = GenUpdateLV(beta, Y, model, cuts, Sigma);
    else
        [ yStar ]   = GenUpdateLV(beta, Y, model, cuts);
    end
    
    %% (2) update cutpoints (if desired)
    if updateCuts == 1 && isOrdinal == 1
        [ b ]           = GenUpdateCuts(yStar, Y, cuts);
        cuts(2:end)     = b;
    end
    
    %% (3) reshape latent variable, project into wavelet-space
	yStar           = reshape(yStar,n,O);

    %%
	if i == 1
        switch basis
            case 'wavelet'
                % (3.1) use WavInitVals to get initial values on first iteration
                [ initial ]     = WavInitVals(yStar, Bspecs, model, MCMCspecs);
                bstar           = initial.bstar;    % 1 x # basis coef
                Vbetans         = initial.Vbetans;  % 1 x # basis coef
                PiMat           = initial.PiMat;    % 1 x # basis coef
                TauMat          = initial.TauMat;   % 1 x # basis coef
                Wv              = initial.Wv;
                Bspecs          = initial.wavespecsy;
                theta           = initial.theta;    % 1 x # basis coef
                D               = initial.D;        % n x # basis coef
                W               = initial.W;
                prior_Theta_a   = initial.prior_Theta_a;    % 1 x # basis coef
                prior_Theta_b   = initial.prior_Theta_b;    % 1 x # basis coef
                propsd_Theta    = initial.propsd_Theta;     % 1 x # basis coef
                K               = Bspecs.K;         % dim of basis coef
                MCMC_alpha      = NaN(B,p*K);
                MCMC_MHflag     = NaN(B,K);
        
                % for sampling tau and pi
                a_tau           = initial.a_tau;    % 1 x levels of wav decomp
                b_tau           = initial.b_tau;    % 1 x levels of wav decomp
                a_pi            = initial.a_pi;     % 1 x levels of wav decomp
                b_pi            = initial.b_pi;     % 1 x levels of wav decomp
                meanop          = initial.meanop;   % # basis coef x levels wd
                expandop        = initial.expandop; % levels wd x # basis coef
                MCMC_tau        = NaN(B,p*Bspecs.J);% Bspecs.J = levels wd
                MCMC_pi         = NaN(B,p*Bspecs.J);
            otherwise
                %% code for B-spline and O'Sullivan spline initialization %
                [ initial ]     = BSplineInitVals(yStar,Bspecs,model);
                Pmat            = initial.Pmat;
                Bspecs          = initial.Bspecs;
                    K           = Bspecs.K;
                    Kp          = Bspecs.Kp;
                    Theta       = Bspecs.Theta;
                sigW            = initial.sigW;
                ThtTh           = initial.ThtTh;
                WdkTheta        = initial.WdkTheta;
                Yvec            = initial.Yvec;
                sig2me          = MCMCspecs.sig2me;
                Wdes            = model.X;
                
                % for sampling lambdas
                lambdaInit.Aw       = initial.Aw;
                lambdaInit.Bw       = initial.Bw;
                lambdaInit.Apsi     = initial.Apsi;
                lambdaInit.Bpsi     = initial.Bpsi;
                lambdaInit.Asig     = MCMCspecs.Asig;
                lambdaInit.Bsig     = MCMCspecs.Bsig;
                lambdaInit.Kp       = Kp;
                
                % initialize matrices and starting values
                BW              = zeros(K,p,B);
                bw              = BW(:,:,1);
                BPSI            = zeros(K,Kp,B);
                bpsi            = BPSI(:,:,1);
                MCMC_C          = zeros(n,Kp,B);
                cmat            = normrnd(0, 0.01, n, Kp);
                MCMC_C(:,:,1)   = cmat;
                MCMC_Sigma      = sig2me*ones(B,1);
                MCMC_lbw        = ones(B,p);
                lambdabw        = MCMC_lbw(1,:);
                MCMC_lpsi       = ones(B,Kp);
                lambdapsi       = MCMC_lpsi(1,:);
                psicur          = bpsi'*Theta';
                pcaefcur        = cmat*psicur;
                MCMC_psi        = NaN(B,Kp*O);
        end
        fprintf('\n Beginning Sampler \n \n');
	else
        switch basis
            case 'wavelet'
                % (3.2) update W and D using new Y*
                [D, ~]      = dwt_rows(yStar, Bspecs);
                W           = GetW(model,D);
            otherwise
                Yvec        = reshape(yStar',[],1);
                psicur      = bpsi'*Theta';
                pcaefcur	= cmat*psicur;
        end
	end
    
    %% (4) update model parameters
    switch basis
        case 'wavelet'
        %% (4.1) Update beta
            [bstar,gamma,alpha]     = UpdateBetaNoOrthog(bstar,Vbetans,PiMat,TauMat,Wv,model,Bspecs,MCMCspecs); 
            Wv.L2                   = Get_L2(bstar,Wv,Bspecs); %#zhu#% beta is updated here, so does L2. 

        %% (4.2) Update theta. q_jk, and s_jk
            [theta,MHflag,Vbetans,Wv]    = UpdateTheta(bstar,theta,Vbetans,Wv,D,W,model,prior_Theta_a,prior_Theta_b,propsd_Theta,Bspecs,MCMCspecs);   

        %% (4.3) Update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat)
            tau     = Update_tau(bstar,gamma,a_tau,b_tau,meanop);
            TauMat  = tau*expandop./Vbetans;
            PI      = Update_pi(gamma,a_pi,b_pi,Bspecs);
            PiMat   = PI*expandop;
        otherwise
        %% (4.1) Update BW            
            bw      = UpdateBW(Yvec,WdkTheta,pcaefcur,sig2me,sigW,lambdabw,Pmat,model);
            
        %% (4.2) Update Bpsi
            betacur     = bw'*Theta';
            fixefcur    = Wdes*betacur;
            bpsi        = UpdateBPsi(Yvec,fixefcur,Kp,sig2me,cmat,ThtTh,lambdapsi,Pmat,Theta);
            
        %% (4.3) Update theta
            cmat        = UpdateC(yStar,cmat,fixefcur,bpsi,sig2me,Kp,Theta,model);
            
        %% (4.4) Update tau
            [ sig2me, lambdabw, lambdapsi ] = UpdateLambdas(yStar,fixefcur,lambdabw,lambdapsi,lambdaInit,bw,bpsi,Theta,cmat,Pmat,Bspecs,model);
    end
    
	%% (5) project beta* and theta into data-space
    switch basis
        case 'wavelet'
            beta            = idwt_rows(bstar,Bspecs);
            [~,Sigma,~]     = GetSigma(theta,Bspecs);
        otherwise
            beta            = bw'*Theta';
            psi             = bpsi'*Theta';
    end
    
    %%  Record MCMC samples %%
    if   (i > burnin) && (mod(i-burnin,thin) == 0)   %% Save MCMC samples of beta in single matrix.
            ii                      = ii+1; % this is the real row number among the B samples.
            is                      = mod(ii-1,B)+1; % this is the row number in the block.

            %% save betas in data-space for easier post-processing
            MCMC_beta(is,:)         = reshape(beta',1,p*O); %#zhu#% each row=[beta(1,j,k),j,k|beta(2,j,k),j,k|...|beta(p,j,k),j,k]
            if updateCuts == 1
                MCMC_cuts(is,:)     = cuts(2:end);
            end
            
            %% save wavelet-space parameters
            switch basis
                case 'wavelet'
                    %% save wavelet-space parameters
                    MCMC_alpha(is,:)        = reshape(alpha',1,p*K);
                    MCMC_pi(is,:)           = reshape(PI,1,p*Bspecs.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values.
                    MCMC_tau(is,:)          = reshape(tau,1,p*Bspecs.J);   % each block is one j, contains p values.

                    %% save variance components
                    MCMC_Sigma(is,:)        = Sigma;
                    MCMC_MHflag(is,:)       = MHflag;
                otherwise
                    %% save spline-based components
                    MCMC_psi(is,:)          = reshape(psi',1,Kp*O);
                    MCMC_lbw(is,:)           = lambdabw;
                    MCMC_lpsi(is,:)          = lambdapsi;
                    MCMC_Sigma(is)           = sig2me;
                    MCMC_C(:,:,is)           = cmat;
            end
    end
    if mod(i, 10) == 0
        fprintf('.')
    end
    if mod(i, MCMCspecs.time_update) == 0
       fprintf('\n %d \n',i),toc;
    end
    
end

fprintf('\n Done with MCMC \n');

%% perform posterior functional inference
[psi, pst]                      = BFDR(MCMC_beta, model.delt, model.alf);
[SimBaS, upper_CI, lower_CI]    = jointband_sbs(MCMC_beta, model.alf);

%% save MCMC samples
res.MCMC_beta   = MCMC_beta;
switch basis
    case 'wavelet'
        res.MCMC_tau    = MCMC_tau;
        res.MCMC_alpha  = MCMC_alpha;
        res.MCMC_pi     = MCMC_pi;
        res.pihat       = reshape(mean(MCMC_pi)',p,Bspecs.J);
        res.tauhat      = reshape(mean(MCMC_tau),p,Bspecs.J);
        res.MCMC_Sigma  = MCMC_Sigma;        
        res.MCMC_MHflag = mean(MCMC_MHflag);
    otherwise
        res.MCMC_psi    = MCMC_psi;
        res.MCMC_lbw    = MCMC_lbw;
        res.MCMC_lpsi   = MCMC_lpsi;
        res.MCMC_Sigma  = MCMC_Sigma;
        res.MCMC_C      = MCMC_C;
end

if updateCuts == 1
    res.MCMC_cuts   = MCMC_cuts;
end
res.model       = model;
res.Bspecs      = Bspecs;
res.MCMCspecs   = MCMCspecs;

%% reshape and store model estimates %%
res.bhat        = reshape(mean(MCMC_beta),O,p)';
res.Q025_bhat   = reshape(quantile(MCMC_beta,0.025),O,p)';
res.Q975_bhat   = reshape(quantile(MCMC_beta,0.975),O,p)';
if updateCuts == 1
    res.cuts        = [cuts(1) mean(MCMC_cuts)];
end

fprintf('\n Done storing model estimates \n');

%% store PFI results %%
res.SimBaS      = reshape(SimBaS,O,p);
res.USimBaS     = reshape(upper_CI,O,p);
res.LSimBaS     = reshape(lower_CI,O,p);
res.psi         = reshape(psi',O,p);
res.pst         = reshape(pst',O,p);

fprintf('\n Done with PFI \n');

%% latent variable and initial values %%
res.yStar       = yStar;    % last yStar
switch basis
    case 'wavelet'
        res.Vbetans     = Vbetans;  % last Vbetans
        res.D           = D;        % last yStar*BasisMatrix
        res.Wv          = Wv;       % last Wv
    otherwise
        res.bw          = bw;       % last basis-space fixed effects
        res.bspi        = bpsi;     % last basis-space FPCA effects
end
res.initial     = initial;

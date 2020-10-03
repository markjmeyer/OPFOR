function [ initial ] = BSplineInitVals( yStar, Bspecs, model )
%BSPLINEINITVALS2 Find inital values for coefficient sampler in OPFM
%
%   Created:    10/12/2018
%   Modified:   10/14/2018
%   By:         MJ Meyer
tic

%%
K       = Bspecs.K;
basis   = model.basis;
Wdes    = model.X;
D       = size(yStar,2);
p       = size(Wdes,2);

%% set spline basis functions %%
fprintf('\n Constructing Basis Functions and Penalty Matrix \n')
switch basis
    case 'B-Spline'
        knots           = [zeros(1,3), linspace(0,1,K-2), ones(1,3)];
        Bord            = Bspecs.degree+1;  % MATLAB parameterizes B-Splines in terms of order, 1 + degree
        Theta           = spcol(knots, Bord, linspace(0,1,D));
        Bspecs.GKPen    = 1;
    case 'O-Spline'
        %% save osullivan.R inputs for flat file communication %%
        cd('/Users/mjm556/Documents/MATLAB')
        save('K.mat','K')
        save('D.mat','D')
        
        %% paths to access R from matlab (only run once) %%
%         path1 = getenv('PATH');
%         path1 = [path1 ':/usr/local/bin'];
%         setenv('PATH', path1)

        %% call R from matlab %%
        ! R CMD BATCH osullivan.R outfile.txt
        
        %% pull R output into matlab %%
        Theta           = load('Theta.mat');
        Theta           = Theta.Theta;
        K               = size(Theta,2);    % update K, will be original K + 2
        Bspecs.K        = K;                % update Bspecs
        Bspecs.GKPen    = 0;
        
        %% clear path %%
        clear path1;
end
Bspecs.Theta    = Theta;

%% construct penalty matrix %%
if Bspecs.GKPen == 1
    %% GK penalty matrix %%
    diff0   = eye(D);
    d2temp  = eye(D) - 2*diag(ones(D-1,1),1) + diag(ones(D-2,1),2);
    diff2   = d2temp(1:(D-2),:);

    P0      = Theta'*(diff0'*diff0)*Theta;
    P2      = Theta'*(diff2'*diff2)*Theta;

    alps    = Bspecs.alps;
    Pmat    = alps*P0 + (1-alps)*P2;
else
    %% O'Sullivan penalty matrix %%
    Omega   = load('Omega.mat');
    Pmat   = Omega.Omega;
end

%% construct matrices %%
Yvec        = reshape(yStar',[],1);
WdkTheta    = kron(Wdes,Theta)';
WtW         = Wdes'*Wdes;
ThtTh       = Theta'*Theta;
sigW        = kron(WtW,ThtTh);
varW        = sigW^(-1);
vecBW       = varW*WdkTheta*Yvec;

%% set hyper-parameters %%
muqBW       = reshape(vecBW, [], p);
Aw          = K/2;
u           = 1:p;
Bw          = arrayfun(@(x) max(1, 0.5*trace((muqBW(:,x)')*Pmat*muqBW(:,x))), u);
Apsi        = K/2;
Bpsi        = K/2;

%% store output %%
initial.Pmat        = Pmat;
initial.Bspecs      = Bspecs;
initial.sigW        = sigW;
initial.ThtTh       = ThtTh;
initial.WdkTheta    = WdkTheta;
initial.Yvec        = Yvec;
initial.Aw          = Aw;
initial.Bw          = Bw;
initial.Apsi        = Apsi;
initial.Bpsi        = Bpsi;
fprintf('\n Starting Values Initialized \n \n'),toc;

end


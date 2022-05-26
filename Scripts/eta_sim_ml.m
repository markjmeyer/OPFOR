%% paths to access R from matlab (only run once) %%
path1 = getenv('PATH');
path1 = [path1 ':/usr/local/bin'];
setenv('PATH', path1)

%% data parameters %%
N       = 40;
T       = 365;
t       = linspace(1,108,T);

%% simulate covariates and build beta(t) %%
rng(6915);
x1      = normrnd(0, 1, N, 1);
b1      = (1/7)*exp(2*(1-(t./127)));

%% exponential %%
nu1     = 0.0049; % 0.0049
w1      = 0.33; % 0.33
a1      = 0.01; % 0.01
CovStrE  = zeros(T);
for i = 1:T
    for j = 1:T
        CovStrE(i,j) = nu1*exp(-0.5*w1*(t(i) - t(j))^2) + ((a1)^2)*t(i)*t(j);
    end
end

%% independence %%
CovStrI         = diag(diag(CovStrE));

%% compound symmetric %%
R               = 0.5*ones(T);
R               = R.*~eye(T) + eye(T);
D               = CovStrI;
CovStrC         = D*R*D;

%% set CovStr %%
CovStr          = NaN(T,T,3);
CovStr(:,:,1)   = CovStrE;
CovStr(:,:,2)   = CovStrI;
CovStr(:,:,3)   = CovStrC;

%% number of simulated data sets %%
ndata       = 200;

%% construct matrices %%
mise        = NaN(ndata, 3);
runt        = NaN(ndata, 3);

%% loop over settings %%
for c = 1:3
    for seed = 1:ndata
        %% set seed %%
        rng(seed);

        %% build errors %%
        E       = NaN(N, T);
        for n = 1:N
            E(n,:) = mvnrnd(zeros(1,T), CovStr(:,:,c));
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

        %% save osullivan.R inputs for flat file communication %%
        cd('~/MATLAB')
        save('Y.mat','Y')
        save('X.mat','x1')

        %% call R from matlab %%
        ! R CMD BATCH mixed.R outfile.txt

        %% pull R output into matlab %%
        et              = load('eta.mat');
        eta             = et.eta;
        rtime           = load('rtime.mat');
        rt              = rtime.rtime;

        %% MISE %%
        mise(seed,c)	= mean(mean((eta - x1*b1).^2));    
        
        %% time %%
        runt(seed,c)    = rt;
    end
end

%% save files %%
save('etamiseML.mat', 'mise')
save('etaruntML.mat', 'runt')

function [ yStar ] = GenUpdateLV( beta, Y, model, cuts, Sigma )
%UPDATELV Update Latent Variable for OPWAVFM using truncated normals.
%           Normals truncated based on values of cuts. Y levels must be
%           of the form Y = {0, 1, ..., L} with 0 as the lowest level.
%           Binary also allowed. Non binary is assumed ordinal.
%
%   Created: 09/26/2016
%   By: MJ Meyer

%% get parameters %%
X           = model.X;
N           = model.n;
O           = model.O;
Ylevels     = unique(Y);

%%
if nargin < 5
    Sigvec  = 1; % for when we don't sample Sigma
else
    Sigvec  = repmat(Sigma, N, 1)';
end

%% find E(Y*)
EY  = X*beta;

%% vectorize Y and EY, set up yStar %%
Yvec    = reshape(Y,1,N*O);
EYvec   = reshape(EY,1,N*O);
yStar   = NaN(size(Yvec));
Nvec    = size(Yvec,2);

%% sample from truncated normals %%

%% Y = 0 %
rand0   = 0 + (normcdf(cuts(1),EYvec,Sigvec)-0).*rand(1,Nvec);
norm0   = norminv(rand0,EYvec,Sigvec);
yStar(Yvec == Ylevels(1))       = norm0(Yvec == Ylevels(1));

%% 0 < Y < L %
if length(Ylevels) > 2
    for c = 1:Ylevels(end-1)
        randc   = normcdf(cuts(c),EYvec,Sigvec) + (normcdf(cuts(c+1),EYvec,Sigvec)-normcdf(cuts(c),EYvec,Sigvec)).*rand(1,Nvec);
        normc   = norminv(randc,EYvec,Sigvec);
        yStar(Yvec == c)        = normc(Yvec == c);
    end
end

%% Y = L %
randL   = normcdf(cuts(end),EYvec,Sigvec) + (1-normcdf(cuts(end),EYvec,Sigvec)).*rand(1,Nvec);
normL   = norminv(randL,EYvec,Sigvec);
yStar(Yvec == Ylevels(end))     = normL(Yvec == Ylevels(end));

end


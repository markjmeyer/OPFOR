function [ b ] = GenUpdateCuts( yStar, Y, cuts )
%UPDATECUTS Update cutpoints (a, b) for OPWAVFM first fixing a
%           for identifiability then sampling b from Unif(c1,c2)
%
%   Created 09/26/2016
%   By: MJ Meyer

%% vectorize Y %%
N           = size(Y,1);
O           = size(Y,2);
Yvec        = reshape(Y,1,N*O);
Ylevels     = unique(Y);
b           = NaN(1,length(cuts)-1);
g           = [cuts inf];

%% loop over cut points to sample %%
for j = 2:length(cuts)
    %%
    l   = Ylevels(j);
    
    %% find c1, c2 using yStar, Yvec, and a (usually set a = 0) %%
    c1      = max(max(yStar(Yvec == l)), g(j-1));
    c2      = min(min(yStar(Yvec == l+1)), g(j+1));

    %% sample b|y*,y,a ~ U(c1,c2) %%
    b(j-1)  = c1 + (c2-c1)*rand;
    g(j)    = b(j-1);
end

end


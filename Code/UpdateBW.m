function [ bw ] = UpdateBW( Yvec, WdkTheta, pcaefcur, sig2me, sigW, lambdabw, Pmat, model )
%UpdateBW Update BW for spline-based model
%
%   Created:    10/14/2018
%   Modified:   10/14/2018
%
%   By:         MJ Meyer
    
    p           = model.p;
    meancur     = reshape(pcaefcur,[],1);
    presc       = (1/sig2me)*sigW + kron(diag(1./lambdabw),Pmat);
    sigma       = presc^(-1);
    mu          = (1/sig2me)*sigma*WdkTheta*(Yvec - meancur);
    bw          = reshape(mvnrnd2(mu,sigma),[],p);

end


function [ bpsi ] = UpdateBPsi( Yvec, fixefcur, Kp, sig2me, cmat, ThtTh, lambdapsi, Pmat, Theta )
%UNTITLED3 Update BPsi for spline-based model
%
%   Created:    10/14/2018
%   Modified:   10/14/2018
%
%   By:         MJ Meyer

    meancur     = reshape(fixefcur',[],1);
    presc       = (1/sig2me)*kron(cmat'*cmat, ThtTh) + kron(diag(1./lambdapsi),Pmat);
    sigma       = presc^(-1);
    CkTht       = kron(cmat,Theta)';
    mu          = (1/sig2me)*sigma*CkTht*(Yvec - meancur);
    bpsi        = reshape(mvnrnd2(mu,sigma),[],Kp);

end


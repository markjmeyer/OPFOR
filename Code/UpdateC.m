function [ cmat ] = UpdateC( yStar, cmat, fixefcur, bpsi, sig2me, Kp, Theta, model )
%UpdateC Update C matrix of FPCA basis functions
%
%   Created:    10/14/2018
%   Modified:   10/14/2018
%   By:         MJ Meyer

    n           = model.n;
    psicur      = bpsi'*Theta';
    ppT         = psicur*psicur';
    for c = 1:n
        presc       = (1/sig2me)*ppT + eye(Kp);
        sigma       = presc^(-1);
        mu          = (1/sig2me)*sigma*psicur*(yStar(c,:) - fixefcur(c,:))';
        cmat(c,:)   = mvnrnd(mu,sigma);
    end

end


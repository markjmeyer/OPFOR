function [ sig2me, lambdabw, lambdapsi ] = UpdateLambdas( yStar, fixefcur, lambdabw, lambdapsi, lambdaInit, bw, bpsi, Theta, cmat, Pmat, Bspecs, model  )
%UpdateLambdas Update Lambda smoothing parameters
%   for splined-based models
%
%   Created:    10/12/2018
%   Modified:   10/14/2018
%
%   By:         MJ Meyer

    Asig        = lambdaInit.Asig;
    Bsig        = lambdaInit.Bsig;
    Aw          = lambdaInit.Aw;
    Bw          = lambdaInit.Bw;
    Apsi        = lambdaInit.Apsi;
    Bpsi        = lambdaInit.Bpsi;
    Kp          = lambdaInit.Kp;
    

    n           = model.n;
    D           = model.O;
    K           = Bspecs.K;
    p           = model.p;
    psicur      = bpsi'*Theta';
    pcaefcur    = cmat*psicur;
    Ycur        = fixefcur + pcaefcur;
    Yres        = reshape((yStar - Ycur),[],1);
    YrtYr       = (Yres')*Yres;
    apost       = Asig + (n*D)/2;
    bpost       = Bsig + 1/2*YrtYr;
    sig2me      = 1/gamrnd(apost, 1/bpost);
    for j = 1:p
        apost           = Aw + K/2;
        bpost           = Bw(j) + (1/2)*(bw(:,j)')*Pmat*bw(:,j);
        lambdabw(j)     = 1/gamrnd(apost,1/bpost);
    end
    for k = 1:Kp
        apost           = Apsi + K/2;
        bpost           = Bpsi + (1/2)*(bpsi(:,k)')*Pmat*bpsi(:,k);
        lambdapsi(k)    = 1/gamrnd(apost,1/bpost);
    end


end


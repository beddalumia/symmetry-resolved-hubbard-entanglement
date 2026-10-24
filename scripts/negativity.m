%% Compute negativity for a two-orbital RDM in usual EDIPACK ordering
%
%   >> [logNB,logNF,NB,NF] = negativity(RDMij)
%
%      - logNB: Logarithmic negativity (standard "bosonic" partial transpose)
%      - logNF: Logarithmic negativity (twisted fermionic  partial transpose)
%      - NB:    "Original negativity"  (standard "bosonic" partial transpose)
%      - NF:    "Original negativity"  (twisted fermionic  partial transpose)
%
%  © Gabriele, Bellomia, 2025
function [logNB,logNF,NB,NF] = negativity(RDMij)
    if sum(sum(imag(RDMij)))>1e-8
    disp(imag(RDMij))
    error
    end
    % Preprocess the RDM
    [FiRDM, TiRDM,tindex] = partial_transpose(RDMij);
    % Get the eigenvalues
    p = eig(TiRDM);
    % Test for negativity
    NB = -sum(p(p<0));
    NB_= (trace(sqrtm(TiRDM'*TiRDM))-1)/2;
    if(abs(NB-NB_)>1e-8);
        warning('Negativity definitions mismatch: %g vs %g', NB, NB_);
    end
    p = eig(FiRDM);
    % Test for twisted vs untwisted negativity
    NF = -sum(p(p<0)) + (trace(FiRDM) - 1)/2;
    %NF_ = (trace(sqrtm(FiRDM'*FiRDM))-1)/2;
    NF_ = (sum(svd(FiRDM))-1)/2;
        if(abs(NF-NF_)>1e-8);
        warning('Negativity definitions (Fermionic) mismatch: %g vs %g', NF, NF_);
    end
    % Measure the entanglement
    logNB = log2(2*NB+1);
    logNF = log2(2*NF_+1);
 end

 

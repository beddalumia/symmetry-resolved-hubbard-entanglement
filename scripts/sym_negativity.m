%% Compute symmetry-resolved negativity for a two-orbital RDM in usual EDIPACK ordering
%
%   >> [Nud,Nft,Nhd] = sym_negativity(RDMij)
%
%  © Gabriele Bellomia, 2025

function [Nud,Nft,Nhd] = sym_negativity(RDMij)
    
% Build charge-jump blocks RDM_q (q: charge jump on each site)

    % q = 0 -> spin-singlets, fixed charge, fluctuating spin
    RDM_0 = zeros(4,4);
    RDM_0(1,1) = RDMij(4,4);
    RDM_0(2,3) = RDMij(7,10);
    RDM_0(3,2) = RDMij(10,7);
    RDM_0(4,4) = RDMij(13,13);
    
    % q = 1 -> hopping amplitudes
    % > unimplemented (todo...) 
    
    % q = 2 -> holon-doublon binding (related to pairing?)
    RDM_2 = zeros(4,4);
    RDM_2(1,1) = RDMij(1,1);
    RDM_2(2,3) = RDMij(6,11);
    RDM_2(3,2) = RDMij(11,6);
    RDM_2(4,4) = RDMij(16,16);
    
% Build partial transposes (no need to enforce twisted phases under P-SSR)
    
    % q = 0
    RDM_0(1,4) = -RDM_0(2,3); % Fermionic sign in the
    RDM_0(4,1) = -RDM_0(3,2); % partial transpose! :O
    RDM_0(3,2) = 0;
    RDM_0(2,3) = 0;
    
    % q = 1 [here we need to take care of twisted phases]
    % TODO: build this directly instead of using the full partial transpose

    % q = 2
    RDM_2(1,4) = RDM_2(2,3);
    RDM_2(4,1) = RDM_2(3,2);
    RDM_2(3,2) = 0;
    RDM_2(2,3) = 0;
    
    % Evaluate transposed spectra...
    p0 = eig(RDM_0);
    p2 = eig(RDM_2);
    % ...and their negativity :)
    Nud = -sum(p0(p0<0));
    Nft = NaN; % Unimplemented (TODO)
    Nhd = -sum(p2(p2<0));
    
    % > Let's leverage the additivity of negativities to get the q = 1 term!
    [FT, BT] = partial_transpose(RDMij);
    Ntot = 0.5*(trace(sqrtm(FT'*FT))-1);
    Ntot_bosonic = 0.5*(trace(sqrtm(BT'*BT))-1);
    if abs(Ntot_bosonic - Ntot) < 1e-12
       warning("You might want to investigate why bosonic and fermionic PPT coincide");
    end
    Nft = Ntot - Nhd - Nud;
    
end
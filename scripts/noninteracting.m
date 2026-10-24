%% Compute Two-Site RDMs for the noninteracting square lattice, in the large size limit
clear all
L = 100; % Nsites = L^2;
t = 1/4; % hopping amplitude
kspan = (1-mod(L,2)-round((L-1)/2):1:round((L-1)/2)) * 2*pi/L;
for ix = 1:L
    for iy = 1:L
        k{(ix-1)*L+iy} = [kspan(ix),kspan(iy)];
        E((ix-1)*L+iy) = -2*t*(cos(kspan(ix))+cos(kspan(iy)));
        Esurf(ix,iy) = E((ix-1)*L+iy);
    end
end
% imagesc(abs(Esurf)<1e-12)
[~,indices] = sort(E); sorted_k = k(indices);

% pre-compute Slater-Condon rules for 2 orbitals
sc_matrix = slater_condon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]

% Filling of the band η = N/(2*L)
global eta
fillings = [0.5]
Eud = zeros(4,length(fillings));
Ehd = zeros(4,length(fillings));
logNud = zeros(4,length(fillings));
logNft = zeros(4,length(fillings));
logNhd = zeros(4,length(fillings));
logN = zeros(4,length(fillings));
for d=1:length(fillings)
    eta = fillings(d) 
    N = round(2*L^2*eta);
    kx = zeros(round(N/2),2);
    ky = zeros(round(N/2),2);
    for n = 1:round(N/2)
        for s = 1:2
            kx(n,s) = sorted_k{n}(1);
            ky(n,s) = sorted_k{n}(2);
        end
    end
    
    %% Numerical 1-body density matrix: vectorized code
    disp("Computing the ONE-BODY DENSITY MATRIX")

    % 1. Create a 1D vector of all coordinates (replaces the ix, iy, jx, jy loops)
    i_vec = (1:L^2).'; 
    ix_vec = floor((i_vec - 1) / L) + 1;
    iy_vec = mod(i_vec - 1, L) + 1;

    % 2. Extract k-vectors as row vectors (for matrix broadcasting)
    kx_up = kx(:,1).';  
    ky_up = ky(:,1).';  

    kx_dw = kx(:,2).';  
    ky_dw = ky(:,2).';  

    % 3. Compute phase matrices via outer products (Size: L^2 x N_k)
    % This replaces the innermost 'n' loops
    phase_up = ix_vec * kx_up + iy_vec * ky_up;
    phase_dw = ix_vec * kx_dw + iy_vec * ky_dw;

    A_up = exp(1i * phase_up);
    A_dw = exp(1i * phase_dw);

    % 4. Compute the OBDM blocks using matrix multiplication (A * A^dagger)
    % The ' operator automatically applies the complex conjugate transpose
    OBDM_up = (A_up * A_up') / L^2;
    OBDM_dw = (A_dw * A_dw') / L^2;

    % 5. Assemble the final 2x2 block diagonal OBDM matrix
    OBDM = zeros(2*L^2, 2*L^2);
    OBDM(1:L^2, 1:L^2)         = OBDM_up;
    OBDM(L^2+1:end, L^2+1:end) = OBDM_dw;
    
     if L <= 30
        OBDM_vectorized = OBDM;
        %% Numerical 1-body density matrix: loops
        OBDM = zeros(2*L^2,2*L^2); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
        disp("Recomputing the ONE-BODY DENSITY MATRIX to check vectorized code")
        for ix = 1:L
            for iy = 1:L
                i = (ix-1)*L + iy;
                for jx = 1:L
                    for jy = 1:L
                        j = (jx-1)*L + jy;
                        for n = 1:length(kx(:,1))
                            % Spin-up
                            OBDM(i,j) = OBDM(i,j) + exp(1i*((ix-jx)*kx(n,1)+(iy-jy)*ky(n,1))) / L^2;
                        end
                        for n = 1:length(kx(:,2))
                            % Spin-dw
                            OBDM(i+L^2,j+L^2) = OBDM(i+L^2,j+L^2) + exp(1i*((ix-jx)*kx(n,2)+(iy-jy)*ky(n,2))) / L^2;
                        end
                        % TODO: harder to vectorize in 2D, but we gotta do it :(
                    end
                end
            end
        end
        err = norm(OBDM-OBDM_vectorized);
        assert(err<1e-12,"Vectorized code is not working",err)
     end

    %% Analytical benchmark (valid only at half filling)
    OBDM_true = zeros(4,4); % [1,2,...,Nsites]_up[1,2,...,Nsites]_dw
    disp("Building analytical formula in the thermodynamic limit to check size scaling")
    for ix = 1:2
        for iy = 1:1
            i = (ix-1) + iy;
            for jx = 1:2
                for jy = 1:1
                    j = (jx-1) + jy; 
                    if ix==jx && iy==jy
                        OBDM_true(i,j) = 0.5;
                    else
                        dd = (ix-jx)^2-(iy-jy)^2;
                        if dd > 0
                            OBDM_true(i,j) = -((-1)^(ix-jx)-(-1)^(iy-jy))/(pi^2*dd);
                        end
                    end
                    if ix==jx && iy==jy
                        OBDM_true(i+2,j+2) = 0.5;
                    else
                        dd = (ix-jx)^2-(iy-jy)^2;
                        if dd > 0
                            OBDM_true(i+2,j+2) = -((-1)^(ix-jx)-(-1)^(iy-jy))/(pi^2*dd);
                        end
                    end
                end
            end
        end
    end

    %% Build two-orbital RDMs
    OBDM_ij = zeros(4,4);
    counter = 0;
    disp("Computing the THE TWO-SITE DENSITY MATRIX")
    for i = 2%L
        for j = 1%[2,3,11,100]

            % Restrict the one-body density matrix to two orbitals
            % -> spin up
            OBDM_ij(1,1) = OBDM(i,i); 
            OBDM_ij(1,2) = OBDM(i,j); 
            OBDM_ij(2,1) = OBDM(j,i);
            OBDM_ij(2,2) = OBDM(j,j);
            % -> spin down
            OBDM_ij(3,3) = OBDM(i+L^2,i+L^2); 
            OBDM_ij(3,4) = OBDM(i+L^2,j+L^2); 
            OBDM_ij(4,3) = OBDM(j+L^2,i+L^2);
            OBDM_ij(4,4) = OBDM(j+L^2,j+L^2);
            % To check how close to 1 the eigvals are
            %  disp(eig(OBDM_ij))
            %  disp(sum(eig(OBDM_ij)))
            %  disp ********************

            disp("Peschel's formula")
            RDM = RDM0(OBDM_ij,2,sc_matrix);
            RDM_true = RDM0(OBDM_true,2,sc_matrix);
            
            %% Apply SSR and measure entanglement
            counter = counter + 1;
            [~,logN(counter,d)] = negativity(RDM);
            [Nud,Nft,Nhd] = sym_negativity(RDM);
            logNud(counter,d) = log2(2*Nud+1);
            logNft(counter,d) = log2(2*Nft+1);
            logNhd(counter,d) = log2(2*Nhd+1);
            [Eud(counter,d),Ehd(counter,d)] = sym_ree(RDM);

        end
    end

    writematrix(RDM_true,"reduced_density_matrix_2sites_noninteracting.dat","Delimiter","\t")
    writematrix(RDM,"reduced_density_matrix_2sites_noninteracting_size.dat","Delimiter","\t")
    assert(norm(RDM-RDM_true)<1e-3)

    %% Eq. 27: nearest-neighbor dimer in the infinite noninteracting lattice
    RDMi_true = RDM0(OBDM_true([1,3],[1,3]),1,slater_condon(1)); % single site
    S_true = vonNeumann(RDM_true);
    I_true = 2*vonNeumann(RDMi_true) - S_true;
    [Nud,Nft,Nhd] = sym_negativity(RDM_true);
    fprintf("\nNearest neighbors in the infinite noninteracting lattice:\n")
    fprintf("  <cdg_i c_j> = %.4f (2/pi^2 = %.4f)\n",OBDM_true(1,2),2/pi^2)
    fprintf("  N^ud_ij     = %.4f bit\n",log2(2*Nud+1))
    fprintf("  N^hd_ij     = %.4f bit\n",log2(2*Nhd+1))
    fprintf("  N^Ft_ij     = %.4f bit\n",log2(2*Nft+1))
    fprintf("  I_ij        = %.4f bit\n",I_true)
    fprintf("  S_ij        = %.4f bit\n",S_true)

end

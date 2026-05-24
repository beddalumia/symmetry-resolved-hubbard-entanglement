set(0,'defaulttextinterpreter','latex')
addpath ../lib/qcmp-lab/
QcmP.plot.import_colorlab

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../half-filled')

[u,udir] = QcmP.post.get_list('U');

addpath(HERE)

hop_cdmft = zeros(size(u));
Iij_cdmft = zeros(size(u));
Nft_cdmft = zeros(size(u));
Nud_cdmft = zeros(size(u));
Nhd_cdmft = zeros(size(u));
for i=1:length(u)
    cd(udir(i))
    RDM1 = load('reduced_density_matrix_1sites.dat');
    RDM2 = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    cd ..
    fprintf("CDMFT @ U = %g\n\n",u(i))
    disp("Local OBDM")
    OBDM1 = get_obdm(RDM1);
    disp("2-site OBDM")
    OBDM2 = get_obdm(RDM2);
    hop_cdmft(i) = OBDM2(1,2);
    Iij_cdmft(i) = 2*vonNeumann(RDM1)-vonNeumann(RDM2);
    [Nud_cdmft(i),Nft_cdmft(i),Nhd_cdmft(i)] = sym_negativity(RDM2);
end

cd CRISBreference

[uRISB,uRISBdir] = QcmP.post.get_list('U');

hop_crisb = zeros(size(uRISB));
Iij_crisb = zeros(size(u));
Nft_crisb = zeros(size(u));
Nud_crisb = zeros(size(u));
Nhd_crisb= zeros(size(u));
for i=1:length(uRISB)
    cd(uRISBdir(i))
    RDM1 = load('reduced_density_matrix_i1l1.dat');
    RDM2 = load('reduced_density_matrix_i1l1_i2l1.dat');
    cd ..
    fprintf("CRISB U = %g\n\n",uRISB(i))
    disp("Local OBDM")
    OBDM1 = get_obdm(RDM1);
    disp("2-site OBDM")
    OBDM2 = get_obdm(RDM2);
    hop_crisb(i) = OBDM2(1,2);
    Iij_crisb(i) = 2*vonNeumann(RDM1)-vonNeumann(RDM2);
    [Nud_crisb(i),Nft_crisb(i),Nhd_crisb(i)] = sym_negativity(RDM2);
end

cd ..

rmpath(HERE)
cd(HERE)

%% Plotting
tiledlayout(5,1)

nexttile % OBDM hopping term
plot(u,hop_cdmft,'s:','Color',str2rgb('teal'))
hold on
plot(uRISB,hop_crisb,'.:','Color',str2rgb('teal'))
scatter(0, 2/pi^2, 100,'pentagram','filled','MarkerEdgeColor',str2rgb('k'),'MarkerFaceColor',str2rgb('amber'));
xlim([-0.2,5.2])
xticks([])
xline(1.55,'--','Color',[0.6,0.6,0.6])
xline(2.15,'-.','Color',[0.6,0.6,0.6])
ylim([0,0.25])
ytickformat('%.2f')
ylabel("$\langle c^\dagger_{i\sigma} c_{j\sigma} \rangle$")
legend(["CDMFT","CRISB"],"Location","best","box","off","Interpreter","latex")

nexttile % Mutual Information
plot(u,Iij_cdmft,'s:','Color',str2rgb('matlab4'))
hold on
plot(uRISB,Iij_crisb,'.:','Color',str2rgb('matlab4'))
scatter(0, 0.4879, 100,'pentagram','filled','MarkerEdgeColor',str2rgb('k'),'MarkerFaceColor',str2rgb('amber'));
xlim([-0.2,5.2])
xticks([])
xline(1.55,'--','Color',[0.6,0.6,0.6])
xline(2.15,'-.','Color',[0.6,0.6,0.6])
ylim([0.4,1.1])
yticks([0.4,0.6,0.8,1])
ytickformat('%.2f')
ylabel("$I_{\langle ij \rangle}$ [bit]")
legend(["CDMFT","CRISB"],"Location","best","box","off","Interpreter","latex")

nexttile % Hopping negativity
plot(u,log2(2*Nft_cdmft+1),'s:','Color',str2rgb('amber'))
hold on
plot(uRISB,log2(2*Nft_crisb+1),'.:','Color',str2rgb('amber'))
scatter(0, 0.4388, 100,'pentagram','filled','MarkerEdgeColor',str2rgb('k'),'MarkerFaceColor',str2rgb('amber'));
xlim([-0.2,5.2])
xticks([])
xline(1.55,'--','Color',[0.6,0.6,0.6])
xline(2.15,'-.','Color',[0.6,0.6,0.6])
ylim([0.33,0.63])
yticks([0.3,0.4,0.5,0.6])
ytickformat('%.2f')
ylabel("$\mathcal{N}^{\mathrm{F}t}_{\langle ij \rangle}$ [bit]")
legend(["CDMFT","CRISB"],"Location","best","box","off","Interpreter","latex")

nexttile % AFM negativity
plot(u,log2(2*Nud_cdmft+1),'s:','Color',str2rgb('cornflower'))
hold on
plot(uRISB,log2(2*Nud_crisb+1),'.:','Color',str2rgb('cornflower'))
scatter(0, 0, 100,'pentagram','filled','MarkerEdgeColor',str2rgb('k'),'MarkerFaceColor',str2rgb('amber'));
xlim([-0.2,5.2])
xticks([])
xline(1.55,'--','Color',[0.6,0.6,0.6])
xline(2.15,'-.','Color',[0.6,0.6,0.6])
ylim([0,0.6])
yticks([0,0.2,0.4,0.6])
ytickformat('%.2f')
ylabel("$\mathcal{N}^{\uparrow\downarrow}_{\langle ij \rangle}$ [bit]")
legend(["CDMFT","CRISB"],"Location","best","box","off","Interpreter","latex")

nexttile % Holon-Doublon negativity
plot(u,log2(2*Nhd_cdmft+1),'s:','Color',str2rgb('rose pink'))
hold on
plot(uRISB,log2(2*Nhd_crisb+1),'.:','Color',str2rgb('rose pink'))
scatter(0, 0, 100,'pentagram','filled','MarkerEdgeColor',str2rgb('k'),'MarkerFaceColor',str2rgb('amber'));
xlim([-0.2,5.2])
xlabel("$U/D$")
xline(1.55,'--','Color',[0.6,0.6,0.6])
xline(2.15,'-.','Color',[0.6,0.6,0.6])
ylim([0,0.06])
yticks([0,0.02,0.04,0.06])
ytickformat('%.2f')
ylabel("$\mathcal{N}^\mathrm{hd}_{\langle ij \rangle}$ [bit]")
legend(["CDMFT","CRISB"],"Location","best","box","off","Interpreter","latex")

%% contains

function obdm = get_obdm(RDM)
% Compute one-body density matrix via Slater Condon rules
% and the closest free state (Gaussian density matrix). 

    % Get number of single-body modes
    N = size(RDM,1); assert(N==size(RDM,2),"Not a square matrix!");
    nmodes = round(log(N)/log(4)); % N = 4^nmodes

    % Pre-compute all the one-body Slater-Condon excitations
    sc_matrix = slater_condon(nmodes);

    % Compute the /spinful/ one-body density matrix ❬c†_i c_j❭
    obdm = zeros(2*nmodes,2*nmodes);
    for n = 1:2*nmodes
        for m = 1:2*nmodes
            obdm(n,m) = trace(RDM * squeeze(sc_matrix(n,m,:,:)));
        end
    end
    % logging
    disp(obdm)

end


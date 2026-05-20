set(0,'defaulttextinterpreter','latex')
addpath ../lib/qcmp-lab/

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../half-filled')

[u,udir] = QcmP.post.get_list('U');

addpath(HERE)
for i=1:length(u)
    cd(udir(i))
    RDM1 = load('reduced_density_matrix_1sites.dat');
    RDM2 = QcmP.post.get_Hloc('reduced_density_matrix_2sites.dat');
    cd ..
    fprintf(udir(i))
    disp("Local OBDM")
    OBDM1 = get_obdm(RDM1);
    disp("2-site OBDM")
    OBDM2 = get_obdm(RDM2);
end

cd(HERE)

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


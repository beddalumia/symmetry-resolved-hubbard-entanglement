%% Next-nearest-neighbor superselected entanglement (Sec. 4.1)
%  The antiferromagnetic and holon-doublon sectors of the next-nearest-neighbor
%  density matrices, i.e. the plaquette diagonals (1,4) and (2,3), are expected
%  to be separable at every computed point. Here we check it, both at half
%  filling and in the doped regime.

addpath ../lib/qcmp-lab/

HERE = fileparts(mfilename('fullpath'));
addpath(HERE)

tol = 1e-12; % same PPT tolerance as in sym_ree
datasets = {'half-filled','U'; 'doped','xmu'};
pairs = {'i1l1_i4l1','i2l1_i3l1'};

for d = 1:size(datasets,1)

    cd(fullfile(HERE,'..',datasets{d,1}))
    [x,xdir] = QcmP.post.get_list(datasets{d,2});

    Nud = zeros(length(x),length(pairs)); Nft = Nud; Nhd = Nud; Eud = Nud; Ehd = Nud;
    for i = 1:length(x)
        cd(xdir(i))
        for k = 1:length(pairs)
            RDM = rdm_load(['reduced_density_matrix_',pairs{k},'.dat']);
            [Nud(i,k),Nft(i,k),Nhd(i,k)] = sym_negativity(RDM);
            [Eud(i,k),Ehd(i,k)] = sym_ree(RDM);
        end
        cd('..')
    end

    fprintf('\n%s: next-nearest neighbors (1,4) and (2,3), %d points\n',datasets{d,1},length(x))
    fprintf('  max N^ud = %.1e bit,  max N^hd = %.1e bit\n',max(log2(2*Nud(:)+1)),max(log2(2*Nhd(:)+1)))
    fprintf('  max E^ud = %.1e bit,  max E^hd = %.1e bit\n',max(Eud(:)),max(Ehd(:)))
    fprintf('  N^Ft between %.3f and %.3f bit\n',min(log2(2*Nft(:)+1)),max(log2(2*Nft(:)+1)))
    entangled = any(Nud > tol | Nhd > tol, 2);
    if any(entangled)
        fprintf('  superselected entanglement found at %s = %s\n',datasets{d,2},mat2str(x(entangled)'))
    else
        fprintf('  no superselected entanglement at any point\n')
    end

end

rmpath(HERE)
cd(HERE)

%% Utilities

function rdm = rdm_load(fname)
   try
      rdm = QcmP.post.get_Hloc(fname);
   catch
      rdm = load(fname);
   end
end

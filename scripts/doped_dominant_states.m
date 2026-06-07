set(0,'defaulttextinterpreter','latex')
addpath ../lib/qcmp-lab/

%% Data wrangling

[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../doped')

[mu,mudir] = QcmP.post.get_list('xmu');

QcmP.post.observables_line('xmu','site001');
dens = load('dens__site001.txt');

RDM12 = cell(size(mu));
p1 = zeros(size(mu));
p2 = zeros(size(mu));
p3 = zeros(size(mu));
p4 = zeros(size(mu));
p5 = zeros(size(mu));
p6 = zeros(size(mu));
for i = 1:length(mu)
    cd(mudir(i))
    try
        RDM12{i} = QcmP.post.get_Hloc('reduced_density_matrix_i1l1_i2l1.dat');
    catch
        RDM12{i} = load('reduced_density_matrix_i1l1_i2l1.dat');
    end
    probs = eig(RDM12{i});
    probs = sort(probs);
    p1(i) = probs(end);
    p2(i) = probs(end-1);
    p3(i) = probs(end-2);
    p4(i) = probs(end-3);
    p5(i) = probs(end-4);
    p6(i) = probs(end-5);
    cd('..')
end

cd(HERE)


%% Actual graphics

QcmP.plot.import_colorlab

scatter(1-dens(mu>-0.84),p1(mu>-0.84),'v')
hold on
scatter(1-dens(mu>-0.84),p2(mu>-0.84),'+')
scatter(1-dens(mu>-0.84),p3(mu>-0.84),'x')
scatter(1-dens(mu>-0.84),p4(mu>-0.84),'o')
scatter(1-dens(mu>-0.84),p5(mu>-0.84),'s')
scatter(1-dens(mu>-0.84),p6(mu>-0.84),'^')

stop

tiledlayout(1,3)

nexttile

% Symmetry-resolved entanglement
plot(1-dens,En,':+','LineWidth',0.5,'Color',str2rgb('sky')); hold on
plot(1-dens,Ed,':x','LineWidth',0.5,'Color',str2rgb('strawberry'))

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
ytickformat('%.2f')
%ylim([0,1.05]);
legend(["$E^{\uparrow\downarrow}_{\langle ij \rangle}$",...
        "$E^\mathrm{hd}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','north outside','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])

nexttile

% Symmetry-resolved negativity
plot(1-dens,log2(2*N0+1),':+','LineWidth',0.5,'Color',str2rgb('cornflower')); hold on
plot(1-dens,log2(2*N1+1),':o','LineWidth',0.5,'Color',str2rgb('amber'))
plot(1-dens,log2(2*N2+1),':x','LineWidth',0.5,'Color',str2rgb('rose pink'))

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
ytickformat('%.2f')
%ylim([0,1.05]);
legend(["$\mathcal{N}^{\,\uparrow\downarrow}_{\langle ij \rangle}$",...
        "$\mathcal{N}^{\mathrm{F}t}_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{hd}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','north outside','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])

nexttile

% Full entanglement
plot(1-dens,MI,':s','LineWidth',0.5,'Color',str2rgb('matlab4')); hold on
plot(1-dens,logN,':d','LineWidth',0.5,'Color',str2rgb('lilac'));

% Axes
xlabel("$\delta = 1-n$")
ylabel("[bit]")
ytickformat('%.2f')
%ylim([0,1.05]);
legend(["$I_{\langle ij \rangle}$",...
        "$\mathcal{N}^\mathrm{F}_{\langle ij \rangle}$"],...
    "Interpreter",'latex','Location','north outside','Orientation','horizontal')
legend('boxoff')
xlim([-0.02,0.4])
ylim([0,1])

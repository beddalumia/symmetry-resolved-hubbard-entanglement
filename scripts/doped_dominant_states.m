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
p0 = zeros(size(mu));
p1 = zeros(size(mu));
p2 = zeros(size(mu));
p3 = zeros(size(mu));
pu = zeros(size(mu));
pd = zeros(size(mu));
ps = zeros(size(mu));
ph = zeros(size(mu));

% Triplet states in EDIpack basis
t1 = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
t2 = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0];
t3 = [0 0 0 0 0 0 1 0 0 -1 0 0 0 0 0 0]/sqrt(2);

% Singlet state in EDIpack basis
ss = [0 0 0 0 0 0 1 0 0 1 0 0 0 0 0 0]/sqrt(2);

% Holon-doublon binding in EDIpack basis
hd = [0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0]/sqrt(2);

% Holon-spin binding in EDIpack basis
up = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0]/sqrt(2);
dw = [0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0]/sqrt(2);


for i = 1:length(mu)
    cd(mudir(i))
    try
        RDM12{i} = QcmP.post.get_Hloc('reduced_density_matrix_i1l1_i2l1.dat');
    catch
        RDM12{i} = load('reduced_density_matrix_i1l1_i2l1.dat');
    end
    [V, D] = eig(RDM12{i}); probs = diag(D);
    % Dominant state (spin-singlet+holon-doublon)
    p0(i) = probs(end); % it is sorted, of course
    % Decomposition of dominant state in singlet and holon-doublon
    ps(i) = abs(ss * V(:,end))^2; 
    ph(i) = abs(hd * V(:,end))^2;
    % Find triplet states (via overlaps)
    o1 = t1 * V; [~,max1] = max(abs(o1)); p1(i) = probs(max1);
    o2 = t2 * V; [~,max2] = max(abs(o2)); p2(i) = probs(max2);
    o3 = t3 * V; [~,max3] = max(abs(o3)); p3(i) = probs(max3);
    % Find holon-spin states (via overlap)
    ou = up * V; [~,maxu] = max(abs(ou)); pu(i) = probs(maxu);
    od = dw * V; [~,maxd] = max(abs(od)); pd(i) = probs(maxd);

    cd('..')
end

cd(HERE)


%% Actual graphics

QcmP.plot.import_colorlab

plot(1-dens(mu>-0.84&dens<1),p0(mu>-0.84&dens<1),'d:','Color',str2rgb('cornflower'),'MarkerFaceColor',str2rgb('light periwinkle'))
hold on
plot(1-dens(mu>-0.84&dens<1),p3(mu>-0.84&dens<1),'s:','Color',str2rgb('hospital green'),'MarkerFaceColor',str2rgb('pale green'))
plot(1-dens(mu>-0.84&dens<1),p1(mu>-0.84&dens<1),'x:','Color',str2rgb('forest'))
plot(1-dens(mu>-0.84&dens<1),p2(mu>-0.84&dens<1),'+:','Color',str2rgb('teal'))
plot(1-dens(mu>-0.84&dens<1),pu(mu>-0.84&dens<1),'^:','Color',str2rgb('amber'),'MarkerFaceColor',str2rgb('lemon'))
plot(1-dens(mu>-0.84&dens<1),pd(mu>-0.84&dens<1),'v:','Color',str2rgb('amber'),'MarkerFaceColor',str2rgb('lemon'))
legend(["$C^{\uparrow\downarrow}\frac{|\uparrow\rangle\otimes|\downarrow\rangle-|\downarrow\rangle\otimes|\uparrow\rangle}{\sqrt{2}} \enspace+\enspace C^\mathrm{hd}\frac{|\bullet\rangle\otimes|\uparrow\downarrow\rangle+|\uparrow\downarrow\rangle\otimes|\bullet\rangle}{\sqrt{2}}$",...
    "$\frac{|\uparrow\rangle\otimes|\downarrow\rangle+|\downarrow\rangle\otimes|\uparrow\rangle}{\sqrt{2}}$",...
    "$|\!\!\uparrow\rangle\otimes|\!\!\uparrow\rangle$",...
    "$|\!\!\downarrow\rangle\otimes|\!\!\downarrow\rangle$",...
    "$\frac{|\bullet\rangle\otimes|\uparrow\rangle+|\uparrow\rangle\otimes|\bullet\rangle}{\sqrt{2}}$",...
    "$\frac{|\bullet\rangle\otimes|\downarrow\rangle+|\downarrow\rangle\otimes|\bullet\rangle}{\sqrt{2}}$"],...
    'Interpreter','latex',...
    'Location','eastoutside',...
    'box','off');
xlabel("$\delta = 1-n$")
xticks([0,0.1,0.2,0.3])
yticks([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7])

% Create an inset
% Position vector format: [left, bottom, width, height] (from 0 to 1)
inset_ax = axes('Position', [0.245  0.67  0.19  0.24]);
% Plot singlet and hd weights (dominant state)
plot(1-dens(mu>-0.84&dens<1),ps(mu>-0.84&dens<1),'o:','Color',str2rgb('sky'),'MarkerFaceColor',str2rgb('light sky blue'))
hold(inset_ax, 'on')
plot(1-dens(mu>-0.84&dens<1),ph(mu>-0.84&dens<1),'s:','Color',str2rgb('strawberry'),'MarkerFaceColor',str2rgb('pale pink'))
set(inset_ax,'Ylim',[0.01,5])
set(inset_ax,'Yscale','log')
xticks([0,0.1,0.2,0.3])
xticklabels([]);
legend(["$\bigl|C^{\uparrow\downarrow}\bigr|^2$","$\bigl|C^\mathrm{hd}\bigr|^2$"],...
    'Interpreter','latex',...
    'Location', 'northoutside',...
    'Orientation','horizontal',...
    'box','off')
%xlabel("$\delta = 1-n$")
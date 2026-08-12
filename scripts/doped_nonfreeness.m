set(0,'defaulttextinterpreter','latex')
addpath ../lib/qcmp-lab/

%% LIBRARIES
QcmP.plot.import_colorlab

%% Dirty path selector
[filepath,name,extension] = fileparts(mfilename('fullpath'));
HERE = erase(filepath,name);

cd('../doped')

[mu,mudir] = QcmP.post.get_list('xmu')

si = QcmP.post.eentropy_line('xmu','1sites');
su = zeros(size(si));
sd = zeros(size(si));

nup = load('nup_1_site001.txt');
ndw = load('ndw_1_site001.txt');
nud = nup + ndw;

for i=1:length(mu)
   RDMu=[(1-nup(i)),0;0,nup(i)]; su(i) = vonNeumann(RDMu);
   RDMd=[(1-ndw(i)),0;0,ndw(i)]; sd(i) = vonNeumann(RDMd);
end

intra = su+sd-si;


figure("Name",'Doped Local Nonfreeness')
hold on
% Data plotting
plot(mu,si,'-s','Color',str2rgb('mauve'),...
   'MarkerSize',5,'MarkerEdgeColor',str2rgb('mauve'),...
   'MarkerFaceColor',str2rgb('light rose'),'Linewidth',0.5);
plot(mu,nud,'-o','Color',str2rgb('dark sky blue'),...
   'MarkerSize',4.5,'MarkerEdgeColor',str2rgb('dark sky blue'),...
   'MarkerFaceColor',str2rgb('pale sky blue'),'Linewidth',0.5)
plot(mu,intra,'-d','Color',str2rgb('grass'),...
   'MarkerSize',3.7,'MarkerEdgeColor',str2rgb('grass'),...
   'MarkerFaceColor',str2rgb('yellowish green'),'Linewidth',0.5);
xlim([-1.01,-0.35]);ylim([0,2.25]); box on;
xlabel('$\mu/D$','Interpreter','latex');
xtickformat('%.1f')
ytickformat('%.2f')
% Phase separators
annotation('textbox',[0.15,0.77,0.1,0.05],'string','Fermi liquid',...
   'Interpreter','latex','EdgeColor','none')
xline(-0.715,':')
annotation('textbox',[0.48,0.77,0.1,0.05],'string','pseudogap metal',...
   'Interpreter','latex','EdgeColor','none')
xline(-0.485,':')
annotation('textbox',[0.75,0.77,0.1,0.05],'string','Mott insulator',...
   'Interpreter','latex','EdgeColor','none')
% Outside legend
legend(["$S_i$ [bit]",...
        "$\langle n_i \rangle$",...
        "$I(\,\uparrow\,:\,\downarrow\,)$ [bit]"],...
        "Interpreter","latex","box","off","orientation","horizontal",...
        "location","northoutside")


%close all
cd(HERE);

%% Utilities

function E = vonNeumann(RDM)
   p = eig(RDM);
   E = -sum(p.*log2(p))
end
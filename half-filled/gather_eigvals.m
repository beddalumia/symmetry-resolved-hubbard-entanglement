
[U,UDIRs] = QcmP.post.get_list('U');
d1 = zeros(size(U));
d2 = zeros(size(U));
d3 = zeros(size(U));
d4 = zeros(size(U));

for i = 1:length(U)
    cd(UDIRs(i))
    tab = readtable('eigenvalues_list.ed','FileType','text');
    E = tab.x_X; E = sort(E(E<0));
    for j = 1:length(E)
        x(j) = U(i);
        y(j) = E(j);
    end
    scatter(x,y); hold on
    d1(i) = E(2)-E(1);
    d2(i) = E(3)-E(1);
    d3(i) = E(4)-E(1);
    d4(i) = E(5)-E(1);
    cd ..
end

% plot(4*U,d1); hold on
% plot(4*U,d2,'--');
% plot(4*U,d3,'-.');
% plot(4*U,d4,':');

xlabel('$U/t$','Interpreter','latex');
% ylabel('$E_n - E_0$','Interpreter','latex');
ylabel('All eigenvalues','Interpreter','latex');
title('$2\times2$ plaquette - $N_\mathrm{bath}=8$','Interpreter','latex');
% legend(["$n=1$","$n=2$","$n=3$","$n=4$"],'Interpreter','latex');


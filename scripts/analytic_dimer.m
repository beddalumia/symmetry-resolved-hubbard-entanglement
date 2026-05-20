set(0,'defaulttextinterpreter','latex')
addpath ../lib/qcmp-lab/
QcmP.plot.import_colorlab

Uspan = 0:0.1:50;
t = 1;

global N
global E

% Pre-compute Slater-Condon rules for 2 orbitals
sc_matrix = slater_condon(2); % 4D array [2*nmodes,2*nmodes,4^nmodes,4^nmodes]


for i = 1:length(Uspan)

    U = Uspan(i)
    N = ( 2 + 0.5*( ( U+sqrt(16*t^2+U^2) ) / 2 )^2 / t^2 )^(-0.5);
    E = 0.5 * ( U + sqrt(16*t^2+U^2) );

    GSij = groundstate(t,U);
    DMij = GSij'*GSij;
    RDMi = diag([N^2,0.5-N^2,0.5-N^2,N^2]);
    % Measures of full entanglement
    I(i) = 2*vonNeumann(RDMi) - vonNeumann(DMij);
    [B(i),F(i)] = negativity(DMij);
    % Symmetry-resolved entanglement
    [En(i),Ed(i)] = build_symE(DMij);
    p1 = RDMi(1,1); p2 = RDMi(2,2); p3 = RDMi(3,3); p4 = RDMi(4,4); 
    % Double-check with formulas for P-SSR and N-SSR local EE
    % As taken from J. Chem. Theory Comput. 2021, 17, 1, 79–95
    pEE(i) = (p1+p4).*log2(p1+p4) + (p2+p3).*log2(p2+p3) + I(i)/2;
    nEE(i) = (p2+p3).*log2(p2+p3) - p2.*log2(p2) - p3.*log2(p3);
    assert(abs(nEE(i)-En(i))<1e-12)
    assert(abs(pEE(i)-En(i)-Ed(i))<1e-12)
    %e1(i) = -p1.*log2(p1); e2(i) = -p2.*log2(p2); e3(i) = -p3.*log2(p3); e4(i) = -p4.*log2(p4);
    [E0(i),E1(i),E2(i)] = sym_negativity(DMij);

    % One-body density matrix
    obdm = zeros(4,4);
    obdm_ = zeros(4,4);
    for n = 1:4
        for m = 1:4
            obdm(n,m) = trace(DMij * squeeze(sc_matrix(n,m,:,:)));
            % Cross-check
            for istate = 1:16
                for jstate = 1:16
                    obdm_(n,m) = obdm_(n,m) + GSij(istate)*GSij(jstate)*sc_matrix(n,m,istate,jstate);
                end
            end
        end
    end
    

    % Dimer nonfreeness
    N2(i) = vonNeumann(obdm) + vonNeumann(eye(4)-obdm) - vonNeumann(DMij);

end

% Special case: atomic limit
% > the singlet (S1) is degenerate with the triplet states (T1,T2,T3) for t=0

S1 = ([0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0] + ...
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0] ) / sqrt(2); % Be very careful with fermi signs!

T1 = ([0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0] - ...
        [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0] ) / sqrt(2); % Be very careful with fermi signs!
T2 = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0]; % |up up>
T3 = [0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0]; % |dw dw>

DMij = 0.25 * ( S1'*S1 + T1'*T1 + T2'*T2 + T3'*T3);

RDMi = diag([0,0.5,0.5,0]);
% Proof of no correlation in the t=0 limit
I_atomic = 2*vonNeumann(RDMi) - vonNeumann(DMij);
assert(not(I_atomic))

figure
plot(Uspan,I./2,'LineWidth',3,'Color',str2rgb('matlab4'));
hold on
plot(Uspan,En,'Color',str2rgb('sky blue'),'LineWidth',3);
plot(Uspan,I./2-En-Ed,'Color',str2rgb('goldenrod'),'LineWidth',3);
plot(Uspan,Ed,'Color',str2rgb('strawberry'),'LineWidth',3);
plot(Uspan,F,'--','Color',str2rgb('lilac'),'LineWidth',2);
plot(Uspan,log2(2*E0+1),'--','Color',str2rgb('light teal'),'LineWidth',2)
plot(Uspan,log2(2*E1+1),'--','Color',str2rgb('wheat'),'LineWidth',2)
plot(Uspan,log2(2*E2+1),'--','Color',str2rgb('pink'),'LineWidth',2)
%plot(Uspan,En+Ed,'Color',str2rgb('purple'));
legend([ ...
    "full entanglement", ...
    "antiferromagnetic entanglement",...
    "hopping entanglement",...
    "holon-doublon entanglement",...
    "full fermionic negativity", ...
    "antiferromagnetic negativity",...
    "hopping negativity",...
    "holon-doublon negativity"],"Interpreter","latex","box","off")
xlabel('$U/t$','Interpreter','latex')
ylabel('[bit]')



function gs = groundstate(t,U)

    global N
    global E

    gs = N * (  [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0]+...
                [0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0]+...
        0.5*E/t*[0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0]+...
        0.5*E/t*[0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0]); % Be very careful with fermi signs!

end

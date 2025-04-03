%Linearização CSTR
%Mat(h)euses
%% Definição das variáveis
%{
A = 0; %Área para a troca de calor 
CA = 0; %Concentração de área no reator
CAf = 0; %Concentração de área no fluxo de alimentação
Cp = 0; %Capacidade de calor em energia pora massa*temperatura
F = 0; %Fluxo volumétrico
k0 = 0; %fator pré exponencial
R = 0; %Constante ideal de gás
r = 0; %Taxa da reação por unidade de volume
t = 0; %Tempo
T = 0; %Temperatura do reator
Tf = 0; %Temperatura de alimentação
Tj = 0; %Temperatura da jaqueta
Tref = 0; %Temperatura de referência
U = 1/(t*A*T); %Coeficiente geral de transferência de calor energia/(tempo*area*temperatura)
V = 0; %Volume do reator
deltaE = 0; %Energia de ativação
deltaH = 0; %Calor da reação (Negativo)
p = 0; %densidade
%}
syms A CA CAf Cp F k0 R r t T Tf Tj Tref U V deltaE deltaH p;
%% Balança geral de materiais

%dVp = 0;
%Fout = F;
%Fin = F;
%dV = 0;
%dCa = (F*CAf - F*CA - r*V)/V; %Assumindo o valor constante do volume no reator.

%% Balanço de energia
%dT = (F*p*Cp*(Tf-T) + (-deltaH)*V*r-U*A*(T-Tj))/(V*p*Cp);
%Assumindo o constante: volume, capacidade de calor e densidade.

%% Forma das variáveis de estado das equações dinâmicas
%Podemos então escrever que a média por unidade de volume é dada por:
r = k0*exp(-deltaE/(R*T))*CA;

dCa = F/V*(CAf-CA)-r;
dT = F/V*(Tf-T)+(-deltaH/(p*Cp)*r)-U*A/(V*p*Cp)*(T-Tj);

%% Matrizes para linearização Fornecidas no Bequette
f1 = dCa;
f2 = dT;
%definição dos valores:
%Caso 2:
A11 = diff(f1, CA); %o que seria referenciado ao CAs que é atribuído no bequette, onde o subscrito s denota o valor estacionário, não adotado aq.
A12 = diff(f1, T);
A21 = diff(f2, CA);
A22 = diff(f2, T);

V = 1;  % m^3
deltaE = 11843;  % J/mol
pCp = 500;  % J/Km^3
R = 1.987;  % L/molK
UA = 150;  % J/Kh
k0 = 9703 * 3600;  % h^-1
deltaH = -5960;  % J/mol

% Ponto de operação
F = 1;  % m^3/h
CAf = 10;  % kgmol/m^3
Tf = 298;  % K
Tj = 298;  % K
T = 311.2;  % K
CA = 8.564;  % kgmol/m^3

%Podemos fazer com a função diff do matlab, ou simplesmente escrever o resultado como se segue
%Aqui adotaremos o subscrito s
ks = k0*exp(-deltaE/(R*T));
ks_linha = k0*exp(-deltaE/(R*T))*deltaE/(R*T^2);
%Somente então
ks_linha = ks*(deltaE/(R*T^2));

A11 = -F/V - ks;
A12 = -CA*ks_linha;
A21 = (-deltaH)*ks/(pCp);
A22 = -F/V - UA/(V*pCp) + (-deltaH)*CA*ks_linha/(pCp);
A = [A11 A12; A21 A22];

B11 = (CA*F - CA)/V;
B21 = (Tf-T)/V;
B12 = F/V;
B22 = 0;
B13 = 0;
B23 = F/V;
B14 = 0;
B24 = UA/(V*pCp);

B = [B11 B12 B13 B14; B21 B22 B23 B24];

C = [1 0];
D = [0 0 0 0];

G = ss(A,B,C,D);
[num_u, den_u] = ss2tf(A,B,C,D,1);
g_u = tf(num_u, den_u);

Polos = roots(den_u);
tempo_constante = -1./Polos;

h1 = 0.2;
Gh1 = c2d(G, h1, 'zoh');

h2 = 0.5;
Gh2 = c2d(G, h2, 'zoh');

h3 = 1.5;
Gh3 = c2d(G, h3, 'zoh');

dF = 0.01;
dCAf = 0;
dTf = 0;
dTj = 0; 


%% Identificação
%Mesmo procedimento do laboratório de digital

% Entradas utilizadas
du1 = 0.01; du2 = 0; du3 = 0; du4 = 0;

% Realizar simulação não-linear
sim('modeloNL_2022a.slx');
simulacao_nao_linear = ans.simNL;
t = simulacao_nao_linear.time;
Ca_degrau_F = simulacao_nao_linear.signals.values(:,1);
F_degrau = simulacao_nao_linear.signals.values(:,3);

% Estimar parâmetros com dados do modelo simulado
[K11, T11, L11] = modeloDegrau(Ca_degrau_F - Ca_degrau_F(1), F_degrau-F, t(end)-t(end-1));
G11 = tf(K11, [T11 1], 'iodelay', L11);
CaF_Simulado = lsim(G11, F_degrau - F, t);
emqG11 = mean((Ca_degrau_F - (CaF_Simulado+Ca_degrau_F(1))).^2);

% Realizar simulação do modelo discreto e identificação
sim('modeloDiscreto_2022a.slx')
simData = out.simDiscreto;
output_02s_Ca = simData.signals.values(:,2);  % Ta = 0.2s
time = simData.time;
[K11d, T11d, L11d] = modeloDegrau(output_02s_Ca - output_02s_Ca(1), F_degrau-F, time(end)-time(end-1));
G11d = tf(K11d, [T11d 1], 'iodelay', L11d);
CaFd_Simulado = lsim(G11d, F_degrau - F, time);

figure(2)
plot(t, Ca_degrau_F,'--', t, CaF_Simulado+Ca_degrau_F(1),t,output_02s_Ca, t, CaFd_Simulado+Ca_degrau_F(1) );
title("Resposta ao degrau para a saída Ca e a entrada F e a identificação realizada")
xlabel("t (hora)");
ylabel("Ca (kgmol/m³");
legend('Modelo não linear', 'Modelo G11 identificado', 'Modelo linearizado discretizado','Modelo G11 identificado discreto', 'Location', 'east');
grid on;


%% Repetindo o processo para G12
du1 = 0; du2 = 0.007; du3 = 0; du4 = 0;
sim('modeloNL_2022a.slx');
simulacao_nao_linear = ans.simNL;
t = simulacao_nao_linear.time;
Cb_degrau_F = simulacao_nao_linear.signals.values(:,1);
F_degrau = simulacao_nao_linear.signals.values(:,4); % Valor degrau CAf

%{
figure(3)
plot(t,Ca_degrau_F)
title("Resposta ao degrau para a saída Ca e a entrada CAf, tendo variação de 0.007 kgmol por metro cúbico")
xlabel("t (hora)");
ylabel("Ca (kgmol/m³");
grid on;
%}
[K12, T12, L12] = modeloDegrau(Cb_degrau_F - Cb_degrau_F(1), F_degrau-CAf, t(end)-t(end-1));
G12 = tf(K12, [T12 1], 'iodelay', L12);
Cb_F_simulado = lsim(G12, F_degrau - CAf, t);
emqG12 = mean((Cb_degrau_F - (Cb_F_simulado+Cb_degrau_F(1))).^2);


% Realizar simulação do modelo discreto e identificação
sim('modeloDiscreto_2022a.slx');
simData = out.simDiscreto;
output_02s_Ca = simData.signals.values(:,2);  % Ta = 0.2s
time = simData.time;
[K12d, T12d, L12d] = modeloDegrau(output_02s_Ca - output_02s_Ca(1), F_degrau-CAf, time(end)-time(end-1));
G12d = tf(K12d, [T12d 1], 'iodelay', L12d);
CaFd_Simulado = lsim(G12d, F_degrau - CAf, time);

%%

figure(4)
plot(t, Cb_degrau_F,'--', t, Cb_F_simulado+Cb_degrau_F(1), t, output_02s_Ca, t, CaFd_Simulado+Cb_degrau_F(1));
title("Resposta ao degrau para a saída Ca e a entrada CAf e a identificação realizada")
xlabel("t (hora)");
ylabel("Ca (kgmol/m³");
legend('Modelo não linear', 'Modelo G12 identificado','Modelo linearizado discretizado','Modelo G12 identificado discreto', 'Location', 'east');
grid on;

%% Resposta ao degrau para a saída Ca e a entrada T, tendo variação de 1 K

%G13
du1 = 0; du2 = 0; du3 = 1; du4 = 0;
sim('modeloNL_2022a.slx');
simulacao_nao_linear = ans.simNL;
t = simulacao_nao_linear.time;
Ca_degrau_Caf = simulacao_nao_linear.signals.values(:,1);
CaF_degrau = simulacao_nao_linear.signals.values(:,5);

[K13, T13, L13] = modeloDegrau(Ca_degrau_Caf - Ca_degrau_Caf (1), CaF_degrau - Tf, t(end)-t(end-1));
G13 = tf(K13, [T13 1], 'iodelay', L13);
CaF_Simulado = lsim(G13, CaF_degrau - Tf, t);
emqG13 = mean((Ca_degrau_Caf - (CaF_Simulado+Ca_degrau_Caf (1))).^2);

% Realizar simulação do modelo discreto e identificação
sim('modeloDiscreto_2022a.slx');
simData = out.simDiscreto;
output_02s_Ca = simData.signals.values(:,2);  % Ta = 0.2s
time = simData.time;
[K13d, T13d, L13d] = modeloDegrau(output_02s_Ca - output_02s_Ca(1), CaF_degrau - Tf, time(end)-time(end-1));
G13d = tf(K13d, [T13d 1], 'iodelay', L13d);
CaFd_Simulado = lsim(G13d, CaF_degrau - Tf, time);



%%
figure(6)
plot(t,Ca_degrau_Caf , '--', t,CaF_Simulado+Ca_degrau_Caf (1), t, output_02s_Ca, t, CaFd_Simulado+output_02s_Ca(1))

title('Resposta ao degrau para a saída Ca e a entrada Tf e a identificação realizada');
grid on;
legend('Modelo não linear', 'Modelo G13 identificado','Modelo linearizado discretizado','Modelo G13 identificado discreto', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');

%% G14 Resposta ao degrau para saída Ca e entrada Tj com variação de 1 K

du1 = 0; du2 = 0; du3 = 0; du4 = 1;
sim('modeloNL_2022a.slx');
simNL = ans.simNL;
t = simNL.time;
Cb_degrau_Caf = simNL.signals.values(:,1);
CaF_degrau = simNL.signals.values(:,6);

[K14, T14, L14] = modeloDegrau(Cb_degrau_Caf - Cb_degrau_Caf(1), CaF_degrau - Tj, t(end)-t(end-1));
G14 = tf(K14, [T14 1], 'iodelay', L14);
Cb_F_simulado = lsim(G14, CaF_degrau - Tj, t);
emqG14 = mean((Cb_degrau_Caf - (Cb_F_simulado+Cb_degrau_Caf(1))).^2);

% Realizar simluação do modelo discreto e identificação
sim('modeloDiscreto_2022a.slx');
simData = out.simDiscreto;
output_02s_Ca = simData.signals.values(:,2);  % Ta = 0.2s
time = simData.time;
[K14d, T14d, L14d] = modeloDegrau(output_02s_Ca - output_02s_Ca(1), CaF_degrau - Tj, time(end)-time(end-1));
G14d = tf(K14d, [T14d 1], 'iodelay', L14d);
CaFd_Simulado = lsim(G14d, CaF_degrau - Tj, time);


%% 

figure(8)
plot(t,Cb_degrau_Caf, '--', t,Cb_F_simulado+Cb_degrau_Caf(1), t, output_02s_Ca, t, CaFd_Simulado+output_02s_Ca(1))
title('Resposta ao degrau para saída Ca e a entrada Tj e a identificação realizada');
grid on;
legend('Modelo não linear', 'Modelo G14 identificado','Modelo linearizado discretizado','Modelo G14 identificado discreto', 'Location', 'Best');
xlabel('t (h)');
ylabel('Ca (kgmol/m^3)');


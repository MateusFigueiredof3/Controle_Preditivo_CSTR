%% Processo
% I) Desenvolver simulador/sistema
%   - a) Modelo
%   - b) Linearização*
%   - c) Discretização do modelo linear**
%   - d) Simular (sistema linearizado discreto)
% II) Resposta ao degrau + Identificação do modelo contínuo*
% III) DMC sem restrições
% IV) Identificar modelo discretizado**
% V) Comparações dos itens com * e ** (devem ser equivalentes)

%% Processo - Simulação do STH

%Parametros

clc, clear, close all;

V = 283.168e-3;
Vj = 28.319e-3;
rhoCp = 2.165e3;
rhojCpj = 2.165e3;
UA = 3.065;

Fs = 4.72e-4;
Fjs = 7.079e-4;
Tis = 10; %ºC
Tjis = 93.33; %ºC
Ts = 51.67; %ºC
Tjs = 65.56; %ºC


%% Linearização - Modelo Espaço de Estados

A11 = -Fs/V - UA/V/rhoCp;
A12 = UA/V/rhoCp;
A21 = UA/Vj/rhojCpj;
A22 = -Fjs/Vj - UA/Vj/rhojCpj;

A = [A11, A12; A21, A22];

B11 = 0; B12 = (Tis - Ts)/V; B13 = Fs/V; B14 = 0;
B21 = (Tjis - Tjs)/Vj; B22 = 0; B23 = 0; B24 = Fjs/Vj;

Bu = [B11; B21]; % Em relação a entrada Fj
Bd = [B12, B13, B14; B22, B23, B24]; % Em relação as perturbações F, Ti e Tji

B = [Bu, Bd]; % Modelo em espaço de estados agrupando entradas e perturbações
C = [1 0];
D = [0 0 0 0];

G = ss(A,B,C,D);

%% Linearização - Modelo TF
[numu, denu] = ss2tf(A, B, C, D, 1);
gu = tf(numu,denu)
[numd1, dend1] = ss2tf(A, B, C, D, 2);
gd1 = tf(numd1,dend1)
[numd2, dend2] = ss2tf(A, B, C, D, 3);
gd2 = tf(numd2,dend2)
[numd3, dend3] = ss2tf(A, B, C, D, 4);
gd3 = tf(numd3,dend3)

%% Discretização do modelo linear
% Modelo Espaço de Estados discreto com Ta = 50, 150, 300s

h1 = 50;
Gh1 = c2d(G, h1, 'zoh');

h2 = 150;
Gh2 = c2d(G, h2, 'zoh');

h3 = 300;
Gh3 = c2d(G, h3, 'zoh');


%% Discretização do modelo linear - Simulação
% Com base na maior constante de tempo do processo, o tempo de amostragem deve ser Ta <= 313/2 ~= 150

dFj = 4.72e-5; dF = 0; dTi = 0; dTij = 0; 

open('modeloDiscreto.slx')
sim('modeloDiscreto.slx')

simDiscreto = ans.simDiscreto;
t = simDiscreto.time;

Tcont = simDiscreto.signals.values(:,1);
Tdisc50 = simDiscreto.signals.values(:,2);
Tdisc150 = simDiscreto.signals.values(:,3);
Tdisc300 = simDiscreto.signals.values(:,4);

plot(t,Tcont,'--', t,Tdisc50, t,Tdisc150, t,Tdisc300);
grid on;
xlabel('t (s)');
title('Saídas do sistema contínuo e discretizados com Ta = 50, 150 e 300 s');
legend('Contínuo', 'Ta = 50 s', 'Ta = 150 s', 'Ta = 300 s', 'Location', 'Best');
ylabel('T (°C)');


%% Simulação dos sistemas linear e não linear com variaçãona entrada
dFj = 4.72e-5; dF = 0; dTi = 0; dTij = 0; 

open('modeloNL.slx');
sim('modeloNL.slx');
simNL = ans.simNL;
t = simNL.time;
y_NL = simNL.signals.values(:,1);

open('modeloLinear.slx');
outLinear = sim('modeloLinear.slx');
t_lin = outLinear.simLinear.time;
y_lin = outLinear.simLinear.signals.values;

open('modelotf.slx');
dadosTf = sim('modelotf.slx');
t_tf = dadosTf.simTf.time;
y_tf = dadosTf.simTf.signals.values;

plot(t, y_NL, t_lin, y_lin, t_tf, y_tf)
title('Saídas do sistema não linear e linear');
legend('Não linear', 'Linear(SS)','Linear (TF)', 'Location', 'Best');
grid on;
xlabel('t (s)');
ylabel('T (°C)');
%% Simulação dos sistemas linear com variação 10x maior na entrada 1
% Ao aumentar a variação na entrada por 10x, é esperado que a variação na
% saída aumente em 10x no sistema linear. Com a maior variação no ponto de
% operação, o erro entre as saídas do sistema não linear e linear será
% maior
dFj = 4.72e-4; dF = 0; dTi = 0; dTij = 0; 

sim('modeloNL.slx');
simNL = ans.simNL;
t = simNL.time;
y_NL = simNL.signals.values(:,1);

outLinear = sim('modeloLinear.slx');
t_lin = outLinear.simLinear.time;
y_lin = outLinear.simLinear.signals.values;


dadosTf = sim('modelotf.slx');
t_tf = dadosTf.simTf.time;
y_tf = dadosTf.simTf.signals.values;

plot(t, y_NL, t, y_lin, t_tf, y_tf)
title('Saídas do sistema não linear e linear');
legend('Não linear', 'Linear(SS)','Linear (TF)', 'Location', 'Best');
grid on;
xlabel('t (s)');
ylabel('T (°C)');

%% Resposta ao degrau + Identificação do modelo contínuo - G11
% Teste do degrau aplicado na vazão do fluido da jaqueta (Fj) e identificação de G11 = K11/(T11s+1)exp(-sL11)

dFj = 4.72e-5; dF = 0; dTi = 0; dTij = 0; 
h = dFj;

%Usando modeo Linear para identificação
open('modeloLinear.slx');
sim('modeloLinear.slx');
simLinear = ans.simLinear;

t = simLinear.time;
TstepFj = simLinear.signals.values(:,1);
u = h*ones(size(t));

% Identificação G11
[G011, T11, L11] = parametrosFOPTD(TstepFj - TstepFj(1), h, t(end)-t(end-1)) %t = 5
%Usar TstepFj - TstepFj(1) para retirar no offset

G11 = tf(G011, [T11 1], 'iodelay', L11);
TFjsim = lsim(G11, u, t);
emqG11 = mean((TstepFj - (TFjsim+TstepFj(1))).^2)

figure
plot(t,TstepFj, '--', t,TFjsim+TstepFj(1), '-')
title('Comparação entre o teste do degrau e simulação de G11 identificado');
grid on;
xlabel('t (s)');
legend('Modelo linear', 'Modelo G11 identificado', 'Location', 'Best');
ylabel('T (°C)')

%% Resposta ao degrau + Identificação do modelo contínuo - G12
dFj = 0; dF = 4.72e-5; dTi = 0; dTij = 0; 
h = dF;

%Usando modeo NÃO Linear para identificação
sim('modeloNL.slx');
simNL = ans.simNL;

t = simNL.time;
TstepFi = simNL.signals.values(:,1);
Fistep = simNL.signals.values(:,3);
u = h*ones(size(t));

% Identificação G11
[G012, T12, L12] = parametrosFOPTD(TstepFi - TstepFi(1), h, t(end)-t(end-1)) %t = 5

G12 = tf(G012, [T12 1], 'iodelay', L11);
TFisim = lsim(G12, u, t);
emqG12 = mean((TstepFi - (TFisim+TstepFi(1))).^2)

figure
plot(t,TstepFi, '--', t,TFisim+TstepFi(1), '-')
title('Comparação entre o teste do degrau e simulação de G12 identificado');
grid on;
xlabel('t (s)');
legend('Modelo não linear', 'Modelo G12 identificado', 'Location', 'Best');
ylabel('T (°C)')

%% Resposta ao degrau + Identificação do modelo contínuo - G13
dFj = 0; dF = 0; dTi = 1; dTij = 0; 
h = dTi;

%Usando modeo NÃO Linear para identificação
sim('modeloNL.slx');
simNL = ans.simNL;

t = simNL.time;
TstepTi = simNL.signals.values(:,1);
Tistep = simNL.signals.values(:,4);
u = h*ones(size(t));

% Identificação G11
[G013, T13, L13] = parametrosFOPTD(TstepTi - TstepTi(1), h, t(end)-t(end-1)) %t = 5

G13 = tf(G013, [T13 1], 'iodelay', L11);
TTisim = lsim(G13, u, t);
emqG13 = mean((TstepTi - (TTisim+TstepTi(1))).^2)

figure
plot(t,TstepTi, '--', t,TTisim+TstepTi(1), '-')
title('Comparação entre o teste do degrau e simulação de G13 identificado');
grid on;
xlabel('t (s)');
legend('Modelo não linear', 'Modelo G13 identificado', 'Location', 'Best');
ylabel('T (°C)')

%% Resposta ao degrau + Identificação do modelo contínuo - G14

dFj = 0; dF = 0; dTi = 0; dTij = 1; 
h = dTij;

%Usando modeo NÃO Linear para identificação
sim('modeloNL.slx');
simNL = ans.simNL;

t = simNL.time;
TstepTji = simNL.signals.values(:,1);
Tjistep = simNL.signals.values(:,5);
u = h*ones(size(t));

% Identificação G11
[G014, T13, L13] = parametrosFOPTD(TstepTji - TstepTji(1), h, t(end)-t(end-1)) %t = 5

G14 = tf(G014, [T13 1], 'iodelay', L11);
TTjisim = lsim(G14, u, t);
emqG14 = mean((TstepTji - (TTjisim+TstepTji(1))).^2)

figure
plot(t,TstepTji, '--', t,TTjisim+TstepTji(1), '-')
title('Comparação entre o teste do degrau e simulação de G14 identificado');
grid on;
xlabel('t (s)');
legend('Modelo não linear', 'Modelo G14 identificado', 'Location', 'Best');
ylabel('T (°C)')


%% DMC sem restrições


%% Identificação modelo discretizado

dFj = 4.72e-5; dF = 0; dTi = 0; dTij = 0; 
h = dFj;

sim('modeloDiscreto.slx')

simDiscreto = ans.simDiscreto;
t = simDiscreto.time;

Tcont = simDiscreto.signals.values(:,1);
Tdisc50 = simDiscreto.signals.values(:,2);
Tdisc150 = simDiscreto.signals.values(:,3);
Tdisc300 = simDiscreto.signals.values(:,4);

u = h*ones(size(t));

% Identificação G discreto - 50s
[Gd50, Td50, Ld50] = parametrosFOPTD(Tdisc50 - Tdisc50(1), h, t(end)-t(end-1)) %t = 5
GDisc50 = tf(Gd50, [Td50 1], 'iodelay', L11);
DiscretoSim50 = lsim(GDisc50, u, t);
emqGDisc50 = mean((Tdisc50 - (DiscretoSim50+Tdisc50(1))).^2)

% Identificação G discreto - 150s
[Gd150, Td150, Ld150] = parametrosFOPTD(Tdisc50 - Tdisc50(1), h, t(end)-t(end-1)) %t = 5
GDisc150 = tf(Gd150, [Td150 1], 'iodelay', Ld150);
DiscretoSim150 = lsim(GDisc150, u, t);
emqGDisc150 = mean((Tdisc50 - (DiscretoSim150+Tdisc50(1))).^2)

% Identificação G discreto - 300s
[Gd300, Td300, Ld300] = parametrosFOPTD(Tdisc300 - Tdisc300(1), h, t(end)-t(end-1)) %t = 5
GDisc300 = tf(Gd300, [Td300 1], 'iodelay', Ld300);
DiscretoSim300 = lsim(GDisc50, u, t);
emqGDisc300 = mean((Tdisc300 - (DiscretoSim300+Tdisc300(1))).^2)

%Plots
subplot(1,3,1)
plot(t,Tdisc50, '--', t,DiscretoSim50+Tdisc50(1), '-')
grid on;
ylabel('T (°C)');
legend('Identificação', 'Ta = 50 s', 'Location', 'Best');

subplot(1,3,2)
plot(t,Tdisc150, '--', t,DiscretoSim150+Tdisc150(1))
grid on;
xlabel('t (s)');
title('Identificação dos sistema discretizados com Ta = 50, 150 e 300 s');
legend('Identificação', 'Ta = 150 s', 'Location', 'Best');

subplot(1,3,3)
plot(t,Tdisc300, '--', t,DiscretoSim300+Tdisc300(1), '-')
grid on;
legend('Identificação', 'Ta = 300 s', 'Location', 'Best');

%% Comparação das Identificações

figure
plot(t,TFjsim+TstepFj(1), '-',t,Tdisc50, ':.', t,DiscretoSim50+Tdisc50(1), '--')
title('Comparação');
grid on;
xlabel('t (s)');
legend('Modelo G11 identificado','Discretização G_d50s','Modelo G_d50s identificado', 'Location', 'Best');
ylabel('T (°C)')

figure
plot(t,Tcont, '--', t,TFjsim+TstepFj(1), '-',t,Tdisc50, t,DiscretoSim50+Tdisc50(1))
legend('Modelo Linear', 'Modelo G11 identificado','Discretização G_d50s','Modelo G_d50s identificado', 'Location', 'Best');

%EMQ entre Modelo Linear - Identificaão G11
emqML_G11 = mean((Tcont - (TFjsim+TstepFj(1))).^2)

%EMQ entre Discreto 50s - Identificaão Ts = 50s
emqG50_Gd50_I = mean((Tdisc50 - (DiscretoSim50+Tdisc50(1))).^2)
%EMQ entre real - Discretizaçao Ts = 50s
emqML_Gd50 = mean((Tcont - (Tdisc50)).^2)

%EMQ entre real - Identificaão Ts = 50s
emqML_I_Gd50 = mean((Tcont - (DiscretoSim50+Tdisc50(1))).^2)

%EMQ entre Identificaão G11 - Identificaão Ts = 50s
emqML_I_Gd50 = mean((TFjsim+TstepFj(1) - (DiscretoSim50+Tdisc50(1))).^2)

emqML_Gd150 = mean((Tcont - (Tdisc150)).^2)
emqML_Gd300 = mean((Tcont - (Tdisc300)).^2)




%Linearização CSTR
%Mat(h)euses
%% Definição das variáveis
%{
A = 0; %Área para a troca de calor 
CA = 0; %Concentração de área no reator
CAf = 0; %Concentração de área no fluxo de alimentação
Cp = 0; %Capacidade de calor em energia pora massa*temperaturahttps://www.skyscanner.com.br/hotels/search?entity_id=46994211&checkin=2025-05-28&checkout=2025-05-29&adults=2&rooms=1&sort=price&currency=BRL
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
% A11 = diff(f1, CA); %o que seria referenciado ao CAs que é atribuído no bequette, onde o subscrito s denota o valor estacionário, não adotado aq.
% A12 = diff(f1, T);
% A21 = diff(f2, CA);
% A22 = diff(f2, T);

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
du1 = 0.0; du2 = 0.00; du3 = 1; du4 = 0;

%% Extrai os dados da simulação
simData = out.simDiscreto;

% Verifica a estrutura dos dados (para debug)
% whos simData
% disp(simData)

% Extrai o vetor de tempo (assumindo que está na primeira coluna)
time = simData.time;

% Extrai as saídas (assumindo a ordem: contínuo, Ta=0.2s, Ta=0.6s, Ta=1s)
output_cont = simData.signals.values(:,1);  % Sistema contínuo
output_02s  = simData.signals.values(:,2);  % Ta = 0.2s
output_06s  = simData.signals.values(:,3);  % Ta = 0.6s
output_1s   = simData.signals.values(:,4);  % Ta = 1s

% Configurações do gráfico
figure;
hold on;
grid on;
set(gcf, 'Color', 'w');  % Fundo branco

% Plota a resposta contínua (referência)
plot(time, output_cont, 'LineWidth', 1, 'DisplayName', 'Contínuo');

% Plota as respostas discretas com marcadores
stairs(time, output_02s,  'LineWidth', 1, 'DisplayName', 'T_a = 0.2s');
stairs(time, output_06s, 'LineWidth', 1, 'DisplayName', 'T_a = 0.6s');
stairs(time, output_1s,  'LineWidth', 1, 'DisplayName', 'T_a = 1s');



% Configurações adicionais
title('Comparação de Tempos de Amostragem no CSTR');
xlabel('Tempo (s)');
ylabel('Concentração CA (kgmol/m^3)');
legend('Location', 'best');

% Ajusta limites para melhor visualização
xlim([0 time(end)]);

hold off;
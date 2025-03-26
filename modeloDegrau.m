function [K, T, L] = modeloDegrau(y, u, Ta)
%CALCULAG1 Identifica os parâmetros a1 e b1 de um modelo do tipo
%                   K*exp(-L*s)
%            G(s) = -----------
%                     Ts + 1
%          
%   Entradas:
%       - y : Vetor de dados da saída do processo
%       - u : Vetor de dados da entrada do processo
%   Saídas:
%       - K : Ganho do processo
%       - T : Constante de tempo
%       - L : Atraso

if(length(y) ~= length(u))
    error('Vetores de entrada devem ter o mesmo tamanho');
end

Theta = zeros(3,1); % = [G0 G0*L -T1]'
R = zeros(3,3);
f = Theta;

h = max(u); % Amplitude do degrau

for tau = 0:length(y)-1 % Considera o período de amostragem = 1
   phi = [h*tau*Ta; -h; y(tau+1)];  
   R = R + phi*(phi');
   A = Ta*sum(y(1:tau+1));
   f = f + phi*A;
end 

Theta = R\f;
K = Theta(1);
T = -Theta(3);
L = abs(Theta(2)/K);

end

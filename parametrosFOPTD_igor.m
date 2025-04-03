function [G0, T, L] = parametrosFOPTD(y, h, DeltaT)
% Programa para estima��o dos par�metros de um modelo de primeira ordem com
%   atraso por meio do m�todo dos m�nimos quadrados
% O tempo de amostragem � assumido como unit�rio
%
% CALCULAG1 Identifica os par�metros a1 e b1 de um modelo do tipo
%                       K
%            G(s) = -----------*exp(-L*s)
%                     Ts + 1
%          
%   Entradas:
%       - y : Vetor de dados da sa�da do processo
%       - u : Vetor de dados da entrada do processo
%   Sa�das:
%       - K : Ganho do processo
%       - T : Constante de tempo
%       - L : Atraso

% Os vetores de entrada u e y precisam ter o mesmo comprimento
% if(length(y) ~= length(u))
%     error('Vetores de entrada devem ter o mesmo tamanho');
% end

Theta = zeros(3,1); % = [G0 G0*L -T1]'
R = zeros(3,3);
f = Theta;

% h = max(u); % Amplitude do degrau

for k = 0:length(y)-1 % Considera o per�odo de amostragem = 1
   phi = [h*k*DeltaT; -h; y(k+1)];  
   R = R + phi*(phi');
   A = DeltaT*sum(y(1:k+1));
   f = f + phi*A;
end 

Theta = inv(R)*f;  %R\f
G0 = Theta(1);
T = -Theta(3);
L = abs(Theta(2)/G0);

end

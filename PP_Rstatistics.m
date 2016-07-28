clear all
close all

global lambda1 lambda2 lambda3
global data

% Dados normais
a = normrnd(120, 1, 50, 1);
b = normrnd(130, 1, 50, 1);
c = normrnd(140, 1, 50, 1);
d = normrnd(130, 1, 50, 1);
e = normrnd(120, 1, 150, 1);
f = normrnd(134, 1, 250, 1);
g = normrnd(120, 1, 50, 1);

data = vertcat(a, b, c, d, e, f, g);
time = (1:length(data))';

% Constantes de suavização exponencial
lambda1 = 0.2;
lambda2 = 0.1;
lambda3 = 0.1;

% Plotar Série temporal
plot(time, data)
    title('Série Temporal')
    xlabel('tempo')
    ylabel('dados')

% Criar vetores para solução da estatística R
X       = zeros((length(data)-9), 1);
Xf      = zeros((length(data)-9), 1);
Varf    = zeros((length(data)-9), 1);
Deltaf  = zeros((length(data)-9), 1);
R       = zeros((length(data)-9), 1);
Rc      = zeros((length(data)-9), 1);
timeX   = 1:(length(data)-9);

% Amostragem inical de tamanho 10
X(1)        = mean(data(1:10));
Xf(1)       = mean(data(1:10));
Varf(1)     = var(data(1:10));
Deltaf(1)   = 2*var(data(1:10));
R(1)        = ((2 - lambda1)*Varf(1)) / Deltaf(2);

% Solucionar estatística R
for i = 2:(length(data)-9)
    X(i)        = data(i+9);
    Varf(i)     = lambda2*(X(i) - Xf(i-1))^2 + (1 - lambda2)*Varf(i-1);
    Xf(i)       = lambda1*X(i) + (1 - lambda1)*Xf(i-1);
    Deltaf(i)   = lambda3*(X(i) - X(i-1))^2 + (1 - lambda3)*Deltaf(i-1);
    R(i)        = ((2 - lambda1)*Varf(i)) / Deltaf(i);
end

% Verificar se o estado é estacionário
%   R < 2 : dados em estado estacionário
%   R >=2 : dados não estacionários
for i = 1:length(R)
    if R(i)<2
        Rc(i)=0;
    else
        Rc(i)=mean(X);
    end
end

% Plotar série temporal e estatística Rc
plot(timeX, X, timeX, Rc);
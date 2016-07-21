clear all
close all

% Carregar dados
datarefluxflow;

plot(refluxflow(:,1), refluxflow(:,2))
    title('Série Temporal')
    xlabel('tempo')
    ylabel('dados')

% Intervalo de amostragem
si = 5;

% Constantes de suavização exponencial
lambda1 = 0.2;
lambda2 = 0.1;
lambda3 = 0.1;

% Criar vetores para solução da estatística R
X       = zeros(length(refluxflow)/si, 1);
Xf      = zeros(length(refluxflow)/si, 1);
Varf    = zeros(length(refluxflow)/si, 1);
Deltaf  = zeros(length(refluxflow)/si, 1);
R       = zeros(length(refluxflow)/si, 1);
Rc      = zeros(length(refluxflow)/si, 1);
Rcr     = zeros(length(refluxflow), 1);
timeX   = 1:(length(refluxflow));

% Amostragem inical de tamanho si
X(1)        = mean(refluxflow(1:si,2));
Xf(1)       = mean(refluxflow(1:si,2));
Varf(1)     = var(refluxflow(1:si,2));
Deltaf(1)   = 2*var(refluxflow(1:si,2));
R(1)        = ((2 - lambda1)*Varf(1)) / Deltaf(2);

% Solucionar estatística R
for i = 2:length(X)
    X(i)        = refluxflow((si+si*(i-1)),2);
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

for i = 1:length(Rc)
    Rcr((i+((si-1)*(i-1))):(i*si)) = Rc(i);
end

% Plotar série temporal e estatística Rc
plot(refluxflow(:,1), refluxflow(:,2), refluxflow(:,1), Rcr)
    title('Série Temporal')
    xlabel('tempo')
    ylabel('dados')
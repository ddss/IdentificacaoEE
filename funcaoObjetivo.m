function [ fobj ] = funcaoObjetivo( pc, serie, uyy, tipofobj, setN, PA )
% Function para avaliar a fun????o objetivo do problema de identifica????o de estados estacion??rios
%
% Entradas:
% pc: ?? um vetor linha de 0 e 1, identificando se o ponto de corte est?? ativo (0) ou inativo (1). 
% Neste, deve ser exclu??do as extremidades, ou seja, ele tem 2 posi????es a
% menos do que a serie (Ex.: A posi????o n representa a posi????o n+1 na s??rie)
%
% serie: ?? a s??rie (vetor linha) de dados hist??ricos, na qual deseja-se identificar os
% estados estacion??rios
%
% uyy: vetor linha contendo a incerteza dos pontos
%
% tipofobj: escolha da fun??ao objetivo
% 1 - R2ajustado com pondera??ao do numero de pontos em EE
% 2 - R2ajustado
% 3 - R2ajustado com phi
% setN: com avaliar o n?mero de par?metros
% setN = 1 -> N = 2*(sum(pc)+1);
% setN = 2 -> N = sum(pc)+2*(sum(pc)+1);

% TESTE:
% serie = [1,1,1,2,2,2,2,3,3,3];
% os pontos de corte est??o na posicao 3 e 7 da s??rie, logo, devem estar
% na posi????o 2 e 6 de pc.
% pc    =   [0,1,0,0,0,1,0,0]; 
% uyy   = ones(1,length(serie)).^2;

[ Residuo,NE,~,~,~,~,~,~,~,~,~,phi] = estimacao( serie, uyy, pc, PA, false );

% numero de parametros para teste:
if setN == 1
    N = 2*(sum(pc)+1);
else
    N = sum(pc)+2*(sum(pc)+1);
end

SSE = sum(Residuo.^2);
SST = sum((serie-mean(serie)).^2);

NEprojeto = 30;

if N*NE/NEprojeto < length(serie) % impedir NaN
    if tipofobj == 1
        fobj = (SSE/(length(serie)-N))/(SST/(length(serie) - 1)) + ((length(serie)-NE)/length(serie))^2;
    elseif tipofobj == 2
        fobj = (SSE/(length(serie)-N))/(SST/(length(serie) - 1));
    else
        fobj = (SSE/(length(serie)-N*NE/NEprojeto))/(SST/(length(serie) - 1))/(phi+eps);
    end
else
    if tipofobj == 1
        fobj = (SSE/10^(-10))/(SST/(length(serie) - 1)) + ((length(serie)-NE)/length(serie))^2;
    elseif tipofobj == 2
        fobj = (SSE/10^(-10))/(SST/(length(serie) - 1));
    else
        fobj = (SSE/10^(-12))/(SST/(length(serie) - 1))/phi;
    end
end
% valor da fun????o objetivo
%fobj = N*log(var(Residuo))+ 2*(N+1);
% fobj = sum(Residuo.^2)+2*(N-1);
% pensar no R2ajustado
end
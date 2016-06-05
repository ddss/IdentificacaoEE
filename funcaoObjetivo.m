function [ fobj ] = funcaoObjetivo( pc, serie, uyy )
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
% TESTE:
% serie = [1,1,1,2,2,2,2,3,3,3];
% os pontos de corte est??o na posicao 3 e 7 da s??rie, logo, devem estar
% na posi????o 2 e 6 de pc.
% pc    =   [0,1,0,0,0,1,0,0]; 
% uyy   = ones(1,length(serie)).^2;

[ Residuo,~,~,~,~,~] = estimacao( serie, uyy, pc, false );

N = 2*(sum(pc)+1);

SSE = sum(Residuo.^2);
SST = sum((serie-mean(serie)).^2);

if N < length(serie) % impedir NaN
    fobj = (SSE/(length(serie)-N))/(SST/(length(serie) - 1));
else
    fobj = (SSE/10^(-10))/(SST/(length(serie) - 1));
% valor da fun????o objetivo
%fobj = N*log(var(Residuo))+ 2*(N+1);
% fobj = sum(Residuo.^2)+2*(N-1);
% pensar no R2ajustado
end
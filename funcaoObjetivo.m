function [ fobj ] = funcaoObjetivo( pc, serie, uyy )
% Function para avaliar a função objetivo do problema de identificação de estados estacionários
%
% Entradas:
% pc: é um vetor linha de 0 e 1, identificando se o ponto de corte está ativo (0) ou inativo (1). 
% Neste, deve ser excluído as extremidades, ou seja, ele tem 2 posições a
% menos do que a serie (Ex.: A posição n representa a posição n+1 na série)
%
% serie: é a série (vetor linha) de dados históricos, na qual deseja-se identificar os
% estados estacionários
%
% uyy: vetor linha contendo a incerteza dos pontos
%
% TESTE:
% serie = [1,1,1,2,2,2,2,3,3,3];
% os pontos de corte estão na posicao 3 e 7 da série, logo, devem estar
% na posição 2 e 6 de pc.
% pc    =   [0,1,0,0,0,1,0,0]; 
% uyy   = ones(1,length(serie)).^2;

[ Residuo,~,~,~ ] = estimacao( serie, uyy, pc );

N = length(find(pc==1));

SSE = sum(Residuo.^2);
SST = sum((serie-mean(serie)).^2);
fobj = (SSE/(length(serie)-N))/(SST/(length(serie) - 1));
% valor da função objetivo
% fobj = N*log(var(Residuo))+2*(N+1);
% fobj = sum(Residuo.^2)+2*(N-1);
% pensar no R2ajustado
end
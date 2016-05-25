function [ fobj ] = funcaoObjetivo( pc, serie, uyy )
% Function para avaliar a função objetivo do problema de identificação de estados estacionários
%
% Entradas:
% pc: é um vetor linha de 0 e 1, identificando se o ponto de corte está ativo (0) ou inativo (1). 
% Neste, deve ser excluído as extermidades, ou seja, ele tem 2 posições a
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

% gerando um vetor para identificar as amostras (posições)
amostras = 1:length(serie);

% convertando o vetor uyy na matrix covariância
Uyy = diag(uyy.^2);

% Criando o vetor que indica os pontos de corte ativos (as extremidas
% sempre estão ativas
pontosAtivos = [1 find(pc==1)+1 length(serie)];

N = length(pontosAtivos)-2; % número de pontos de corte

% Avaliar as retas (Regressão de com múltiplos pontos de corte - RLMPC);
Residuo = zeros(1,length(serie)+length(pontosAtivos)-2);

for pos = 1:length(pontosAtivos)-1
    % obter os dados da reta
    dadosReta  = serie(pontosAtivos(pos):pontosAtivos(pos+1))';
    % obter os dados de x
    xDummy     = [amostras(pontosAtivos(pos):pontosAtivos(pos+1));ones(1,length(dadosReta))]';
    % matriz covariância
    Uyy_aux    = Uyy(pontosAtivos(pos):pontosAtivos(pos+1),pontosAtivos(pos):pontosAtivos(pos+1));
    % estimação dos parâmetros - WLS
    parametros = (xDummy'/(Uyy_aux)*xDummy)\xDummy'/(Uyy_aux)*dadosReta;
    
    % salvando os resíduos
    if pos == 1
        Residuo(pontosAtivos(pos):pontosAtivos(pos+1)) = (dadosReta - xDummy*parametros);
    else
        Residuo(pontosAtivos(pos)+1:pontosAtivos(pos+1)+1) = (dadosReta - xDummy*parametros);
    end
end

% valor da função objetivo
% fobj = N*log(var(Residuo))+2*(N+1);
fobj = sum(Residuo.^2) + 2*(N+1);
% pensar no R2ajustado
end
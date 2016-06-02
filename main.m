% Programa principal
% Dada uma serie, este program
%% Inicialização
clear all
close all
clc

%% Obtenção dos dados
serie = [1,1,1,2,2,2,2,3,3,3];
%pc    =   [0,1,0,0,0,1,0,0]; 
uyy   = ones(1,length(serie)).^2;

%% Otimização

% número de variáveis de decisão
nvars  = length(serie)-2;

% Limite inferior
LB = zeros(1,nvars);

% Limite superior
UB = ones(1,nvars);

% definido que todas as variáveis de decisão
% são números inteiros
IntCon = 1:nvars;
% definindo o vetor de opções
options= [];

% Algoritmo genético
%                          ga(fitnessfcn                        ,nvars,...
%                             A,b,[],[],LB,UB,nonlcon,IntCon,options)
[x,fval,exitflag,output] = ga(@(pc) funcaoObjetivo(pc,serie,uyy),nvars,...
                              [],[],[],[],LB,UB,[],IntCon,options);
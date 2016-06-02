% Programa principal
% Dada uma serie, este program
%% Inicialização
clear all
close all
clc

%% Obtenção dos dados
serie = [1.1,0.9,1.1,0.9,2,2.1,1.8,2,1.8,3.1,3,2.9,3,4,5,6,7];
%pc    =   [0,1,0,0,0,1,0,0]; 
uyy   = ones(1,length(serie)).^2;

nRetas = 5;
%% Otimização

% número de variáveis de decisão
nvars  = nRetas;

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
[pontosCorte,fval,exitflag,output] = ga(@(pc) funcaoObjetivo(pc,serie,uyy),nvars,...
                              [],[],[],[],LB,UB,[],IntCon,options);
                          
[ ~,retas,parametros,Uparametros ] = estimacao( serie, uyy, pontosCorte, true );

retas
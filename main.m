% Programa principal
% Dada uma serie, este programa identifica os Estados Estacion?rios
%% Inicializa??o
clear all
close all
clc

%% Obten??o dos dados

serie = [1.02,0.95,1.01,0.97,2,2.01,1.89,2,1.89,3.01,3,2.95,3,4,5,6,7,4.1,4.08,4.12,4.06,4.09];
%serie = [1 1 1 2 2 2 3 3 3] ;
%pc    =   [0,1,0,0,0,1,0,0]; 
uyy   = 0.01*ones(1,length(serie)).^2;

nRetas = 12;
%% Otimiza??o

% n?mero de vari?veis de decis?o
nvars  = length(serie)-2;

% Limite inferior
LB = zeros(1,nvars);

% Limite superior
UB = ones(1,nvars);

% definido que todas as vari?veis de decis?o
% s?oo n?meros inteiros
IntCon = 1:nvars;
% definindo o vetor de op??es
options= [];

% Algoritmo gen?tico
%                          ga(fitnessfcn                        ,nvars,...
%                             A,b,[],[],LB,UB,nonlcon,IntCon,options)
[pontosCorte,fval,exitflag,output] = ga(@(pc) funcaoObjetivo(pc,serie,uyy),nvars,...
                              [],[],[],[],LB,UB,@(pc) restricao(pc,nRetas),IntCon,options);
                          
[ residuo,retas,pontosAtivos,parametros,Uparametros,Residuos ] = estimacao( serie, uyy, pontosCorte, true );

%% Figuras
amostras = 1:length(serie);

figure()
ax = subplot(1,1,1);
hold(ax,'on')
plot(amostras,serie)
for pos = 1:length(pontosAtivos)-1
    x = [amostras(pontosAtivos(pos):pontosAtivos(pos+1));ones(1,length(retas{pos}))]';
    y = x*parametros{pos};
    plot(x(:,1),y,'--');
end

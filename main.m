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
uyy   = 0.001*ones(1,length(serie)).^2;

nRetas = 12;

PA = 0.95;

projeto = 'Teste';
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
                          
[ residuo,retas,pontosAtivos,parametros,Uparametros,Residuos,FuncaoObjetivo ] = estimacao( serie, uyy, pontosCorte, true );

%% C?lculo da regi?o de abrang?ncia e verifica??o das retas candidatas
% Aqui est? apenas se calculando os pontos extremos da elipse
% Considerando que os par?metros seguem uma distribui??o normal

posCandidatasEE = 1;
for pos = 1:length(pontosAtivos)-1
    
    Fisher = finv(PA,2,(length(retas{pos})-2));
    raioEllip = FuncaoObjetivo{pos}*(2/(length(retas{pos})-2)*Fisher);

    parametro_aux = parametros{pos};
        
    invUparametros = inv(Uparametros{pos});
    
    fator = invUparametros(1,1)/(invUparametros(1,2) + eps); % eps evita NaN quando a covari?ncia ? zero.
    delta = sqrt(raioEllip/(fator^2*invUparametros(2,2)-2*fator*invUparametros(1,2)+invUparametros(1,1)));
    coordenadas_x = [parametro_aux(1)+delta,parametro_aux(1)-delta];
    coordenadas_y = [parametro_aux(2)-delta*fator,parametro_aux(2)+delta*fator];

    fator = invUparametros(2,2)/(invUparametros(1,2) + eps); % eps evita NaN quando a covari?ncia ? zero.
    delta = sqrt(raioEllip/(fator^2*invUparametros(1,1)-2*fator*invUparametros(1,2)+invUparametros(2,2)));
    coordenadas_y = [coordenadas_y [parametro_aux(2)+delta,parametro_aux(2)-delta]];
    coordenadas_x = [coordenadas_x [parametro_aux(1)-delta*fator,parametro_aux(1)+delta*fator]];

    % Obten??o das posi??es das retas candidatas a EE
    % - verificar se e ellipse do par?metro a, cruza o zero.
    if and(any(coordenadas_x>0), any(coordenadas_x<0))
        CandidatasEE(posCandidatasEE) = pos;
        posCandidatasEE = posCandidatasEE+1;
    end
end

%% Testes estat?sticos
stat_test = {};
stat_test.residuo = {};
% p-valor do ttest
stat_test.residuo.media = ones(1,length(CandidatasEE));
% p-valor Ljung-Box Q-test
stat_test.serie.autocorrLjung = ones(1,length(CandidatasEE));
% aleatoriedade
stat_test.serie.random = ones(1,length(CandidatasEE));

for pos = CandidatasEE
    % RES?DUOS
    % Teste para verificar se a m?dia do res?duo de regress?o ? zero
    [~, stat_test.residuo.media(pos)] = ttest(Residuos{pos});
   
    % S?RIE DE DADOS
    desvio = retas{pos} - mean(retas{pos});
    % autocorrela??o
    [~,stat_test.serie.autocorrLjung(pos)] = lbqtest(desvio);
    % aleatoriedade
    if length(desvio)>5
        [~,stat_test.serie.random(pos)] = vratiotest(desvio);
    else
        stat_test.serie.random(pos) = NaN;
    end
end

%% Configura??es de diret?rio:
% criando folder
for pos = CandidatasEE
    folder = strcat('./',projeto,'/','reta_',num2str(pos));
    if not(exist(folder,'dir'))
        mkdir(folder)
    end
end

%% Figuras

amostras = 1:length(serie);

% S?rie de Dados
figure()
ax = subplot(1,1,1);
hold(ax,'on')
plot(amostras,serie,'.','MarkerSize',15)
for pos = 1:length(pontosAtivos)-1
    x = [amostras(pontosAtivos(pos):pontosAtivos(pos+1));ones(1,length(retas{pos}))]';
    y = x*parametros{pos};
    plot(x(:,1),y,'--','LineWidth',1.5);
end
xlabel('Amostra','FontSize',12)
ylabel('Serie','FontSize',12)
set(ax,'FontSize',12)
saveas(gcf, strcat('./',projeto,'/','geral.png'))

for pos = CandidatasEE
    folder = strcat('./',projeto,'/','reta_',num2str(pos));
    
    % RES?DUOS

    % BOXPLOT
    
    fig = figure('Visible','off');
    ax = subplot(1,1,1);
    boxplot(Residuos{pos})
   
    set(ax,'FontSize',12)
    
    title(num2str(pos))
    saveas(gcf,strcat(folder,'/residuos-boxplot.png'))
    close(fig)
    
    % DISPERS?O
    
    fig = figure('Visible','off');
    ax = subplot(1,1,1);
    hold(ax,'on')
    
    x = amostras(pontosAtivos(pos):pontosAtivos(pos+1));
    meanRes = mean(Residuos{pos});
    stdRes  = std(Residuos{pos});
    
    fatorK = tinv(PA,length(Residuos{pos})-1);

    p1 = plot(x,Residuos{pos},'.','MarkerSize',15);  
    p2 = plot([x(1) x(end)],[meanRes meanRes],'r--','LineWidth',2);
    p3 = plot([x(1) x(end)],[meanRes+fatorK*stdRes meanRes+fatorK*stdRes],'r-.','LineWidth',2);
    p4 = plot([x(1) x(end)],[meanRes-fatorK*stdRes meanRes-fatorK*stdRes],'r-.','LineWidth',2);
    
    leg = legend([p1,p2,p3],{'res?duos','m?dia','intervalo abrang?ncia'});
    
    set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')
    
    xlabel('Amostra','FontSize',12)
    ylabel('Res?duo','FontSize',12)
    
    set(ax,'FontSize',12)
    
    title(num2str(pos))
    saveas(gcf,strcat(folder,'/residuos-tendencia.png'))
    close(fig)
    
    % S?RIE
    
    % BOXPLOT
    fig = figure('Visible','off');
    ax = subplot(1,1,1);
    boxplot(retas{pos})
   
    set(ax,'FontSize',12)
    
    title(num2str(pos))
    saveas(gcf,strcat(folder,'/dados-boxplot.png'))
    close(fig)
    
    % DISPERS?O
    fig = figure('Visible','off');
    ax = subplot(1,1,1);
    hold(ax,'on')
    
    x = amostras(pontosAtivos(pos):pontosAtivos(pos+1));
    meanReta = mean(retas{pos});
    stdReta  = std(retas{pos});
    
    fatorK = tinv(PA,length(retas{pos})-1);

    p1 = plot(x,retas{pos},'.','MarkerSize',15);  
    p2 = plot([x(1) x(end)],[meanReta meanReta],'r--','LineWidth',2);
    p3 = plot([x(1) x(end)],[meanReta+fatorK*stdReta meanReta+fatorK*stdReta],'r-.','LineWidth',2);
    p4 = plot([x(1) x(end)],[meanReta-fatorK*stdReta meanReta-fatorK*stdReta],'r-.','LineWidth',2);
    
    leg = legend([p1,p2,p3],{'res?duos','m?dia','intervalo abrang?ncia'});
    
    set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')
    
    xlabel('Amostra','FontSize',12)
    ylabel('Dados','FontSize',12)
    
    set(ax,'FontSize',12)
    
    title(num2str(pos))
    saveas(gcf,strcat(folder,'/dados-tendencia.png'))
    close(fig)
    
    % AUTOCORRELA??O
    fig = figure('Visible','off');
    ax = subplot(1,1,1);
    autocorr(retas{pos});

    xlabel('Lag - Amostra','FontSize',12)
    ylabel('Autocorrelacao','FontSize',12)
    
    set(ax,'FontSize',12)

    title(num2str(pos))
    saveas(gcf,strcat(folder,'/dados-autocorr.png'))
    close(fig)
    
end

%% Relat?rios
for pos = CandidatasEE
    
    folder = strcat('./',projeto,'/','reta_',num2str(pos));
   
    fileID = fopen(strcat(folder,'/','Testes.txt'),'w');
    
    fprintf(fileID,'Testes estatisticos para avaliar estacionaridade \n');
    fprintf(fileID,'Intervalo de amostra: %d a %d \n',pontosAtivos(pos),pontosAtivos(pos+1));
    fprintf(fileID,'M?dia dos dados: %0.3g \n',mean(retas{pos}));    
    fprintf('\n');
    fprintf(fileID,'=============RES?DUOS========== \n');
    
    if stat_test.residuo.media(pos)>(1-PA)
        h = 'n?o se pode rejeitar a hip?tese de que a m?dia seja zero.'; 
    else
        h = 'a m?dia n?o ? zero.';
    end
    
    fprintf(fileID,' - M?dia zero(ttest): %0.3g \n',stat_test.residuo.media(pos));
    fprintf(fileID,' - - %s \n',h);
    fprintf('\n');
    fprintf(fileID,'==============S?RIE=========== \n');
    
    if stat_test.serie.autocorrLjung(pos)>(1-PA)
        h = 'n?o se pode rejeitar a hip?tese de que n?o h? autocorrela??o na s?rie'; 
    else
        h = 'h? autocorrela??o na s?rie.';
    end
    
    fprintf(fileID,' - Autocorrela??o - Ljung-Box Q-test: %.3g \n',...
        stat_test.serie.autocorrLjung(pos));
    fprintf(fileID,' - - %s \n',h);
    
    if stat_test.serie.random(pos)>(1-PA)
        h = 'n?o se pode rejeitar a hip?tese de que a s?rie seja aleat?ria'; 
    else
        h = 'a s?rie n?o ? aleat?ria.';
    end
    
    fprintf(fileID,' - Aleatoriedade: %.3g \n',...
        stat_test.serie.random(pos));
    fprintf(fileID,' - - %s \n',h);
    fclose(fileID);
end



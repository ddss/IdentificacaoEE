% Programa principal
% Dada uma serie, este programa identifica os Estados Estacion?rios
%% Inicializa??o
clear all
close all
clc

%% Obten??o dos dados

%serie = [1.02,0.95,1.01,0.97,2,2.01,1.89,2,1.89,3.01,3,2.95,3,4,5,6,7,4.1,4.08,4.12,4.06,4.09];
%serie = [1 1 1 2 2 2 3 3 3] ;
serie = [229.604599000000;226.687622100000;227.820114100000;229.907409700000;228.141418500000;230.153366100000;231.065383900000;229.407714800000;230.797912600000;229.163162200000;231.249160800000;232.265884400000;231.252777100000;231.539047200000;233.166824300000;233.332122800000;232.809249900000;233.270858800000;230.996933000000;232.833908100000;232.798996000000;233.840606700000;234.772781400000;232.139968900000;232.471939100000;234.406753500000;232.941848800000;231.381149300000;233.731506300000;233.150878900000;231.146347000000;232.640945400000;232.733856200000;233.160491900000;230.095550500000;231.057312000000;230.029037500000;233.827560400000;233.516861000000;233.714141800000;231.120559700000;232.040817300000;235.269226100000;233.175415000000;228.865325900000;232.330902100000;230.694442700000;231.959671000000;230.746398900000;232.896026600000;231.356857300000;235.484481800000;234.790420500000;231.641677900000;235.882095300000;233.969863900000;234.795013400000;232.874740600000;235.759841900000;234.014358500000;234.339447000000;231.794479400000;229.508087200000;230.994903600000;229.290115400000;229.806564300000;229.635940600000;228.726211500000;228.784912100000;230.039459200000;227.657623300000;228.108779900000;227.601760900000;226.256195100000;225.171173100000;228.148468000000;225.396652200000;224.723388700000;226.428482100000;222.853759800000;221.912399300000;224.512466400000;221.798858600000;222.577880900000;221.570877100000;219.529144300000;218.106185900000;220.487564100000;218.010437000000;219.559188800000;222.021499600000;216.721023600000;216.899734500000;221.097946200000;223.148620600000];
serie = serie';

uyy   = 1*ones(1,length(serie));

nRetas = 5;
%serie = [1 1 1 2 2 2 3 3 3] ;

PA = 0.90;

setN = 2;

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
options= gaoptimset('UseParallel','always','TolFun',1e-12,'TolCon',1e-12,...
    'PopulationSize',500);
            %gaoptimset('Generations',1000,...
            %        'PopulationSize',35,'TolFun',1e-7,'TolCon',1e-7,...
            %        'UseParallel','always');PA

% Algoritmo gen?tico
%                          ga(fitnessfcn                        ,nvars,...
%                             A,b,[],[],LB,UB,nonlcon,IntCon,options)
[pontosCorte,fval,exitflag,output] = ga(@(pc) funcaoObjetivo(pc,serie,uyy,setN, PA),nvars,...
                              [],[],[],[],LB,UB,@(pc) restricao(pc,nRetas),IntCon,options);
                          
[ residuo,~, retas,pontosAtivos,parametros,Uparametros,Residuos,FuncaoObjetivo, CandidatasEE ] = estimacao( serie, uyy, pontosCorte, PA, true );

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
    x = [1:length(retas{pos});ones(1,length(retas{pos}))]';
    y = x*parametros{pos};
    plot(amostras(pontosAtivos(pos):pontosAtivos(pos+1)),y,'--','LineWidth',1.5);
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



% Programa principal
% Dada uma serie, este programa identifica os Estados Estacion?rios
%% Inicializa??o
%matlabpool close  % force Paralelo

clear all
close all
clc


%matlabpool open

% Habilitar o arquivo de teste
series_teste

%% Obten??o dos dados

%serie = [1.02,0.95,1.01,0.97,2,2.01,1.89,2,1.89,3.01,3,2.95,3,4,5,6,7,4.1,4.08,4.12,4.06,4.09];
%serie = [1 1 1 2 2 2 3 3 3] ;
%serie = [229.604599000000;226.687622100000;227.820114100000;229.907409700000;228.141418500000;230.153366100000;231.065383900000;229.407714800000;230.797912600000;229.163162200000;231.249160800000;232.265884400000;231.252777100000;231.539047200000;233.166824300000;233.332122800000;232.809249900000;233.270858800000;230.996933000000;232.833908100000;232.798996000000;233.840606700000;234.772781400000;232.139968900000;232.471939100000;234.406753500000;232.941848800000;231.381149300000;233.731506300000;233.150878900000;231.146347000000;232.640945400000;232.733856200000;233.160491900000;230.095550500000;231.057312000000;230.029037500000;233.827560400000;233.516861000000;233.714141800000;231.120559700000;232.040817300000;235.269226100000;233.175415000000;228.865325900000;232.330902100000;230.694442700000;231.959671000000;230.746398900000;232.896026600000;231.356857300000;235.484481800000;234.790420500000;231.641677900000;235.882095300000;233.969863900000;234.795013400000;232.874740600000;235.759841900000;234.014358500000;234.339447000000;231.794479400000;229.508087200000;230.994903600000;229.290115400000;229.806564300000;229.635940600000;228.726211500000;228.784912100000;230.039459200000;227.657623300000;228.108779900000;227.601760900000;226.256195100000;225.171173100000;228.148468000000;225.396652200000;224.723388700000;226.428482100000;222.853759800000;221.912399300000;224.512466400000;221.798858600000;222.577880900000;221.570877100000;219.529144300000;218.106185900000;220.487564100000;218.010437000000;219.559188800000;222.021499600000;216.721023600000;216.899734500000;221.097946200000;223.148620600000];

serie = [EE1 EE2 EE3 EE5 EE4];

uyy   = 1*ones(1,length(serie));

nRetas = 9;
%serie = [1 1 1 2 2 2 3 3 3] ;

PA = 0.95;

tipofobj = 3;
setN = 2;

projeto = 'Teste-phi-mod';

repetirOtimizacao = 10;

%% Otimiza??o

% numero de variaveis de decisao
% -3 -> pois os limites sempre est?o ativos e o segundo ponto nao pode ser
% ponto de corte
% metade -> para evitar retas de 1 ponto (estimacao preenche o resto)
nvars  = length(serie)-3;

% Limite inferior
LB = zeros(1,nvars);

% Limite superior
UB = ones(1,nvars);

% definido que todas as vari?veis de decis?o
% s?oo n?meros inteiros
IntCon = 1:nvars;
% definindo o vetor de op??es
options= gaoptimset('UseParallel','always','TolFun',1e-12,'TolCon',1e-12,...
    'PopulationSize',300,'Generations',500);
            %gaoptimset('Generations',1000,...
            %        'PopulationSize',35,'TolFun',1e-7,'TolCon',1e-7,...
            %        'UseParallel','always');PA

% Algoritmo gen?tico
%                          ga(fitnessfcn                        ,nvars,...
%                             A,b,[],[],LB,UB,nonlcon,IntCon,options)


for cont = 1:repetirOtimizacao
    cont
    [pontosCorte(cont,:),fval(cont),~, ~] = ga(@(pc) funcaoObjetivo(pc,serie,uyy, tipofobj, setN, PA),nvars,...
                              [],[],[],[],LB,UB,@(pc) restricao(pc,nRetas),IntCon,options);
end

[fval,indice] = min(fval);

pontosCorte = pontosCorte(indice,:);

[ residuo,~, retas,pontosInicioAtivos, pontosFimAtivos, parametros,Uparametros,Residuos,FuncaoObjetivo, CandidatasEE,IndiceEstimacao,~ ] = estimacao( serie, uyy, pontosCorte, PA, true );

%% Testes estat?sticos
stat_test = {};
stat_test.residuo = {};
% p-valor do ttest
stat_test.residuo.media = ones(1,length(pontosInicioAtivos));
% p-valor Ljung-Box Q-test
stat_test.serie.autocorrLjung = ones(1,length(pontosInicioAtivos));
% aleatoriedade
stat_test.serie.random = ones(1,length(pontosInicioAtivos));
% normalidade - Kolmogorov-Smirnov
stat_test.serie.normks = ones(1,length(pontosInicioAtivos));
% normalidade - Lilliefors
stat_test.serie.normlil = ones(1,length(pontosInicioAtivos));


for pos = 1:length(pontosInicioAtivos)
    
    % RES?DUOS
    % Teste para verificar se a m?dia do res?duo de regress?o ? zero
    [~, stat_test.residuo.media(pos)] = ttest(Residuos{pos});
   
    % S?RIE DE DADOS
    desvio = retas{pos} - mean(retas{pos});
    % autocorrela??o
    if length(desvio)>2
        [~,stat_test.serie.autocorrLjung(pos)] = lbqtest(desvio);
    else
        stat_test.serie.autocorrLjung(pos) = NaN;
    end
    % aleatoriedade
    if length(desvio)>5
        [~,stat_test.serie.random(pos)] = vratiotest(desvio);
    else
        stat_test.serie.random(pos) = NaN;
    end
    % normalidade - Kolmogorov-Smirnov 
    if length(desvio)>2
        [~,stat_test.serie.normks(pos)] = kstest(desvio);
    else
        stat_test.serie.normks(pos) = NaN;
    end
    
    % normalidade - Lilliefors
    if length(desvio)>5
        [~,stat_test.serie.normlil(pos)] = lillietest(desvio);
    else
        stat_test.serie.normlil(pos) = NaN;
    end

end

%% Configuracoes de diretorio:
% criando folder
for pos = 1:length(retas)
    folder = strcat('./',projeto,'/','reta_',num2str(pos));
    if not(exist(folder,'dir'))
        mkdir(folder)
    end
end

%% Figuras

pontosEllip = {};

amostras = 1:length(serie);

% S?rie de Dados
figure()
ax = subplot(1,1,1);
hold(ax,'on')
plot(amostras,serie,'.','MarkerSize',15)
coresretas = repmat(['m','g','k','c'],1,ceil(length(pontosInicioAtivos)/4));

for pos = 1:length(retas)
    if IndiceEstimacao{pos}
        x = [1:length(retas{pos});ones(1,length(retas{pos}))]';
        y = x*parametros{pos};
        plot(amostras(pontosInicioAtivos(pos):pontosFimAtivos(pos)),y,'--','LineWidth',1.5,'Color',coresretas(pos));
    end
end
xlabel('Amostra','FontSize',12)
ylabel('Serie','FontSize',12)
set(ax,'FontSize',12)
saveas(gcf, strcat('./',projeto,'/','geral.png'))

for pos = 1:length(retas)
    if IndiceEstimacao{pos}
        folder = strcat('./',projeto,'/','reta_',num2str(pos));

        % ======= RES?DUOS ========

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

        x = amostras(pontosInicioAtivos(pos):pontosFimAtivos(pos));
        meanRes = mean(Residuos{pos});
        stdRes  = std(Residuos{pos});

        fatorK = tinv(PA,length(Residuos{pos})-1);

        p1 = plot(x,Residuos{pos},'.','MarkerSize',15);  
        p2 = plot([x(1) x(end)],[meanRes meanRes],'r--','LineWidth',2);
        p3 = plot([x(1) x(end)],[meanRes+fatorK*stdRes meanRes+fatorK*stdRes],'r-.','LineWidth',2);
        p4 = plot([x(1) x(end)],[meanRes-fatorK*stdRes meanRes-fatorK*stdRes],'r-.','LineWidth',2);

        leg = legend([p1,p2,p3],{'residuos','media','intervalo abrangencia'});

        set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')

        xlabel('Amostra','FontSize',12)
        ylabel('Res?duo','FontSize',12)

        set(ax,'FontSize',12)

        title(num2str(pos))
        saveas(gcf,strcat(folder,'/residuos-tendencia.png'))
        close(fig)

        % ========= S?RIE ===========

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

        x = amostras(pontosInicioAtivos(pos):pontosFimAtivos(pos));
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

        % ========= REGIAO ===========

        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        hold(ax,'on')

        Fisher = finv(PA,2,(length(retas{pos})-2));
        aspect = FuncaoObjetivo{pos}*(2/(length(retas{pos})-2)*Fisher);
        pontosEllip_aux = covellipse(parametros{pos},Uparametros{pos},aspect);

        pontosEllip{pos} = pontosEllip_aux;

        pltellip = plot(pontosEllip_aux(:,1),pontosEllip_aux(:,2),'b-','LineWidth',2);
        p_aux = parametros{pos};
        plot(p_aux(1),p_aux(2),'ro','MarkerSize',6)

        pltcoefang = plot([0 0],[min(pontosEllip_aux(:,2)),max(pontosEllip_aux(:,2))],'k-.','LineWidth',1.5);
        pltmedia = plot([min(pontosEllip_aux(:,1)) max(pontosEllip_aux(:,1))],[mean(retas{pos}),mean(retas{pos})],'k-.','LineWidth',2);

        xlabel('Coef. Angular','FontSize',12)
        ylabel('Coef. Linear','FontSize',12)

        set(ax,'FontSize',12)

        leg = legend([pltellip,pltcoefang,pltmedia],{'regiao abrangencia','zero do coef. angular','media dos dados'});

        title(num2str(pos))
        saveas(gcf,strcat(folder,'/regiao-abrangencia.png'))
        close(fig)
    end
end

% ========= REGIAO COMPARACAO ===========
    
fig = figure('Visible','off');
ax = subplot(1,1,1);
hold(ax,'on')
    
for pos = CandidatasEE
    if IndiceEstimacao{pos}
        pontosEllip_aux = pontosEllip{pos};

        plot(pontosEllip_aux(:,1),pontosEllip_aux(:,2),'LineWidth',2);
        p_aux = parametros{pos};
        plot(p_aux(1),p_aux(2),'ro','MarkerSize',6)
    end
end

limites_x = get(gca,'xlim');

xlim([limites_x(1) - (1/10)*(limites_x(2)-limites_x(1)) limites_x(2)])

for pos = CandidatasEE
    
    p_aux = parametros{pos};

    x_aux = limites_x(1) - (1/15)*(limites_x(2)-limites_x(1));
    
    text(x_aux,p_aux(2),strcat('reta  ',num2str(pos)))

    plot([x_aux p_aux(1)],[p_aux(2) p_aux(2)],'--')
end

limites = get(gca,'ylim');

plot([0 0],[limites(1),limites(2)],'k-.','LineWidth',2);

xlabel('Coef. Angular','FontSize',12)
ylabel('Coef. Linear','FontSize',12)
    
set(ax,'FontSize',12)

saveas(gcf,strcat('./',projeto,'/','regioes-comparacao.png'))
close(fig)

%% Relat?rios
for pos = 1:length(pontosInicioAtivos)
    
    folder = strcat('./',projeto,'/','reta_',num2str(pos));
   
    fileID = fopen(strcat(folder,'/','Testes.txt'),'w');
    
    fprintf(fileID,'Testes estatisticos para avaliar estacionaridade \n');
    fprintf(fileID,'Intervalo de amostra: %d a %d \n',pontosInicioAtivos(pos),pontosFimAtivos(pos));
    fprintf(fileID,'M?dia dos dados: %0.3f \n',mean(retas{pos}));    
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
    
    if stat_test.serie.normks(pos)>(1-PA)
        h = 'n?o se pode rejeitar a hip?tese de que a s?rie seja normal'; 
    else
        h = 'a s?rie n?o ? ?? normal.';
    end
        
    fprintf(fileID,' - Normalidade: %.3g \n',...
        stat_test.serie.normks(pos));
    fprintf(fileID,' - - %s \n',h);
    
    if stat_test.serie.normlil(pos)>(1-PA)
        h = 'n?o se pode rejeitar a hip?tese de que a s?rie seja normal'; 
    else
        h = 'a serie nao e normal.';
    end
        
    fprintf(fileID,' - Normalidade: %.3g \n',...
        stat_test.serie.normlil(pos));
    fprintf(fileID,' - - %s \n',h);
    
    fprintf(fileID,'\n');
    
    if isempty(find(CandidatasEE==pos, 1));
        fprintf(fileID,'Nao e candidata a EE pelo criterio da regiao de abrangencia');
    else
        fprintf(fileID,'E candidata a EE pelo criterio da regiao de abrangencia');
    end
    
    fclose(fileID);
end
%matlabpool close


% Programa principal
% Dada uma serie, este programa identifica os Estados Estacionarios
%% Inicializacao
clear all
close all
clc

matlabpool open

% Habilitar o arquivo de teste
series_teste

%% Obtencao dos dados

%serie = [1.02,0.95,1.01,0.97,2,2.01,1.89,2,1.89,3.01,3,2.95,3,4,5,6,7,4.1,4.08,4.12,4.06,4.09];
%serie = [1 1 1 2 2 2 3 3 3] ;
%serie = [229.604599000000;226.687622100000;227.820114100000;229.907409700000;228.141418500000;230.153366100000;231.065383900000;229.407714800000;230.797912600000;229.163162200000;231.249160800000;232.265884400000;231.252777100000;231.539047200000;233.166824300000;233.332122800000;232.809249900000;233.270858800000;230.996933000000;232.833908100000;232.798996000000;233.840606700000;234.772781400000;232.139968900000;232.471939100000;234.406753500000;232.941848800000;231.381149300000;233.731506300000;233.150878900000;231.146347000000;232.640945400000;232.733856200000;233.160491900000;230.095550500000;231.057312000000;230.029037500000;233.827560400000;233.516861000000;233.714141800000;231.120559700000;232.040817300000;235.269226100000;233.175415000000;228.865325900000;232.330902100000;230.694442700000;231.959671000000;230.746398900000;232.896026600000;231.356857300000;235.484481800000;234.790420500000;231.641677900000;235.882095300000;233.969863900000;234.795013400000;232.874740600000;235.759841900000;234.014358500000;234.339447000000;231.794479400000;229.508087200000;230.994903600000;229.290115400000;229.806564300000;229.635940600000;228.726211500000;228.784912100000;230.039459200000;227.657623300000;228.108779900000;227.601760900000;226.256195100000;225.171173100000;228.148468000000;225.396652200000;224.723388700000;226.428482100000;222.853759800000;221.912399300000;224.512466400000;221.798858600000;222.577880900000;221.570877100000;219.529144300000;218.106185900000;220.487564100000;218.010437000000;219.559188800000;222.021499600000;216.721023600000;216.899734500000;221.097946200000;223.148620600000];
%serie = serie';
serie = [EE1 EE2 EE3 EE4 EE5];
%serie = serieTEST;

uyy   = 1*ones(1,length(serie));

nRetas = 15;

NEprojeto = 30; % numero de pontos minimos para ser EE

PA = 0.95;

tipofobj = 3;
setN = 2;

projeto = 'Teste-phi-mod';

repetirOtimizacao = 2;

%% Otimizacao

% numero de variaveis de decisao
% -3 -> pois os limites sempre estao ativos e o segundo ponto nao pode ser
% ponto de corte
% metade -> para evitar retas de 1 ponto (estimacao preenche o resto)
nvars  = length(serie)-3;

% Limite inferior
LB = zeros(1,nvars);

% Limite superior
UB = ones(1,nvars);

% definido que todas as variaveis de decisao
% sao numeros inteiros
IntCon = 1:nvars;
% definindo o vetor de opcoes
options= gaoptimset('UseParallel','always','TolFun',1e-12,'TolCon',1e-12,...
    'PopulationSize',300,'Generations',1000);
            %gaoptimset('Generations',1000,...
            %        'PopulationSize',35,'TolFun',1e-7,'TolCon',1e-7,...
            %        'UseParallel','always');PA

% Algoritmo genetico
%                          ga(fitnessfcn                        ,nvars,...
%                             A,b,[],[],LB,UB,nonlcon,IntCon,options)
residuo            = {};
retas              = {};
pontosInicioAtivos = {};
pontosFimAtivos    = {};
parametros         = {};
Uparametros        = {};
Residuos           = {};
FuncaoObjetivo     = {};
CandidatasEE       = {};
IndiceEstimacao    = {};
for cont = 1:repetirOtimizacao
    cont
    [pontosCorte(cont,:),fval(cont),~, ~] = ga(@(pc) funcaoObjetivo(pc,serie,uyy, tipofobj, setN, NEprojeto, PA),nvars,...
                              [],[],[],[],LB,UB,@(pc) restricao(pc,nRetas),IntCon,options);

    % obtendo os resultados das retas
    [ residuo{cont},~, retas{cont},pontosInicioAtivos{cont}, pontosFimAtivos{cont}, parametros{cont},Uparametros{cont},Residuos{cont},FuncaoObjetivo{cont}, CandidatasEE{cont}, IndiceEstimacao{cont},~ ] = estimacao( serie, uyy, pontosCorte(cont,:), PA, NEprojeto, true );
end
%% Escolha do ponto de trabalho para os testes, graficos e relatorios

[~,indice_otimo] = min(fval);

%% Testes estatisticos
stat_test = {};
stat_test.residuo = {};
% p-valor do ttest
stat_test.residuo.media = ones(1,length(pontosInicioAtivos{indice_otimo}));
% p-valor Ljung-Box Q-test
stat_test.serie.autocorrLjung = ones(1,length(pontosInicioAtivos{indice_otimo}));
% aleatoriedade
stat_test.serie.random = ones(1,length(pontosInicioAtivos{indice_otimo}));
% normalidade - Kolmogorov-Smirnov
stat_test.serie.normks = ones(1,length(pontosInicioAtivos{indice_otimo}));
% normalidade - Lilliefors
stat_test.serie.normlil = ones(1,length(pontosInicioAtivos{indice_otimo}));


for pos = 1:length(pontosInicioAtivos{indice_otimo})
    
    % RESIDUOS
    % Teste para verificar se a media do residuo de regressao e zero
    [~, stat_test.residuo.media(pos)] = ttest(Residuos{indice_otimo}{pos});
   
    % SERIE DE DADOS
    desvio = retas{indice_otimo}{pos} - mean(retas{indice_otimo}{pos});
    % autocorrelacao
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
        [~,stat_test.serie.normks(pos)] = kstest(desvio/std(desvio));
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
folder_retas = {};
for pos = 1:length(retas{indice_otimo})
    folder_retas{pos} = strcat('./',projeto,'/',num2str(indice_otimo),'/','reta_',num2str(pos));
    if not(exist(folder_retas{pos},'dir'))
        mkdir(folder_retas{pos})
    end
end

%% Figuras
% Configuração das figuras
  
pontosEllip = {};

amostras = 1:length(serie);

resolucao = '-r300';

unidadeFig = 'inches';

PaperPosition = [0 0 4*2 3*2];

PaperSize = [8 6];

Driver = '-dpdf';

% Execucao
% Serie de Dados
fig = figure('Visible','off');
ax = subplot(1,1,1);
hold(ax,'on')
plot(amostras,serie,'.','MarkerSize',15)
coresretas = repmat(['m','g','k','c'],1,ceil(length(pontosInicioAtivos{indice_otimo})/4));

for pos = 1:length(retas{indice_otimo})
    if IndiceEstimacao{indice_otimo}{pos}
        x = [1:length(retas{indice_otimo}{pos});ones(1,length(retas{indice_otimo}{pos}))]';
        y = x*parametros{indice_otimo}{pos};
        plot(amostras(pontosInicioAtivos{indice_otimo}(pos):pontosFimAtivos{indice_otimo}(pos)),y,'--','LineWidth',1.5,'Color',coresretas(pos));
    end
end
xlabel('Amostra','FontSize',12)
ylabel('Serie','FontSize',12)
set(ax,'FontSize',12)

set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
print(fig,Driver,strcat('./',projeto,'/','geral.pdf'),resolucao)
close(fig)


for pos = 1:length(retas{indice_otimo})
    if IndiceEstimacao{indice_otimo}{pos}
        % ======= RESIDUOS ========

        % BOXPLOT

        fig = figure('Visible','off');
                
        ax = subplot(1,1,1);
        boxplot(Residuos{indice_otimo}{pos},'Labels',{strcat('reta ',num2str(pos))})
        set(ax,'FontSize',12)

        title(num2str(pos))
        
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/residuos-boxplot.pdf'),resolucao)
        close(fig)

        % DISPERSAO

        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        hold(ax,'on')

        x = amostras(pontosInicioAtivos{indice_otimo}(pos):pontosFimAtivos{indice_otimo}(pos));
        meanRes = mean(Residuos{indice_otimo}{pos});
        stdRes  = std(Residuos{indice_otimo}{pos});

        fatorK = tinv(PA,length(Residuos{indice_otimo}{pos})-1);

        p1 = plot(x,Residuos{indice_otimo}{pos},'.','MarkerSize',15);  
        p2 = plot([x(1) x(end)],[meanRes meanRes],'r--','LineWidth',2);
        p3 = plot([x(1) x(end)],[meanRes+fatorK*stdRes meanRes+fatorK*stdRes],'r-.','LineWidth',2);
        p4 = plot([x(1) x(end)],[meanRes-fatorK*stdRes meanRes-fatorK*stdRes],'r-.','LineWidth',2);

        leg = legend([p1,p2,p3],{'residuos','media','intervalo abrangencia'});

        set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')

        xlabel('Amostra','FontSize',12)
        ylabel('Residuo','FontSize',12)

        set(ax,'FontSize',12)
        
        title(num2str(pos))
        
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/residuos-tendencia.pdf'),resolucao)
        close(fig)
        
        % ========= SERIE ===========

        % BOXPLOT
        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        boxplot(retas{indice_otimo}{pos},'Labels',{strcat('reta ',num2str(pos))})
        ylabel('Quartis');
        
        set(ax,'FontSize',12)
     
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/dados-boxplot.pdf'),resolucao)
        close(fig)
              
        % DISPERSAO
        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        hold(ax,'on')

        x = amostras(pontosInicioAtivos{indice_otimo}(pos):pontosFimAtivos{indice_otimo}(pos));
        meanReta = mean(retas{indice_otimo}{pos});
        stdReta  = std(retas{indice_otimo}{pos});

        fatorK = tinv(PA,length(retas{indice_otimo}{pos})-1);

        p1 = plot(x,retas{indice_otimo}{pos},'.','MarkerSize',15);  
        p2 = plot([x(1) x(end)],[meanReta meanReta],'r--','LineWidth',2);
        p3 = plot([x(1) x(end)],[meanReta+fatorK*stdReta meanReta+fatorK*stdReta],'r-.','LineWidth',2);
        p4 = plot([x(1) x(end)],[meanReta-fatorK*stdReta meanReta-fatorK*stdReta],'r-.','LineWidth',2);

        leg = legend([p1,p2,p3],{'residuos','media','intervalo abrangencia'});

        set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')

        xlabel('Amostra','FontSize',12)
        ylabel('Dados','FontSize',12)

        set(ax,'FontSize',12)

        title(num2str(pos))
        
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/dados-tendencia.pdf'),resolucao)
        close(fig)

        % AUTOCORRELACAO
        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        autocorr(retas{indice_otimo}{pos});

        xlabel('Lag - Amostra','FontSize',12)
        ylabel('Autocorrelacao','FontSize',12)

        set(ax,'FontSize',12)

        title(num2str(pos))
        
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/dados-autocorr.pdf'),resolucao)
        close(fig)

        % HISTOGRAMA
        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        histfit(retas{indice_otimo}{pos},ceil(sqrt(length(retas{indice_otimo}{pos}))));

        xlabel('Serie','FontSize',12)
        ylabel('Densidade de probabilidade','FontSize',12)

        set(ax,'FontSize',12)

        title(num2str(pos))

        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/dados-histfit.pdf'),resolucao)
        close(fig)

        
        % ========= REGIAO ===========

        fig = figure('Visible','off');
        ax = subplot(1,1,1);
        hold(ax,'on')

        Fisher = finv(PA,2,(length(retas{indice_otimo}{pos})-2));
        aspect = FuncaoObjetivo{indice_otimo}{pos}*(2/(length(retas{indice_otimo}{pos})-2)*Fisher);
        pontosEllip_aux = covellipse(parametros{indice_otimo}{pos},Uparametros{indice_otimo}{pos},aspect);

        pontosEllip{pos} = pontosEllip_aux;

        pltellip = plot(pontosEllip_aux(:,1),pontosEllip_aux(:,2),'b-','LineWidth',2);
        plot(parametros{indice_otimo}{pos}(1),parametros{indice_otimo}{pos}(2),'ro','MarkerSize',6)

        pltcoefang = plot([0 0],[min(pontosEllip_aux(:,2)),max(pontosEllip_aux(:,2))],'k:','LineWidth',1.5);
        pltmedia = plot([min(pontosEllip_aux(:,1)) max(pontosEllip_aux(:,1))],[mean(retas{indice_otimo}{pos}),mean(retas{indice_otimo}{pos})],'k-.','LineWidth',2);

        xlabel('Coef. Angular','FontSize',12)
        ylabel('Coef. Linear','FontSize',12)

        set(ax,'FontSize',12)

        leg = legend([pltellip,pltcoefang,pltmedia],{'regiao abrangencia','zero do coef. angular','media dos dados'});

        set(leg,'box','off','Location','SouthOutside','Orientation','horizontal')

        
        title(num2str(pos))
        
        set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
        print(fig,Driver,strcat(folder_retas{pos},'/regiao-abrangencia.pdf'),resolucao)
        close(fig)
    end
end

% ========= REGIAO COMPARACAO ===========
    
fig = figure('Visible','off');
ax = subplot(1,1,1);
hold(ax,'on')
    
for pos = CandidatasEE{indice_otimo}
    if IndiceEstimacao{indice_otimo}{pos}
        pontosEllip_aux = pontosEllip{pos};

        plot(pontosEllip_aux(:,1),pontosEllip_aux(:,2),'LineWidth',2);
        p_aux = parametros{indice_otimo}{pos};
        plot(p_aux(1),p_aux(2),'ro','MarkerSize',6)
    end
end

limites_x = get(ax,'xlim');

xlim([limites_x(1) - (1/10)*(limites_x(2)-limites_x(1)) limites_x(2)+(1/10)*(limites_x(2)-limites_x(1))])

pos_aux = 1; % contador auxiliar para separar a indicacao das retas
for pos = CandidatasEE{indice_otimo}
    p_aux = parametros{indice_otimo}{pos};
    if (-1)^pos_aux == 1 % se pos_aux for par
        x_aux = limites_x(1) - (1/15)*(limites_x(2)-limites_x(1));
    else % caso pos_aux seja impar
        x_aux = limites_x(2) + (1/15)*(limites_x(2)-limites_x(1));
    end
    
    text(x_aux,p_aux(2),strcat('reta  ',num2str(pos)),'FontSize',12)

    plot([x_aux p_aux(1)],[p_aux(2) p_aux(2)],'--')
    pos_aux = pos_aux + 1;
end

limites = get(ax,'ylim');

plot([0 0],[limites(1),limites(2)],'k-.','LineWidth',2);

xlabel('Coef. Angular','FontSize',12)
ylabel('Coef. Linear','FontSize',12)
    
set(ax,'FontSize',12)

set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
print(fig,Driver,strcat('./',projeto,'/','regioes-comparacao.pdf'),resolucao)
close(fig)

% ========= BOXPLOT COMPARACAO ===========
fig = figure('Visible','off');
ax = subplot(1,1,1);

pos_aux = 1;
vetor_retas = [];
group_retas = [];
boxplot_labels = {};
for pos = CandidatasEE{indice_otimo}
    vetor_retas = [vetor_retas;retas{indice_otimo}{pos}];
    group_retas = [group_retas;pos_aux*ones(length(retas{indice_otimo}{pos}),1)];
    
    boxplot_labels{pos_aux} = strcat('reta ',num2str(pos));
    
    pos_aux = pos_aux + 1;
end

boxplot(vetor_retas,group_retas,'Labels',boxplot_labels)
ylabel('Quartis');

set(ax,'FontSize',12)

set(fig,'PaperUnits',unidadeFig,'PaperPosition',PaperPosition,'PaperSize',PaperSize)
print(fig,Driver,strcat('./',projeto,'/','boxplot-comparacao.pdf'),resolucao)
close(fig)


%% Relatorios
% Resumo

folder = strcat('./',projeto);

fileID = fopen(strcat(folder,'/','ResumoGlobal.txt'),'w');

fprintf(fileID,'===================RESUMO GLOBAL ================= \n\n');

fprintf(fileID,'---INFORMACOES GERAIS PROJETO---- \n');
fprintf(fileID,'Numero de pontos da serie: %d \n',length(serie));
fprintf(fileID,'Numero maximo de retas: %d \n',nRetas);
fprintf(fileID,'Probabilidade de Abrangencia: %.2f \n',PA);
fprintf(fileID,'Funcao Objetivo utilizada: %d \n',tipofobj);
fprintf(fileID,'Penalizacao do numero de parametros: %d \n',setN);
fprintf(fileID,'\n');

fprintf(fileID,'---INFORMACOES OTIMIZACAO ---- \n');
fprintf(fileID,'Numero de geracoes: %d \n',options.Generations);
fprintf(fileID,'Tamanho de populacao: %d \n',options.PopulationSize);
fprintf(fileID,'Tolerancia na funcao: %g \n',options.TolFun);
fprintf(fileID,'Tolerancia na restricao: %g \n',options.TolCon);
fprintf(fileID,'Numero de repeticoes da otimizacao: %d \n',length(fval));
fprintf(fileID,'\n');

fprintf(fileID,'---INFORMACOES OTIMIZACAO (RESULTADOS) ---- \n');


fprintf(fileID,'Funcao Objetivo | Pontos de corte ativos: \n');
for pos = 1:length(fval)
    fprintf(fileID,'%d:   %.6g   ->  ',pos,fval(pos));
    fprintf(fileID,'%d ',pontosInicioAtivos{pos});
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');

fprintf(fileID,'Estados Estacionarios identicados (numero da reta) \n');
for pos = 1:length(fval)
    fprintf(fileID,'%d: ',pos);
    fprintf(fileID,'%d ',CandidatasEE{pos});
    fprintf(fileID,'\n');
end
fprintf(fileID,'\n');
fprintf(fileID,'Os relatorios foram salvos para a repeticao numero: %d',indice_otimo);

fclose(fileID);

% Relatorio por reta
for pos = 1:length(pontosInicioAtivos{indice_otimo})
       
    fileID = fopen(strcat(folder_retas{pos},'/','ResumoReta.txt'),'w');
    
    fprintf(fileID,'===========RESUMO DA RETA %d ======= \n',pos);
    fprintf(fileID,'Testes estatisticos para avaliar estacionaridade \n');
    fprintf(fileID,'Intervalo de amostra: %d a %d \n',pontosInicioAtivos{indice_otimo}(pos),pontosFimAtivos{indice_otimo}(pos));
    fprintf(fileID,'Numero de pontos na reta: %d \n',length(retas{indice_otimo}{pos}));
    fprintf(fileID,'Media dos dados: %0.3g \n',mean(retas{indice_otimo}{pos}));    
    fprintf(fileID,'Desvio-padrao dos dados: %0.3g \n',std(retas{indice_otimo}{pos})); 
    fprintf('\n');
    fprintf(fileID,'=============PARAMETROS========== \n');
    fprintf(fileID,'Parametros: Coef. angular | Coef. Linear \n');
    fprintf(fileID,'                %.3g        %.3g     ',parametros{indice_otimo}{pos});
    fprintf(fileID,'\n');
    fprintf(fileID,'Matriz covariancia: \n');
    fprintf(fileID,'[ %.3g  %.3g ] \n[ %.3g  %.3g ]', Uparametros{indice_otimo}{pos});
    fprintf(fileID,'\n');

    fprintf(fileID,'=============RESIDUOS========== \n');
    
    if stat_test.residuo.media(pos)>(1-PA)
        h = 'nao se pode rejeitar a hipotese de que a media seja zero.'; 
    else
        h = 'a media nao e zero.';
    end
    
    fprintf(fileID,' - Media zero(ttest): %0.3g \n',stat_test.residuo.media(pos));
    fprintf(fileID,' - - %s \n',h);
    fprintf('\n');
    fprintf(fileID,'==============SERIE=========== \n');
    
    if stat_test.serie.autocorrLjung(pos)>(1-PA)
        h = 'nao se pode rejeitar a hipotese de que nao ha autocorrelacao na serie'; 
    else
        h = 'ha autocorrelacao na serie.';
    end
    
    fprintf(fileID,' - Autocorrelacao - Ljung-Box Q-test: %.3g \n',...
        stat_test.serie.autocorrLjung(pos));
    fprintf(fileID,' - - %s \n',h);
    
    if stat_test.serie.random(pos)>(1-PA)
        h = 'nao se pode rejeitar a hipotese de que a serie seja aleatoria'; 
    else
        h = 'a serie nao e aleatoria.';
    end
    
    fprintf(fileID,' - Aleatoriedade: %.3g \n',...
        stat_test.serie.random(pos));
    fprintf(fileID,' - - %s \n',h);
    
    if stat_test.serie.normks(pos)>(1-PA)
        h = 'nao se pode rejeitar a hipotese de que a serie seja normal'; 
    else
        h = 'a serie nao e normal.';
    end
        
    fprintf(fileID,' - Normalidade: %.3g \n',...
        stat_test.serie.normks(pos));
    fprintf(fileID,' - - %s \n',h);
    
    if stat_test.serie.normlil(pos)>(1-PA)
        h = 'nao se pode rejeitar a hipotese de que a serie seja normal'; 
    else
        h = 'a serie nao e normal.';
    end
        
    fprintf(fileID,' - Normalidade: %.3g \n',...
        stat_test.serie.normlil(pos));
    fprintf(fileID,' - - %s \n',h);
    
    fprintf(fileID,'\n');
    
    if isempty(find(CandidatasEE{indice_otimo}==pos, 1));
        fprintf(fileID,'Nao e candidata a EE pelo criterio da regiao de abrangencia e autocorrelacao');
    else
        fprintf(fileID,'E candidata a EE pelo criterio da regiao de abrangencia e autocorrelacao');
    end
    
    fclose(fileID);
end
matlabpool close


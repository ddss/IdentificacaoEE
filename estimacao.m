function [ Residuo,dadosReta,parametros,Uparametros ] = estimacao( serie, uyy, pc )
% Função para avaliar as retas: estimação de parâmetros e resíduos
% serie: vetor linha contendo os dados
% uyy: incerteza dos pontos
% pc: pontos de corte ativos

% gerando um vetor para identificar as amostras (posições)
amostras = 1:length(serie);

% convertando o vetor uyy na matrix covariância
Uyy = diag(uyy.^2);

% Criando o vetor que indica as posiçoes do pontos de corte ativos (as extremidades
% sempre estão ativas
pontosAtivos = [1 find(pc==1)+1 length(serie)];

% Avaliar as retas (Regressão de com múltiplos pontos de corte - RLMPC);
Residuo = zeros(1,length(serie)+length(pontosAtivos)-2);

for pos = 1:length(pontosAtivos)-1;
    % obter os dados da reta
    dadosReta  = serie(pontosAtivos(pos):pontosAtivos(pos+1))';
    % obter os dados de x
    xDummy     = [amostras(pontosAtivos(pos):pontosAtivos(pos+1));ones(1,length(dadosReta))]';
    % matriz covariância
    Uyy_aux    = Uyy(pontosAtivos(pos):pontosAtivos(pos+1),pontosAtivos(pos):pontosAtivos(pos+1));
    % estimação dos parâmetros - WLS
    Uparametros = (xDummy'/(Uyy_aux)*xDummy);
    
    parametros   = Uparametros\xDummy'/(Uyy_aux)*dadosReta;
    
    % salvando os resíduos
    if pos == 1
        Residuo(pontosAtivos(pos):pontosAtivos(pos+1)) = (dadosReta - xDummy*parametros);
    else
        Residuo(pontosAtivos(pos)+1:pontosAtivos(pos+1)+1) = (dadosReta - xDummy*parametros);
    end
end

end


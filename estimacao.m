function [ Residuo,dadosReta,parametros ] = estimacao( input_args )
% Função para avaliar as retas: estimação de parâmetros e resíduos


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


function [ Residuo, retas, pontosAtivos, parametros, Uparametros, Residuos, FuncaoObjetivo ] = estimacao( serie, uyy, pc, armazenar )
% Fun????o para avaliar as retas: estima????o de par??metros e res??duos
% serie: vetor linha contendo os dados
% uyy: incerteza dos pontos
% pc: pontos de corte ativos
% armazenar (bool): se retas, parametros e Uparametros devem ser
% armazenados.

% gerando um vetor para identificar as amostras (posi????es)
amostras = 1:length(serie);

% convertando o vetor uyy na matrix covari??ncia
Uyy = diag(uyy.^2);

% Criando o vetor que indica as posi??oes do pontos de corte ativos (as extremidades
% sempre est??o ativas
pontosAtivos = [1 find(pc==1)+1 length(serie)];

% Avaliar as retas (Regress??o de com m??ltiplos pontos de corte - RLMPC);
Residuo = zeros(1,length(serie)+length(pontosAtivos)-2);

% Vetores de armazenamento
retas          = {};
parametros     = {};
Uparametros    = {};
Residuos       = {};
FuncaoObjetivo = {};

for pos = 1:length(pontosAtivos)-1;
    % obter os dados da reta
    dadosReta  = serie(pontosAtivos(pos):pontosAtivos(pos+1))';
    % obter os dados de x
    xDummy     = [amostras(pontosAtivos(pos):pontosAtivos(pos+1));ones(1,length(dadosReta))]';
    % matriz covari??ncia
    Uyy_aux    = Uyy(pontosAtivos(pos):pontosAtivos(pos+1),pontosAtivos(pos):pontosAtivos(pos+1));
    % estima????o dos par??metros - WLS
    invUparametros_reta = (xDummy'/(Uyy_aux)*xDummy);
    
    parametros_reta   = invUparametros_reta\xDummy'/(Uyy_aux)*dadosReta;
    
    res_aux = (dadosReta - xDummy*parametros_reta);
    % salvando os res??duos
    if pos == 1
        Residuo(pontosAtivos(pos):pontosAtivos(pos+1)) = res_aux;
    else
        Residuo(pontosAtivos(pos)+1:pontosAtivos(pos+1)+1) = res_aux;
    end
    
    % armazenar as informa????es
    if armazenar
        retas{pos}       = dadosReta;
        parametros{pos}  = parametros_reta;
        Uparametros{pos} = inv(invUparametros_reta);
        FuncaoObjetivo{pos} = res_aux'/Uyy_aux*res_aux;
        Residuos{pos}    = (dadosReta - xDummy*parametros_reta);
    end
end

end


function [ Residuo, NE, retas, pontosInicioAtivos, pontosFimAtivos, parametros, Uparametros, Residuos, FuncaoObjetivo, CandidatasEE,IndiceEstimacao, phi ] = estimacao( serie, uyy, pc, PA, NEprojeto, armazenar )
% Funcaoo para avaliar as retas: estimacao de parametros e residuos
% serie: vetor linha contendo os dados
% uyy: incerteza dos pontos
% pc: pontos de corte ativos - marcam o inicio das retas
% armazenar (bool): se retas, parametros e Uparametros devem ser
% armazenados.
% phi - penalizacaoo da funcao objetivo com base no numero de pontos
% NEprojeto - numero de pontos minimo para ser EE

% gerando um vetor para identificar as amostras (posicoes)
amostras = 1:length(serie);

% convertando o vetor uyy na matrix covari??ncia
Uyy = diag(uyy.^2);

% Criando o vetor que contem as posicoes de inicio das retas
pontosInicioAtivos = [1 sort(pc(1:pc(end)))]

% Criando o vetor que contem as posicos de fim das retas
pontosFimAtivos    = [sort(pc(1:pc(end)))-1 length(serie)]

% Avaliar as retas (Regressaoo de com multiplos pontos de corte - RLMPC);
Residuo = zeros(1,length(serie));

% Vetores de armazenamento
retas          = {};
parametros     = {};
Uparametros    = {};
Residuos       = {};
FuncaoObjetivo = {};
IndiceEstimacao = {}; % Armazenar se a estimacao foi executada com sucesso

posCandidatasEE = 1;
NE              = 0;
CandidatasEE    = [];

phi = 0;

var_a = ones(1,length(pontosInicioAtivos));

for pos = 1:length(pontosInicioAtivos);
    % obter os dados da reta
    dadosReta  = serie(pontosInicioAtivos(pos):pontosFimAtivos(pos))'
    size(dadosReta)
    % obter os dados de x
    xDummy     = [1:length(dadosReta);ones(1,length(dadosReta))]'
    size(xDummy)
    % matriz covariancia
    Uyy_aux    = Uyy(pontosInicioAtivos(pos):pontosFimAtivos(pos),pontosInicioAtivos(pos):pontosFimAtivos(pos))
    size(Uyy_aux)
    if length(dadosReta)~=1
        % estimacao dos parametros - WLS
        invUparametros_reta = (xDummy'/(Uyy_aux)*xDummy);

        parametros_reta   = invUparametros_reta\xDummy'/(Uyy_aux)*dadosReta;

        res_aux = (dadosReta - xDummy*parametros_reta);
    else
        % estimacao dos parametros - WLS
        invUparametros_reta = [1 0;0 1];

        parametros_reta   = [NaN;NaN];

        res_aux = 1e100;        
    end
    % salvando os residuos
    Residuo(pontosInicioAtivos(pos):pontosFimAtivos(pos)) = res_aux;
    
    % armazenar as informacoes
    
    fobj = res_aux'/Uyy_aux*res_aux;
    
    Uparametros_reta = inv(invUparametros_reta);
    
    var_a(pos) = Uparametros_reta(1,1);
    Fisher(pos) = finv(PA,2,(length(dadosReta)-2));
    aspect(pos) = fobj*(2/(length(dadosReta)-2)*Fisher(pos));
    
    retas{pos}       = dadosReta;
    parametros{pos}  = parametros_reta;    
    Uparametros{pos} = Uparametros_reta;
    
    if armazenar
        FuncaoObjetivo{pos} = fobj;
        Residuos{pos}    = (dadosReta - xDummy*parametros_reta);
        if length(dadosReta) == 1
           IndiceEstimacao{pos} = false;
        else
           IndiceEstimacao{pos} = true;
        end
    end
end

for pos = 1:length(pontosInicioAtivos);
    %% Calculo da regiao de abrangencia e verificacao das retas candidatas
    % Considerando que os parametros seguem uma distribuicao normal
       
    desvio = retas{pos} - mean(retas{pos});
    
    % Teste para avlaliar a autocorrelacao
    if length(desvio)>max(2,ceil(0.3*NEprojeto))
        [~,teste_autocorr] = lbqtest(desvio);
    else
        teste_autocorr = 0;
    end
    
    
    %if and(([0;mean(retas{pos})]-parametros{pos})'/Uparametros{pos}*([0;mean(retas{pos})]-parametros{pos}) <= aspect(pos),var_a(pos)<mean(var_a))
    if and(([0;mean(retas{pos})]-parametros{pos})'/Uparametros{pos}*([0;mean(retas{pos})]-parametros{pos}) <= aspect(pos),teste_autocorr>=(1-PA))
        CandidatasEE(posCandidatasEE) = pos;
        posCandidatasEE = posCandidatasEE+1;
        NE = NE + length(retas{pos});
        
        phi = phi + length(retas{pos})/NEprojeto;
    end

end

end


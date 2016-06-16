function [ Residuo, NE, retas, pontosAtivos, parametros, Uparametros, Residuos, FuncaoObjetivo, CandidatasEE ] = estimacao( serie, uyy, pc, PA, armazenar )
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

posCandidatasEE = 1;
NE              = 0;
CandidatasEE    = [];
for pos = 1:length(pontosAtivos)-1;
    % obter os dados da reta
    dadosReta  = serie(pontosAtivos(pos):pontosAtivos(pos+1))';
    % obter os dados de x
    xDummy     = [1:length(dadosReta);ones(1,length(dadosReta))]';
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
    
    fobj = res_aux'/Uyy_aux*res_aux;
    
    if armazenar
        retas{pos}       = dadosReta;
        parametros{pos}  = parametros_reta;
        Uparametros{pos} = inv(invUparametros_reta);
        FuncaoObjetivo{pos} = fobj;
        Residuos{pos}    = (dadosReta - xDummy*parametros_reta);
    end
    %% C?lculo da regi?o de abrang?ncia e verifica??o das retas candidatas
    % Aqui est? apenas se calculando os pontos extremos da elipse
    % Considerando que os par?metros seguem uma distribui??o normal

    Fisher = finv(PA,2,(length(dadosReta)-2));
    aspect = fobj*(2/(length(dadosReta)-2)*Fisher);

    fator = invUparametros_reta(1,1)/(invUparametros_reta(1,2) + eps); % eps evita NaN quando a covari?ncia ? zero.
    delta = sqrt(aspect/(fator^2*invUparametros_reta(2,2)-2*fator*invUparametros_reta(1,2)+invUparametros_reta(1,1)));
    coordenadas_x = [parametros_reta(1)+delta      ,parametros_reta(1)-delta];
    coordenadas_y = [parametros_reta(2)-delta*fator,parametros_reta(2)+delta*fator];

    fator = invUparametros_reta(2,2)/(invUparametros_reta(1,2) + eps); % eps evita NaN quando a covari?ncia ? zero.
    delta = sqrt(aspect/(fator^2*invUparametros_reta(1,1)-2*fator*invUparametros_reta(1,2)+invUparametros_reta(2,2)));
    coordenadas_y = [coordenadas_y [parametros_reta(2)+delta      ,parametros_reta(2)-delta]];
    coordenadas_x = [coordenadas_x [parametros_reta(1)-delta*fator,parametros_reta(1)+delta*fator]];

    % Obten??o das posi??es das retas candidatas a EE
    % - verificar se e ellipse do par?metro a, cruza o zero.
    if and(and(any(coordenadas_x>0), any(coordenadas_x<0)),and(any(coordenadas_y>mean(dadosReta)), any(coordenadas_y<mean(dadosReta))))
        CandidatasEE(posCandidatasEE) = pos;
        posCandidatasEE = posCandidatasEE+1;
        NE = NE + length(dadosReta);
        
    end
end

end


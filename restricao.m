function [ c,ceq ] = restricao( pc,nRetas )
% Restri??o para a otimiza??o: evitar a solu??o particular
% pc    : vetor contendo os pontos de corte
% nRetas: n?mero m?ximo poss?vel de retas 

ceq = [];

if sum(pc) == 0
    
    c = 0;

else
    
    c = log((sum(pc)+1)/nRetas)*10^10;

end

end

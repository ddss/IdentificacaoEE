function [ c,ceq ] = restricao( pc,nRetas )
% Restri??o para a otimiza??o: evitar a solu??o particular
% pc    : vetor contendo os pontos de corte
% nRetas: n?mero m?ximo poss?vel de retas 

ceq = [];

c = nRetas-1;
for cont = 1:nRetas-1
     
    if pc(cont)<pc(cont+1)
       c = c-1; 
    end
        
end

end

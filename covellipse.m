function [ pontosEllipse ] = covellipse( centro, cov, aspectratio )
% Avalia a regi?o de abrang?ncia ellipsoidal
% ENTRADAS:
% centro = vetor contendo o centro da ellipse
% cov = matriz de covari?ncia
% aspectratio = aspecto da ellipse
% SAIDA:
% vetor contendo pontos pertencentes ? ellipse
% REFERENCIAS
% fonte: http://trunghuyduong.blogspot.com.br/2010/10/contents-input-parameters-of-2-d.html
% fonte: http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

% Avaliacao dos autovalores e autovetores
[auto_vetores, auto_valores] = eig(cov);

% comprimento dos semieixos
eixo_maior = sqrt(auto_valores(2,2)*aspectratio);
eixo_menor = sqrt(auto_valores(1,1)*aspectratio);

% angulo de rota??o de elipse em rela??o ao x
angulo_ellip = atan2(auto_vetores(1,2),auto_vetores(2,2));

if angulo_ellip < 0
    angulo_ellip = angulo_ellip + 2*pi;
end

% valores de angulo para gerar os pontos da ellipse
angulo_avaliar = 0:pi/100:2*pi;

% Avaliacao da ellipse sem rotacao
pontosEllipse(1,:) = eixo_menor*sin(angulo_avaliar);
pontosEllipse(2,:) = eixo_maior*cos(angulo_avaliar);

% Rotacao 
Q = [cos(angulo_ellip)  sin(angulo_ellip)
     -sin(angulo_ellip) cos(angulo_ellip)];

pontosEllipse = Q*pontosEllipse;

pontosEllipse = pontosEllipse';

% translacao
pontosEllipse(:,1) = pontosEllipse(:,1) + centro(1);
pontosEllipse(:,2) = pontosEllipse(:,2) + centro(2);

end


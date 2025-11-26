function [Kel] = matK_elem(S1, S2, S3,Reftri)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) Utilisation d une quadrature a 3 point d ordre 2
%    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end


% calcul de la matrice de raideur
% -------------------------------

if Reftri==1
    sigma=@sigma_1;
elseif Reftri==2
    sigma=@sigma_2;
    %On a besoin de savoir l'indice du triangle mais on ne connait que les coordonnées,
    %donc on pourrait remonter aux sommets puis au triangle associé mais ça
    %serait long et deplus on a pas accès aux tableaux Coordneu etc dans la
    %fonction matK_elem.m
end

B=[x2-x1 x3-x1;y2-y1 y3-y1];
gphi=[-1 -1;1 0;0 1];
det=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
BTinv=[y3-y1 y1-y2;x1-x3 x2-x1]/det; %l'inversion et la transpoée sont appliquées
A1=B*[1/6;1/6]+[x1;y1];
A2=B*[2/3;1/6]+[x1;y1];
A3=B*[1/6;2/3]+[x1;y1]; % quadrature à 3 points de Gauss Legendre

Kel = zeros(3,3);
for i=1:3
  for j=1:3
	% A COMPLETER
    Kel(i,j) = (sigma(A1(1),A1(2))+sigma(A2(1),A2(2))+sigma(A3(1),A3(2)))*dot(BTinv*transpose(gphi(i,:)),BTinv*transpose(gphi(j,:)))*abs(det)/6;
  end % j
end % i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2025

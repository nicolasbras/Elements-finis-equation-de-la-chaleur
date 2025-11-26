% =====================================================
%
% principal_chaleur;
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour 
% 1) l'equation de la chaleur suivante stationnaire, avec condition de
% Dirichlet non homogene
%
% | \alpha T - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2
% |         T = T_\Gamma,   sur le bord
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
%
% 2) l'equation de la chaleur dependant du temps avec condition de 
% Dirichlet non homogene
%
% | dT/dt - div(\sigma \grad T)= S,   dans \Omega=\Omega_1 U \Omega_2 et pour tout t< t_max
% |         T = T_\Gamma,   sur le bord et pour tout t< t_max
% |         T = T_0       dans \Omega et pour t=0  
%
% ou S est la source de chaleur, T_\Gamma la temperature exterieure,
% T_0 est la valeur initiale de la temp?rature
% \alpha > 0 et
% \sigma = | \sigma_1 dans \Omega_1
%          | \sigma_2 dans \Omega_2
% =====================================================
% Donnees du probleme
% ---------------------------------
h=0.005;
system(['gmsh -2 -clmax ' num2str(h) ' -clmin ' num2str(h) ' geomChaleur.geo']);
nom_maillage = 'domaine.msh' ;

validation = 'non';
pb_stationnaire = 'non';
pb_temporel = 'oui';

if strcmp(validation,'oui')
    alpha = 1;
    T_Gamma = 0;
end

if strcmp(pb_stationnaire,'oui')
    alpha = 1;
    T_Gamma = 285;
end

if strcmp(pb_temporel,'oui')
    Tps_initial = 0;
    Tps_final = 1;
    delta_t = 0.01;
    alpha = 1/delta_t;
    N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d'iterations necessaires
    T_Gamma = 0;
end

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de masse
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  
  % calcul des matrices elementaires du triangle l 
  
   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:),Reftri(l));
   % LA ROUTINE matK_elem.m DOIT ETRE MODIFIEE

   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));
    
    % On fait l'assemblage de la matrice globale
    % A COMPLETER
   for i=1:3
       for j=1:3
           KK(Numtri(l,i),Numtri(l,j))=KK(Numtri(l,i),Numtri(l,j))+Kel(i,j);
           MM(Numtri(l,i),Numtri(l,j))=MM(Numtri(l,i),Numtri(l,j))+Mel(i,j);
       end
   end
end % for l

% Matrice EF
% -------------------------
AA = alpha*MM+KK;
% =====================================================
% =====================================================
% Pour le probleme stationnaire et la validation
% ---------------------------------

% Calcul du second membre F
% -------------------------
% A COMPLETER EN UTILISANT LA ROUTINE f.m
FF=zeros(Nbpt,1);
for i=1:Nbpt
     FF(i)=f2(Coorneu(i,1),Coorneu(i,2));
end
for i=1:Nbpt
    if Refneu(i)==1
        FF(i)=T_Gamma_fonction(Coorneu(i,1),Coorneu(i,2));
    end
end

LL = MM*FF(:);

% inversion
% ----------
% tilde_AA ET tilde_LL SONT LA MATRICE EF ET LE VECTEUR SECOND MEMBRE
% APRES PSEUDO_ELIMINATION 
% ECRIRE LA ROUTINE elimine.m ET INSERER L APPEL A CETTE ROUTINE
% A UN ENDROIT APPROPRIE
[tilde_AA,tilde_LL]=elimine(AA,LL,Refneu,Coorneu);

UU = tilde_AA\tilde_LL;

TT = UU-T_Gamma;
% validation
% ----------
if strcmp(validation,'oui')
    UU_exact=zeros(Nbpt,1);
    for i=1:Nbpt
    UU_exact(i) = sin(3*pi*Coorneu(i,1))*sin(pi*Coorneu(i,2));
    end
	% Calcul de l erreur L2
	% A COMPLETER
    erreurL2=sqrt(transpose(UU_exact-UU)*MM*(UU_exact-UU))
	% Calcul de l erreur H1
	% A COMPLETER
    erreurH1=sqrt(transpose(UU_exact-UU)*KK*(UU_exact-UU))
	% attention de bien changer le terme source (dans FF)
end

% =====================================================
% =====================================================
% Pour le probleme temporel
% ---------------------------------
if strcmp(pb_temporel,'oui')

    % on initialise la condition initiale
    % -----------------------------------
    T_initial = condition_initiale(Coorneu(:,1),Coorneu(:,2));

	% solution a t=0
	% --------------
    UU = 310*ones(Nbpt,1);
    TT = 310*ones(Nbpt,1);

	% Boucle sur les pas de temps
	% ---------------------------
    for k = 1:N_t
        LL_k = zeros(Nbpt,1);
        S=zeros(Nbpt,1);
        for i = 1:Nbpt
            x=Coorneu(i,1);
            y=Coorneu(i,2);
            S(i,1) =  f_t(x,y,k*delta_t);
        end
        
        % Calcul du second membre F a l instant k*delta t
        % -----------------------------------------------
		% A COMPLETER EN UTILISANT LA ROUTINE f_t.m et le ou les termes precedents (donne par UU)
		LL_k = MM*UU/delta_t+MM*S;
        AA_k=MM/delta_t+KK;

		% inversion
		% ----------
		% Calculer la matrice intervenant dans le schema, 
		% faire la pseudo elimintaion si nécessaire 
        [tilde_AA_k,tilde_LL_k]=elimine(AA_k,LL_k,Refneu,Coorneu);
        UU = tilde_AA_k\tilde_LL_k;
        TT = UU;

        % visualisation 
        pause(0.05)
        affiche(TT, Numtri, Coorneu, ['Temps = ', num2str(k*delta_t)]);
        axis([min(Coorneu(:,1)),max(Coorneu(:,1)),min(Coorneu(:,2)),max(Coorneu(:,2)),...
            ]);
        caxis([100,320]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%2024%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%2025

% visualisation
% -------------
if ( strcmp(validation,'oui') || strcmp(pb_stationnaire,'oui') )
    affiche(UU, Numtri, Coorneu, sprintf('Dirichlet - %s', nom_maillage));
end
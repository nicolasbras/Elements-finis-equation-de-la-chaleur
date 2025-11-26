function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu,Coorneu)
nbpt=length(Refneu);
tilde_AA=AA;
tilde_LL=zeros(nbpt,1);

haut=find(Refneu~=1);
bas=find(Refneu==1);
ordre=[transpose(haut),transpose(bas)];

lim = length(haut)+1;

for i=1:length(bas)
    tilde_LL(bas(i))=T_Gamma_fonction(Coorneu(i,1),Coorneu(i,2));
    tilde_AA(bas(i),:)=0;
    tilde_AA(:,bas(i))=0;
    tilde_AA(bas(i),bas(i))=1;
end
for i=1:length(haut)
    tilde_LL(haut(i))=LL(haut(i));
end
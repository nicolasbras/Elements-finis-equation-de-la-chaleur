function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu,Coorneu,MM)
nbpt=length(Refneu);
tilde_AA=AA;
tilde_LL=zeros(nbpt,1);

haut=find(Refneu~=1);
bas=find(Refneu==1);
ordre=[transpose(haut),transpose(bas)];

%tilde_AA=AA(ordre,ordre);
%tilde_LL=LL(ordre);

lim = length(haut)+1;
%for i=lim:nbpt
    %tilde_LL(i)=0;
    %tilde_AA(i,:)=0;
    %tilde_AA(:,i)=0;
    %tilde_AA(i,i)=1;
%end

for i=1:length(bas)
    tilde_LL(bas(i))=T_Gamma(Coorneu(bas(i),1),Coorneu(bas(i),2));
    tilde_AA(bas(i),:)=0;
    tilde_AA(:,bas(i))=0;
    tilde_AA(bas(i),bas(i))=1;
end
tilde_LL=MM*tilde_LL(:);
for i=1:length(haut)
    tilde_LL(haut(i))=LL(haut(i));
end
tilde_LL
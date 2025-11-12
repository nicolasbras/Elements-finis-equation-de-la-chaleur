function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu)
nbpt=length(Refneu);
tilde_AA=AA;
tilde_LL=LL;
lim=1;
for i=1:nbpt
    if Refneu(i)==0 && i<nbpt
        for j=i+1:nbpt
            if Refneu(j)==1
                Refneu([i j])=Refneu([j i]);
                tilde_AA([i j],:)=tilde_AA([j i],:);
                tilde_AA(:,[i j])=tilde_AA(:,[j i]);
                tilde_LL([i j])=tilde_LL([j i]);
            end
        end
    end
end
for i=1:nbpt-1
    if Refneu(i)~=0
        lim=lim+1;
    end
end
%for i=lim:nbpt
 %   tilde_LL(i)=0;
    for j=1:nbpt
        tilde_AA(i,j)=0;
        tilde_AA(j,i)=0;
    end
    tilde_AA(i,i)=1;
end
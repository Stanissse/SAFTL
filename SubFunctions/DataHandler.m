function [Dat_TOT,Dat_UAF]=DataHandler(j)
% function splits up data in UAF and TOT -> image per image (index j)

global pxSize DAT middle
 
N_Images = max(DAT(:,1));

for i=1:max(N_Images)
    Num_Spots_Im_(i)= length(DAT(DAT(:,1)==i));
    Num_Spots_Im(i) = sum(Num_Spots_Im_(1:i));
end
    

% Daten Aufteilen

clear DAT_temp
 
if j ==1
    DAT_temp=DAT(1:Num_Spots_Im(j),:);
else
    DAT_temp=DAT(Num_Spots_Im(j-1)+1:Num_Spots_Im(j),:);
end
    
DAT_temp = sortrows(DAT_temp,3);
l_tot=length(DAT_temp(DAT_temp(:,3)<middle*pxSize));

Dat_TOT=DAT_temp(1:l_tot,:);
[d_DAT, ~]=size(DAT_temp);
[d_tot, ~]=size(Dat_TOT);


% trow exception if there aere no spots in UAF
if d_tot<d_DAT
    Dat_UAF=DAT_temp(l_tot+1:d_DAT,:);
elseif d_tot==0
    Dat_UAF = ones(1,5);
    Dat_TOT = 1e+5*ones(1,5);
else
    Dat_UAF = ones(1,5);
end



end
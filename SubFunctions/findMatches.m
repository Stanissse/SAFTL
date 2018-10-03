function [Ind_TOT,Ind_UAF]=findMatches(Dat_UAF,Dat_TOT)
global pxSize dist_max
% Dat_UAF = Data from UAF image
% Dat_TOT = Transformed Data from TOT image
% dist_max = max distance between matching points

Ind_TOT=[];
Ind_UAF=[];

dist = []; 


    function [Ind_TOT,Ind_UAF]=findMatches_(Dat_UAF,Dat_TOT)
    %calculate distances between all points
    [l_tot, ~]=size(Dat_TOT);
    [l_uaf, ~]=size(Dat_UAF);
    
    for i=1:l_tot
        for j=1:l_uaf
            dist(i,j)= norm(Dat_UAF(j,2:3)-Dat_TOT(i,2:3));
        end
    end
    % dist = max distance between points
    % finds matching points
    dist(dist(:)<dist_max)=0;
    [Ind_TOT,Ind_UAF] = find(dist==0);
    end

[Ind_TOT,Ind_UAF]=findMatches_(Dat_UAF,Dat_TOT);

% if length(Ind_UAF)>length(Dat_UAF) fasle matched points are found
% ----> rduce maximal allowed distance between matching points
while(length(Ind_UAF)>length(Dat_UAF))
    dist_max = dist_max-10;
    
    text=['Reduncing maximal allowed distance between matching points to dist_max = ' , num2str(dist_max), ' nm'];
    disp(text)
    [Ind_TOT,Ind_UAF]=findMatches_(Dat_UAF,Dat_TOT);
end

end
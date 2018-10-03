function [out] = deleteZeros(In)

X_= In;
%Y_= In(:,2);

% searchs positions of unphysical ratios
Delete = find(In==0);

% sets unphysical ratios to zero
X_(Delete) = 0;
%Y_(Delete) = 0;

% delete unphysical ratios
X = X_(X_~=0);
%Y = Y_(Y_~=0);

out(:,1)=X;
%out(:,2)=Y;

end
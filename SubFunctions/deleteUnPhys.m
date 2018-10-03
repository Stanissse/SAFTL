function [X, Y, Z] = deleteUnPhys(X_, Y_, Z_)

% searchs positions of unphysical ratios
Delete = find(Z_>=200 | Z_ <= 50);

% sets unphysical ratios to zero
X_(Delete) = 0;
Y_(Delete) = 0;
Z_(Delete) = 0;

% delete unphysical ratios
X = X_(X_~=0);
Y = Y_(Y_~=0);
Z = Z_(Z_~=0);

end
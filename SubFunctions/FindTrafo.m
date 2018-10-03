function [t_concord]=FindTrafo(imagedata)

% Punkte suchen und transformation bestimmen
[aerial_points,ortho_points] = ...
       cpselect(imagedata,imagedata,...
                'Wait',true);

% Transformation
t_concord = fitgeotrans(aerial_points,ortho_points,'projective');

end
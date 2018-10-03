function [Kx,Ky,Kr,pupil,defocus,N_pad,uk,ux,v]=fun_usualsuspects(dia_pupil,lambda_0,ux_desired,NA,RI)
%[Kx,Ky,Kr,pupil,defocus,N_pad,uk,ux,v]=fun_usualsuspects(dia_pupil,lambda_0,ux_desired,NA,RI)
%this function returns the "usually needed variables" of my simulations,
%e.g. a 2D k-space grid for defining an objective pupil, the pupil itself,
%the defocus phase for 1m of spherical defocus, etc..
%inputs:
%lambda_0...vacuum wavelength
%ux_desired....desired unit in x-space (this can very closely, but not
%truly be fulfilled because the grid has to be of integer size)
%outputs:
%Kx,Ky,Kr...k-space grid (in units of k=2*pi*RI/lambda_0)
%pupil...circular mask representing the objective pupil (its 2D projection
%on Kx,Ky)
%N_pad...size to which the pupil function must be padded in order to
%achieve a resolution of ux in x-space (use "embed(E,N,0)" to pad with
%zeros)
%defocus...phase function representing 1m of defocus in x-space
%uk...k-space unit (in units of 2*pi*RI/lambda_0)
%ux...true x-space unit
%v..."zoomed" view axis in x-space; use for instance E(v,v) to see E zoomed
%in around the focal spot

k=2*pi*RI/lambda_0;
uk=2*k*NA/RI/dia_pupil; %unit in k-space
[Kx,Ky,Kr,pupil]=create_coord(dia_pupil,uk,'FFT');
N_pad=round(2*pi/ux_desired/uk); %size of padded field (to achieve targeted resolution in object space)
ux=2*pi/N_pad/uk; %actual resolution in object space
defocus=(sqrt(k^2-Kr.^2)).*pupil;  %if defocus is convex, correct focus stack can be calculated by adding phase factors exp(1i*dz*defocus) before F-Trafo
defocus=defocus-mean2(defocus(pupil)); %defocus function for 1 meter of defocus
ux=2*pi/N_pad/uk; %real resolution in object space (the larger the grid, the closer it will be to ux_desired)

dia_airy=(1.22*lambda_0/NA)/ux; %airy diameter in pixels
win=round(3*dia_airy);
v=round((N_pad+1)/2)-win:round((N_pad+1)/2)+win;



function [SA_out,Defocus,a] =fun_SA_RImismatch(N,n1,n2,NA,lambda_0,rd)
%outputs spherical aberration arising from a RI-mismatch, for z=1m inside
%the material of RI2. Output can be scaled to the existing z-value
%returned defocus phase is defined for medium n2.
%according to paper Jesacher et al. 2010. Opex
%the "defocus-term" is removed analytically
%inputs:
%n1...refractive index of immersion medium
%n2...refractive index of target material
%NA...numerical aperture (e.g. 1.4)
%lambda_0...vacuum wavelength
%a...coefficient of "Defocus" contained in SA
%rd...0 or 1, if set to 1, defocus is removed

%% ---test parameters---
% lambda_0=532e-9;
% k=2*pi/lambda_0;
% NA=1.4;
% n2=1.6; %RI of target material
% n1=1.52; %RI of immersion medium
% N=256;
% rd=1;

k=2*pi/lambda_0;
[tmp,tmp,R,mask]=create_coord(N,2/(N-1),'FFT');

%mean value of spherical defocus
MW=@(RI) 2/3*k*(-(-1 + RI^2/NA^2)^(3/2)+(RI^2/NA^2)^(3/2))*NA;
%defocus function
Def=@(RI) real(k*NA*sqrt(RI^2/NA^2-R.^2)-MW(RI)).*mask; %spherical defocus in medium with refractive index RI (for 1m of defocus)
%Def2=real(k*NA*sqrt(n2^2/NA^2-R.^2)-MW2).*mask; %spherical defocus in medium 2 (for 1m of defocus)

SA=-(Def(n2)-Def(n1)); %spherical aberration phase (complete, i.e including defocus) (for 1m of defocus), (see paper: Jesacher & Booth, Opex 2010)
Defocus=Def(n2);
%the following analytic results were calculated with Mathematica: see file
%"aberrations of a refractive index mismatch.nb"
%a_spher=2*pi*integrate((Def2-MW2)*SA*R,{R,0,1} 


a_spher=-1/72*k^2*pi*(72*n2^2-36*NA^2+(32*(-(-1 + n2^2/NA^2)^(3/2) + n2^3/NA^3)*(n1^3-n1^2*sqrt((n1 - NA)*(n1 + NA))+NA^2*sqrt((n1 - NA)*(n1 + NA))))/NA - (32*(n2^3 - n2^2*sqrt((n2 - NA)*(n2 + NA))+NA^2*sqrt((n2 - NA)*(n2 + NA)))^2)/NA^4+(9/NA^2)*(2*(-n1^3*n2 - n1*n2^3 + n1^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))+sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))*(n2^2-2*NA^2)) - (n1^2 - n2^2)^2*(log((n1 - n2)^2)-log(n1^2 + n2^2 - 2*(NA^2 + sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))))));
%a_spher=(1/(72*NA^4))*k^2*pi*(2*(32*n2^6 + 9*n1*n2^3*NA^2-48*n2^4*NA^2-32*n2^5*sqrt((n2 - NA)*(n2 + NA))+n1^3*(-16*n2^3+9*n2*NA^2+16*n2^2*sqrt((n2 - NA)*(n2 + NA))-16*NA^2*sqrt((n2 - NA)*(n2 + NA)))+2*NA^4*(NA^2+sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))+n2^2*NA^2*(12*NA^2+7*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))+n1^2*(16*n2^3*sqrt((n1 - NA)*(n1 + NA))-16*n2^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))+7*NA^2*sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA)))-16*n2^3*NA^2*(-2*sqrt((n2 - NA)*(n2 + NA))+sqrt(n1^2 - NA^2)))+9*(n1^2 - n2^2)^2*NA^2*(log((n1 - n2)^2)-log(n1^2+n2^2-2*(NA^2 + sqrt((n1 - NA)*(n2 - NA)*(n1 + NA)*(n2 + NA))))));

%def_norm=2*pi*integrate((Def2-MW2)^2*R,{R,0,1}
def_norm=-((k^2*(16*n2^6 - 24*n2^4*NA^2 + 6*n2^2*NA^4 + NA^6 - 16*n2^5*sqrt(n2^2 - NA^2) + 16*n2^3*NA^2*sqrt(n2^2 - NA^2))*pi)/(18*NA^4));
SA_corr=SA-a_spher*Def(n2)/def_norm;
a=a_spher/def_norm; %coefficient of Defocus contained in SA (without normalization)

%imagesc(abs(SA_corr));

%comparison with numerical solution
%     [a_SD, SD, ~]=spherical_defocus_analysis(SA,NA,n1,n2);
%     SA_corr_num=SA-a_SD*SD; %numerical check
%     plot(SA_corr(end/2,:));
%     hold on;
%     plot(SA_corr_num(end/2,:),'r'); %numerical check
%     hold off;
% 
if rd==1;
    SA_out=SA_corr;
else 
    SA_out=SA;
end


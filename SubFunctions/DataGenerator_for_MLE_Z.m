%calculating CRB for SAF-detection and finding optimal PSF for SAF-detection by minimizing the CRB:
%-------------------------------------------------------------
%program calculates BFP-fields for various distances (z) of the dipole from the coverslip
%from these data the CRB for x-y-z expectation values are calculated
tic
clear all;
% close all;

global Z ux uk uz Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz Nx n_photon bg mask

%% calculating BFP fields for SAF

N=128;
lambda_0=680e-9;

defocus=-300e-9;  %defocus of objective lens 

NA=1.67; RI=[1.45 1.45 1.78]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]

d2=0e-9; %thickness of intermediate layer (layer 2)
f=1.8e-3; %focal length of objective
mu=1e-14; %magnitude of dipole (arbitrary value)

ux=115e-9; %resolution in focal space
Nx_temp=50;
Nx=14; %desired simulated field size in pixel
%-----------------------
[SA_out,Defocus,~] =fun_SA_RImismatch(N,RI(3),RI(3),NA,lambda_0,1); %Defocus function refers to refractive index n2
uk=4*pi/lambda_0*NA/N; %unit in pupil space (k-space)
[~,~,R,pupil]=create_coord(N,1,'FFT');

% polaristion dependendent transmission
load('ZernikeCoeff_full.mat'); 
Im_Fit = sum(ZernikeCalc([1 4 5 6 11 12 13 22 37],ZernCoeff,length(pupil),'NOLL'),3);
Im_Fit = Im_Fit/max(Im_Fit(:));
obj_transm=sqrt(Im_Fit);
clear ZernCoeff; 

% obj_transm=polyval([0.0832 0 -0.5199 0 1],asin(2/N*R*NA/RI(3))); %objective intensity transmission function model (from measurement with blue textmarker)
pupil_UAF=circshift(R<=((N/2)*(RI(1)/NA))*1,[0 0]); %pupil containing UAF light
aberrNum = [2 : 37,56];
        
load('coeff_Example.mat');
Zernike_stack=ZernikeCalc(aberrNum,ones(length(Z_aberr),1),pupil,'Noll');
hugo=zeros(1,1,size(Zernike_stack,3));
hugo(1,1,1:length(Z_aberr2))=[0 0 Z_aberr2(3:end)];%zeros(length(Z_aberr2),1);%
aberr=sum(Zernike_stack.*repmat(hugo,[N,N,1]),3);

uz=.005e-6;%.025e-6;
z_vec=(0:uz:500e-9); %simulated dipole distances above layer 2 z-dipole
mask=obj_transm.*exp(1i*aberr+1i*defocus*Defocus);       
        
for dz=1:length(z_vec); %BFP image number        %calculating BFP-fields for all dipole orientations

    % dipole=[0,0] [theta, phi], e.g. [0,0] for z-dipole, [pi/2,0] for x-dipole
    [Ex_Pz,Ey_Pz]=fun_dipole_imaging(N,lambda_0,NA,RI,[0,0],d2,z_vec(dz),f,mu); %z-dipole
    [Ex_Px,Ey_Px]=fun_dipole_imaging(N,lambda_0,NA,RI,[pi/2, 0],d2,z_vec(dz),f,mu); %x-dipole
    [Ex_Py,Ey_Py]=fun_dipole_imaging(N,lambda_0,NA,RI,[pi/2, pi/2],d2,z_vec(dz),f,mu); %y-dipole

    %-----calculating total (SAF+UAF) images-----
    I_xx=abs(czt2(Ex_Px.*mask,uk,ux,Nx_temp)).^2;
    I_yx=abs(czt2(Ey_Px.*mask,uk,ux,Nx_temp)).^2;
    I_xy=abs(czt2(Ex_Py.*mask,uk,ux,Nx_temp)).^2;
    I_yy=abs(czt2(Ey_Py.*mask,uk,ux,Nx_temp)).^2;    
    I_xz=abs(czt2(Ex_Pz.*mask,uk,ux,Nx_temp)).^2;
    I_yz=abs(czt2(Ey_Pz.*mask,uk,ux,Nx_temp)).^2;
    PSF_tot_(:,:,dz)=I_xx+I_yx+I_xy+I_yy+I_xz+I_yz;
    C=PSF_tot_(:,:,dz);
    PSF_tot_(:,:,dz)=PSF_tot_(:,:,dz)/sum(sum(C(:)));
    
%     PSF_tot_(:,:,dz)=PSF_tot_(:,:,dz)/(max(C(:)));
end    

PSF_tot(:,:,:)=PSF_tot_(25-6:25+7,25-6:25+7,:);

save('PSF_SAF_5nm.mat','PSF_tot','uz', 'defocus', 'NA')
toc
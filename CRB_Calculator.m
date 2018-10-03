% calculating CRB for SAF-detection and finding optimal PSF for SAF-detection by minimizing the CRB:
% -------------------------------------------------------------
% program calculates BFP-fields for various distances (z) of the dipole from the coverslip
% from these data the CRB for x-y-z expectation values are calculated

clear all;

% close all;

global Z ux uk uz Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz Nx n_photon bg mask

%% calculating BFP fields for SAF

%---user parameters----

noise='n';  %set to 'y' or 'n'

N=128;
lambda_0=680e-9;

load('coeff_Example.mat')
[~,~,R,pupil]=create_coord(N,1,'FFT');
modes=[2:37,56];
%         Z_aberr=zeros(1,36);
Zernike_stack=ZernikeCalc(modes,ones(length(modes),1),pupil,'Noll');

Zernike_stack=ZernikeCalc(modes,ones(length(Z_aberr),1),pupil,'Noll');
hugo=zeros(1,1,size(Zernike_stack,3));
hugo(1,1,1:length(Z_aberr))=zeros(1,length(Z_aberr));
%         hugo(1,1,1:length(Z_aberr))=Z_aberr2; 
aberr=sum(Zernike_stack.*repmat(hugo,[N,N,1]),3);



% polaristion dependendent transmission
load('ZernikeCoeff_full.mat'); 
Im_Fit = sum(ZernikeCalc([1 4 5 6 11 12 13 22 37],ZernCoeff,length(pupil),'NOLL'),3);
Im_Fit = Im_Fit/max(Im_Fit(:));
obj_transm=sqrt(Im_Fit);
clear ZernCoeff; 


%n_photon=photonN(ii); %number of camera counts in the brightest dipole-image
n_photon=1000;
bg=200; %background-counts

%     defocus_vec=(defocusPos-1)*4e-8-1.e-6  %defocus of objective lens 
NA=1.65; RI=[1.33 1.33 1.78]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]
% NA=1.49; RI=[1.3 1.3 1.52]; %refractive indices; RI=[RI_specimen, RI_intermed., RI_immoil]

d2=0e-9; %thickness of intermediate layer (layer 2)
f=1.8e-3; %focal length of objective
mu=1e-16; %magnitude of dipole (arbitrary value)

ux=117e-9; %resolution in focal space
Nx=19; %desired simulated field size in pixel
%-----------------------

[SA_out,Defocus,~] =fun_SA_RImismatch(N,RI(3),RI(3),NA,lambda_0,1); %Defocus function refers to refractive index n2

uk=4*pi/lambda_0*NA/N; %unit in pupil space (k-space)
[~,~,R,pupil]=create_coord(N,1,'FFT');

%     obj_transm=polyval([0.0832 0 -0.5199 0 1],asin(2/N*R*NA/RI(3))); %objective intensity transmission function model (from measurement with blue textmarker)
pupil_UAF=circshift(R<=((N/2)*(RI(1)/NA))*1,[0 0]); %pupil containing UAF light

uz=10e-9;
z_vec=(50e-9:uz:50e-9); %simulated dipole distances above layer 2
defocus_vec=(-1.2e-6:uz:1.5e-6);

%z-dipole

%calculating BFP-fields for all dipole orientations
for m=1:length(z_vec);
    dipole=[0,0]; %[theta, phi], e.g. [0,0] for z-dipole, [pi/2,0] for x-dipole
    [Ex_Pz(:,:),Ey_Pz(:,:)]=fun_dipole_imaging(N,lambda_0,NA,RI,[0,0],d2,z_vec(m),f,mu); %z-dipole
    [Ex_Px(:,:),Ey_Px(:,:)]=fun_dipole_imaging(N,lambda_0,NA,RI,[pi/2, 0],d2,z_vec(m),f,mu); %x-dipole
    [Ex_Py(:,:),Ey_Py(:,:)]=fun_dipole_imaging(N,lambda_0,NA,RI,[pi/2, pi/2],d2,z_vec(m),f,mu); %y-dipole
end

%% calculating BFP images

clear I_BFP ratio I_BFP

for m=1:length(defocus_vec); %BFP image number

    %user-defined additional pupil mask (incorporates objective transmission profile):
        phase=0*pupil_UAF;
        mask=pupil.*obj_transm.*exp(1i*aberr+1i*defocus_vec(m)*Defocus);

        I_BFP(m,:,:)=(abs(Ex_Px(:,:).*mask).^2+abs(Ex_Py(:,:).*mask).^2+abs(Ey_Px(:,:).*mask).^2+abs(Ey_Py(:,:).*mask).^2+abs(Ex_Pz(:,:).*mask).^2+abs(Ey_Pz(:,:).*mask).^2);

        ratio(m)=sum(sum((1-pupil_UAF).*squeeze(I_BFP(m,:,:))))/sum(sum(pupil_UAF.*squeeze(I_BFP(m,:,:))));
end

%% -----calculating PSF as seen on the camera for different emitter z-positions-----

fits='n'; %perform fits to the PSFs? This will slow down the loop; choose 'y' or 'n'

clear PSF_tot PSF_SAF PSF_UAF Gfit_UAF Gfit_SAF

boundary=zeros(Nx,Nx); %border-mask to estimate background from PSF-images
tmp=1/(4*(Nx-1)); %normalization such that sum(boundary(:)) equals 1
boundary(1,:)=tmp; boundary(end,:)=tmp; boundary(:,1)=tmp; boundary(:,end)=tmp; 

for m=1:length(defocus_vec);
    mask=pupil.*obj_transm.*exp(1i*phase+1i*aberr+1i*defocus_vec(m)*Defocus);
    %I_CCD_x=abs(fftshift(fft2(ifftshift(embed(pupil.*(E_x(:,:,m)),N_pad,0))))).^2;

    %-----calculating total (SAF+UAF) images-----
    I_xx=abs(czt2(Ex_Px(:,:).*mask,uk,ux,Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:).*mask,uk,ux,Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:).*mask,uk,ux,Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:).*mask,uk,ux,Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:).*mask,uk,ux,Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:).*mask,uk,ux,Nx)).^2;
    PSF_tot(:,:,m)=I_xx+I_yx+I_xy+I_yy+I_xz+I_yz;


    if m==1; C_norm=sum(sum(PSF_tot(:,:,1))); end; %normalization to total intensity in first image (m=1)
    %C_norm=sum(sum(PSF_tot(:,:,m)));  %normalization if info is contained in shape-changes

    if strcmp(noise,'y'); %if noise is selected
        tmp=PSF_tot(:,:,m)/C_norm*n_photon+bg; 
        PSF_tot(:,:,m)=poissrnd(tmp,size(tmp,1),size(tmp,2)); %normalization of PSF   
    else
        PSF_tot(:,:,m)=PSF_tot(:,:,m)/C_norm; %normalization of PSF   
    end
    est_offset=sum(sum(boundary.*PSF_tot(:,:,m))); %returns the mean value in the boundary-region
    energy_tot(m)=sum(sum(PSF_tot(:,:,m)-est_offset));


    %-----calculating UAF-images-----
    I_xx=abs(czt2(Ex_Px(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;
    I_yx=abs(czt2(Ey_Px(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;
    I_xy=abs(czt2(Ex_Py(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;
    I_yy=abs(czt2(Ey_Py(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;    
    I_xz=abs(czt2(Ex_Pz(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;
    I_yz=abs(czt2(Ey_Pz(:,:).*mask.*pupil_UAF,uk,ux,Nx)).^2;
    PSF_UAF(:,:,m)=I_xx+I_yx+I_xy+I_yy+I_xz+I_yz;    

    if strcmp(noise,'y'); %if noise is selected
        tmp=PSF_UAF(:,:,m)/C_norm*n_photon+bg; 
        PSF_UAF(:,:,m)=poissrnd(tmp,size(tmp,1),size(tmp,2)); %normalization of PSF   
    else
        PSF_UAF(:,:,m)=PSF_UAF(:,:,m)/C_norm; %normalization of PSF   
    end
    est_offset=sum(sum(boundary.*PSF_tot(:,:,m))); %returns the mean value in the boundary-region
    energy_UAF(m)=sum(sum(PSF_UAF(:,:,m)-est_offset));



end

PSF_tot=imgaussfilt(PSF_tot, 0.125);

z_plot=round([1, 2, 2*length(z_vec)/3, length(z_vec)]);

[CRBx,CRBy,CRBz]=fun_CRB(PSF_tot./repmat(sum(sum(PSF_tot,1),2),[Nx Nx 1]),ux,uz,n_photon,bg);

CRBz=sqrt(CRBz);

CRB_Z=CRBz;

metric1=mean(sqrt((CRBx.*CRBy.*CRBz)).^(1/3));  %"localization volume"
metric2=mean(sqrt(CRBz)); 


%%

figure(4)
plot(defocus_vec(1:end-1)*10^6,CRB_Z*sqrt(2),'b-'); xlabel('Defocus in \mum'); ylabel('precision nm');
grid on;
title('z-CRB for system aberrations');
legend('CRB','MLE','CNN')

%%
figure(9)
plot(defocus_vec(1:end-1)*10^6,sqrt(CRBx))
xlabel('Defocus in \mum'); 
ylabel('precision nm');
title('CRB for system aberrations xy-plane');
legend('CRB','CNN-x','CNN-y','MLE-x','MLE-y')

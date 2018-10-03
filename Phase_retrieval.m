%Phase retrieval to measure pupil field
%VERSION USED FOR OPEX PAPER
% an image z-stack of a fluorescent bead (sub-diffraction-limited size) is
% used to retrieve the pupil function
% fully vectorial simulation

%copy the stack (tif) into the working directory of
%matlab and run the code

clear all; 
% close all;
clc; 
            
%% user parameters 
    
    pol='full';    %is a polarizer used? %type 'x', 'y', or 'full' for defining the passing polarization direction; "unpol" means that no polarizer is used
    NA=1.67;
    RI=[1.45 1.78 1.78]; %RIs of liquid, coverglass and imm.oil
    dia_pupil=64;        %diameter of pupil to be retrieved in pixel
    Nx=40;               %focal plane size of simulated field; don't choose too small, otherwise artefacts will results if the lights "leaks" out of the boundaries at high defocus values
    pixsize_CCD=13e-6;   %camera pixel size
    M=100*200/180;              %magnification of objective lens
    lambda_0=680e-9;     %vacuum wavelength of fluorescence light
    dz=0.2e-6;            %axial interspacing between two adjacent frames of the image stack
    dia_bead=170e-9;     %diameter of test-bead
    f_obj=180e-3/M;      %tube lens focal length / magnification
    ux_desired=pixsize_CCD/M; %length of a pixel in the focal plane
    %ux_desired=100e-9;
    
    k0=2*pi/lambda_0;
    [Kx,Ky,Kr,pupil,defocus,N_pad,uk,ux,v]=fun_usualsuspects(dia_pupil,lambda_0,ux_desired,NA,RI(3));
    
ALPHA=asin(Kr/k0/RI(3)).*pupil; %pupil angular coordinate

% polaristion dependendent transmission
if strcmp(pol,'y');
    load('ZernikeCoeff_y.mat');
elseif strcmp(pol,'x');
    load('ZernikeCoeff_x.mat');
elseif strcmp(pol,'full');
    load('ZernikeCoeff_full.mat');
end
  
Im_Fit = sum(ZernikeCalc([1 4 5 6 11 12 13 22 37],ZernCoeff,length(pupil),'NOLL'),3);
Im_Fit = Im_Fit/max(Im_Fit(:));
T_obj=sqrt(Im_Fit);
clear ZernCoeff; 


% Apo=pupil./sqrt(cos(ALPHA)).*T_obj; %combined effect of objective transmission and aplanatic factor of objective lens E~srqt(cos(theta))
Apo=pupil.*T_obj; %combined effect of objective transmission and aplanatic factor of objective lens E~srqt(cos(theta))

%scalar: initial pupil field
%E_pup=pupil./sqrt(cos(ALPHA)).*T_obj;


%% import raw data / centering of images / removing background

% clear I I2 I3 I5 I6 BFP Mess transmis; 


% Import Fits-Stack
import matlab.io.*
[filename,PathName] = uigetfile('*.tif','Select the FITS file');
%fname = 'single_bead.tif';
%fname = 'UAF_stack_PSspeckBeads_dz1µm_crop.tif';
%fname = 'Stack_B_final.tif';
%fname = 'bead_2018-03-22_stack4.tif';
%fname = 'avg-stack2_SAF.tif';
%fname = 'avg-stack2_UAF_flipped.tif';
fname = [PathName filename];

info = imfinfo(fname);
sx=info(1).Height; %size of an image
sy=info(1).Width;

delta=0; %number of images to be deleted from the beginning AND the end of the stack

kk=0;

num_images = numel(info)-2*delta-1;
for m = 1:num_images
    
    kk=kk+1
    
    I(:,:,m) = double(imread(fname, m+delta));
    Itmp(:,:,m)=I(:,:,m)/trapz(trapz(I(:,:,m)));  %normalizing each frame in the stack individually; this is to compensate for different exposure times
    I_max(m)=max(max(Itmp(:,:,m)));
end

[~, m_focus]=max(I_max); %finding the in-focus image
 %m_focus=9; 

%centering the PSF within the frame
    I_infocus=I(:,:,m_focus);
    [~, max_idx]=max(I_infocus(:));
    [maxrow,maxcol] = ind2sub([sx sy],max_idx); %find max position
    I2=circshift(I,round(([sx sy]+1)/2)-[maxrow, maxcol]); clear I;
    
%cropping to user-defined size
    %I3=embed(I2,[Nx Nx,num_images],0);
    I3=I2;
    [X,Y,R,mask2]=create_coord(sx,1,'exact'); 
   
    
%removing background (method of Liu et al.2013: for each cropped z-plane
%image, the background is the smalles of the mean edge values
    clear BG
    BG(1,1,:)=min([squeeze(mean(I3(1,:,:),2)),squeeze(mean(I3(end,:,:),2)),squeeze(mean(I3(:,1,:),1)),squeeze(mean(I3(:,end,:),1))],[],2);
    %BG=mean2(I2(1,:,idx)); %background 
    %BG=10;
    I4=I3-repmat(BG,[sx sx 1]); %clear I2;
    I4(I4<0)=0; %ensuring positivity
    %I4=I4.*repmat(mask2,[1,1,num_images]); %applying binary circular mask
    
%thresholding to remove noise
    thresh=0.0*max(max(I4(:,:,1)));
    I4(I4<thresh)=0;
  
%     num_images=11;
    
  for m=1:num_images;
    I4(:,:,m)=I4(:,:,m)/trapz(trapz(I4(:,:,m)));  %normalizing each frame in the stack individually; this is to compensate for different exposure times
  end
  %I4(:,:,3:5)=0.8*I4(:,:,3:5); plot(squeeze(sum(sum(I4,1),2)),'-o')

  Image_stack=I4/trapz(I4(:));

    
%% simulate raw data (optional)
% 
% num_images=9;
% dz=200e-9;
% m_focus=3;
% sx=50;
% sy=50;
% aberr_vec=[zeros(1,4) 0.5 -0.3 0.1 -0.1 0 0 0.4 0 0 0 0 0 0 0 0 0 0 0.2];
% aberr=sum(ZernikeCalc(1:length(aberr_vec),aberr_vec',pupil,'Noll'),3);
% 
% clear I4 I5 I6;
% 
% %simulate PSF
% for m = 1:num_images
%     mx=1; %magnitude of x-dipole
%     my=1;
%     mz=1;
%     tmp=fun_dipole_imaging_full(T_obj.*exp(1i*aberr),sx,lambda_0,ux,NA,RI,[mx my mz],(m-m_focus)*dz,'full');
%     I4(:,:,m)=tmp;
%   
% %         %add noise
% %         photon_no=20000;
% %         tmp=tmp/sum(tmp(:))*photon_no;
% %         I4(:,:,m)=poissrnd(tmp);
%     
%     I5(:,:,m)=embed(tmp,dia_pupil,0);
%     I6=I5;
% end

%if raw data are simulated, bead-deconvolution is not required --> skip
%corresponding block    
    
%I5=I4/sum(I4(:)); %image stack to be used for phaes retrieval; should be normalized


%% pre-calculating some data required for simulating bead-images
   
    d2=0; %thickness of intermediate layer (see Axelrod); assumed to be zero in our case
    z=50e-9; %distance of bead from glass; 
    f_obj=1.8e-3; %objective focal length; only important for quantitative correct field magnitudes
    dipoles=[1 1 1]*1e-12; %vector defining dipole strengths [Px Py Pz]

    %fun_dipole_imaging(N,lambda_0,NA,RI,dipole,varargin)
    %calculation of BFP-fields for each dipole
    [Ex_Pz,Ey_Pz,~,pupil]=fun_dipole_imaging(dia_pupil,lambda_0,NA,RI,[0 0],d2,z,f_obj,dipoles(3)); %z-dipole
    [Ex_Px,Ey_Px,~,~]=fun_dipole_imaging(dia_pupil,lambda_0,NA,RI,[pi/2 0],d2,z,f_obj,dipoles(1)); %x-dipole
    [Ex_Py,Ey_Py,~,~]=fun_dipole_imaging(dia_pupil,lambda_0,NA,RI,[pi/2 pi/2],d2,z,f_obj,dipoles(2)); %y-dipole

    E_BFP(:,:,1)=Ex_Px.*T_obj;
    E_BFP(:,:,2)=Ex_Py.*T_obj;
    E_BFP(:,:,3)=Ex_Pz.*T_obj;
    E_BFP(:,:,4)=Ey_Px.*T_obj;
    E_BFP(:,:,5)=Ey_Py.*T_obj;
    E_BFP(:,:,6)=Ey_Pz.*T_obj;

    if strcmp(pol,'x')
        E_BFP(:,:,4:6)=[];
    elseif strcmp(pol,'y')
        E_BFP(:,:,1:3)=[];
    end
    
    figure(9)
    imagesc(sum(abs(E_BFP).^2,3))
    
%defining parameters for chirped z-trafo
    Nk=dia_pupil; %size of input field (in k-space)
    %N_pad=2^nextpow2(Nx+Nk-1); %padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
    N_pad=Nx+Nk-0;
    x=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5))';
    y=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5));
    ux_fft=2*pi/(Nk*uk); %fft resolution without padding
    r=ux/ux_fft; %required r-factor to meet the desired resolution at the given grid size N
    alpha=r*1/Nk;
    kernel=exp(-1i*alpha*pi*x.^2)*exp(-1i*alpha*pi*y.^2); %create quadratic convolution phase kernel; faster method
    F_kernel=fft2(ifftshift(kernel));
           
%% considering effects of the bead-size
%simulated PSF is low-pass filtered (see Hanser et. al., J. of Microsc. 2004)

Nz=size(I4,3);
ukz=2*pi/dz/Nz; %k-space z-coordinate unit
ukx=2*pi/ux/Nx;
kz=((1:Nz)-m_focus)*ukz; %kz-axis
kx=(-Nx/2:Nx/2-1)*ukx;
ky=kx;
[Kx3D,Ky3D,Kz3D]=ndgrid(kx,ky,kz);
Kr3D=sqrt(Kx3D.^2+Ky3D.^2+Kz3D.^2); %3D k-space radial coordinate

z=1/2*Kr3D*dia_bead;
b=3*(sin(z)./z.^3-cos(z)./z.^2); %Fourier Transform of the  bead shape
b(isnan(b))=1; %remove nan in the centre

%the 3F-transform of the simulated 3D-PSFs is multiplied with this function b
%this happens in the sub-routine "fun_phase_retrieval_vec_FAST"

%% compare calculated ideal and measured PSfs

clear I_simu tmp;
figure(4); 
colormap gray;

clear I_scalar I_vec_tot I_vec
I_vec=zeros(Nx,Nx,size(E_BFP,3));

for m=1:num_images;
           
    %-----scalar simulation-----
    I_scalar(:,:,m)=embed(abs(fftshift(fft2(ifftshift(embed(pupil./sqrt(cos(ALPHA)).*T_obj.*exp(1i*defocus*(m-m_focus)*dz),N_pad,0))))).^2,sx,0);
    
    %-----vectorial simulation-----
    for q=1:size(E_BFP,3); %propagating all BFP-fields to camera
        I_vec(:,:,q)=abs(czt2(E_BFP(:,:,q).*exp(1i*defocus*(m-m_focus)*dz),uk,ux,Nx)).^2;
    end    
    I_vec_tot(:,:,m)=sum(I_vec,3);    

    %visualization
    subplot(3,1,1);
    imagesc(squeeze(I_scalar(:,:,m))); title('scalar'); axis equal; axis tight;
    subplot(3,1,2);
    imagesc(squeeze(I_vec_tot(:,:,m))); title(['vect. pol=' pol]); axis equal; axis tight;
    subplot(3,1,3);
    imagesc(squeeze(embed(Image_stack(:,:,m),[Nx Nx],0))); title('measured'); axis equal; axis tight;
    pause(.5); 
    disp(m);      
end

%% show sagittal projections (optional)

figure(5);
subplot(2,1,1);
imagesc((1:Nx)*ux*1e6,((1:num_images)-m_focus)*dz*1e6,squeeze(sum(I4,2))'); axis equal; axis tight;
xlabel('x / µm');
ylabel('z / µm');
title('measurements');
subplot(2,1,2); 
imagesc((1:Nx)*ux*1e6,((1:num_images)-m_focus)*dz*1e6,squeeze(sum(I_vec_tot,2))');  axis equal; axis tight;
xlabel('x / µm');
ylabel('z / µm');
title('vec. simulation');


%% phase retrieval

z_vec=((1:num_images)-m_focus)*dz; %vector of z-positions 0=in-focus

%pre-calculating Zernike mode stack; piston is omitted
% maxmode=37;
% Zernike_stack=ZernikeCalc(2:maxmode,ones(maxmode,1),pupil,'Noll');

modes=[2:37,56];
%modes=[2:4,11,22,37,56];
Zernike_stack=ZernikeCalc(modes,ones(length(modes),1),pupil,'Noll');
          
%--minimum search--
Z_ini=zeros(1,length(modes)); %initial estimates
% Z_ini=;

tic
%tmp=fminsearch(@(Z_ini) fun_phase_retrieval_vec(Z_ini,Image_stack, defocus, z_vec, ux, uk, E_BFP, dia_pupil, Zernike_stack, T_obj,'n'),Z_ini);
options=optimoptions('fminunc');
options.Algorithm='quasi-newton';
options.TypicalX=0.5*ones(length(Z_ini),1); %typical values of Zernike coefs
options.FunctionTolerance=0.0001;
options.StepTolerance=0.0001;
options.MaxFunEvals=1000;

% FiniteDifferenceStepSize=0.0125;
tmp=fminunc(@(Z_ini) fun_phase_retrieval_vec_FAST(Z_ini, Image_stack, kernel, F_kernel, defocus, z_vec, E_BFP, dia_pupil, Zernike_stack,b, 'n'),Z_ini,options);

% % % bound=[2 2 2 0.3*ones(1,6) 2 0.3*ones(1,10) 2 0.3*ones(1,14) 2 2];
% % % tmp=fminsearchbnd(@(Z_ini) fun_phase_retrieval_vec_FAST(Z_ini, Image_stack, kernel, F_kernel, defocus, z_vec, E_BFP, dia_pupil, Zernike_stack,b, 'n'),Z_ini,-bound,+bound);

disp('done');
toc

Z_aberr=zeros(1,length(modes));
Z_aberr=tmp;

phase=sum(ZernikeCalc(modes,[0 0 0 Z_aberr(4:end)]',pupil,'Noll'),3);

% removing spherical defocus
%[a_SD, SD, red_mask]=spherical_defocus_analysis(phase,NA,RI(3),RI(1));
defocus_coefs=ZernikeCalc(modes,defocus,pupil,'Noll'); %zernike coefs of high-NA defocus
a_def=dot(defocus_coefs/norm(defocus_coefs),Z_aberr); 
Z_aberr2=(Z_aberr-a_def*defocus_coefs'/norm(defocus_coefs)); %defocus-free Zernike-coefs


%%
save('coeff_Example.mat','Z_aberr','Z_aberr2')

%% visualization

phase=sum(ZernikeCalc(modes,[0 0 0 Z_aberr2(4:end)]',pupil,'Noll'),3);

%calculating Zernike correction coefs for specific lab setup, e.g. rotation
rotangle=0;
% rotangle=91.5;
figure(33);
subplot(2,1,1);
Z_setup=ZernikeCalc(modes,imrotate(phase,rotangle,'crop'),pupil,'Noll',3)';
bar(modes,[0 0 0 Z_setup(4:end)]); xlabel('Zernike mode no.'); ylabel('rad RMS'); grid on;
title('Zernikes for setup; high-NA defocus removed');

subplot(2,1,2);
imagesc(imrotate(phase,rotangle,'crop'));
title('retr. phase (wo. tilt)'); colormap parula; colorbar; axis equal; axis tight; %re-calculate Zernikes for rotated phase pattern

M=[modes(4:end);Z_setup(4:end)]';
%dlmwrite('corr_coefsB_stackB.txt',M,'delimiter','\t','precision',2,'newline','pc');


%% check of retrieved PSF
 
%load('Z_aberr_vec.mat')  
[error, I_retr]=fun_phase_retrieval_vec_FAST(Z_aberr(1:end), Image_stack/sum(Image_stack(:)), kernel, F_kernel, defocus, z_vec, E_BFP, dia_pupil, Zernike_stack,b,'y');

%% show sagittal projections

figure(5);
subplot(3,1,1);
imagesc((1:Nx)*ux*1e6,((1:num_images)-m_focus)*dz*1e6,squeeze(sum(I_vec_tot,2))'); axis equal; axis tight;
xlabel('x / µm');
ylabel('z / µm');
title('vec.simu');

subplot(3,1,2);
imagesc((1:Nx)*ux*1e6,((1:num_images)-m_focus)*dz*1e6,squeeze(sum(I4,2))'); axis equal; axis tight;
xlabel('x / µm');
ylabel('z / µm');
title('measurements');
subplot(3,1,3);
imagesc((1:Nx)*ux*1e6,((1:num_images)-m_focus)*dz*1e6,squeeze(sum(I_retr,2))');  axis equal; axis tight;
xlabel('x / µm');
ylabel('z / µm');
title('retrieval result');
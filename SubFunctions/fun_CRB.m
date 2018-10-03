function [CR_x, CR_y, CR_z, CR_sig, CR_bg]=fun_CRB(PSF,ux,uz,n_photon,bg)
%calculating the CR-lower-bound for x/y/z localization
%output is in nm^2
%input is a 3D PSF 
%ux...lateral resolution 
%uz...axial resolution
%bg...background level (mean) per pixel
%n_photon...total number of photons in 3DPSF

N=size(PSF); %assumes x,y,z-sizes

%adding noise
    
    PSF_s=PSF/sum(sum(PSF(:,:,1)))*n_photon; %normalization: every z-slice contains n_photon number of photons
    
    shot=PSF_s; %variance of shot noise
    readout=0;     %variance of dark/readout noise
    noise=shot+readout+bg; %variance of noise (here shot noise: variance = mean value)

Dz=(diff(PSF_s,1,3)/uz); %intensity differential along z
Dx=(diff(PSF_s,1,1)/ux);
Dy=(diff(PSF_s,1,2)/ux);
Dbg=1;
Dsig=PSF_s/sum(sum(PSF_s(:,:,1)));


if isempty(Dz); dim=2; else dim=3; end; %if only a 2D PSF was provided 

LI=1:(size(PSF,1)-1);
AI=1:(size(PSF,3)-1); 


if dim==2; %if only a 2D PSF was provided 
    Dx=Dx(LI,LI); %bring tensors into same size
    Dy=Dy(LI,LI);
    Dsig=Dsig(LI,LI);
    noise=noise(LI,LI);

    %calculating Fisher information matrix entries
    FI_xx=squeeze(sum(sum(Dx.^2./noise,1),2))*1e-18; 
    FI_yy=squeeze(sum(sum(Dy.^2./noise,1),2))*1e-18; 
    FI_xy=squeeze(sum(sum(Dx.*Dy./noise,1),2))*1e-18;
    FI_yx=FI_xy;
    
    FI_ss=squeeze(sum(sum(Dsig.^2./noise,1),2)); %in 1/n_photon^2
    FI_bb=squeeze(sum(sum(ones(size(Dsig))*Dbg^2./noise,1),2)); %in 1/n_photon^2
    FI_xs=squeeze(sum(sum(Dx.*Dsig./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_ys=squeeze(sum(sum(Dy.*Dsig./noise,1),2))*1e-9; %in 1/nm*1/n_photon

    FI_xb=squeeze(sum(sum(Dx.*Dbg./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_yb=squeeze(sum(sum(Dy.*Dbg./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_bs=squeeze(sum(sum(Dsig.*Dbg./noise,1),2)); %in 1/n_photon^2


    FI=zeros(4,4); %Fisher matrix
    FI(1,1)=FI_xx;
    FI(2,2)=FI_yy;
    FI(3,3)=FI_ss; 
    FI(4,4)=FI_bb; 
    
    FI(1,2)=FI_xy; FI(2,1)=FI_xy;
    FI(1,3)=FI_xs; FI(3,1)=FI_xs;
    FI(1,4)=FI_xb; FI(4,1)=FI_xb; 
    FI(2,3)=FI_ys; FI(3,2)=FI_sy; 
    FI(2,4)=FI_yb; FI(4,2)=FI_yb; 
    FI(3,4)=FI_bs; FI(4,3)=FI_bs; 
    
    
    
%     % visualizing fisher matrix entries
    %figure(22); colormap gray;
    %imagesc(FI); axis equal; axis tight; colorbar; title('Fisher Matrix');
    
    % calculating CRAMER RAO bounds
    %inverting Fisher matrix and taking diagonal values
        FI_inv=inv(FI);
        CR_x=FI_inv(1,1);
        CR_y=FI_inv(2,2);
        CR_z=0;
        CR_sig=FI_inv(3,3);
        CR_bg=FI_inv(4,4);
  

elseif dim==3; %if a 3D PSF was provided 
    Dx=Dx(LI,LI,AI); %bring tensors into same size
    Dy=Dy(LI,LI,AI);
    Dz=Dz(LI,LI,AI);
    Dsig=Dsig(LI,LI,AI);
    noise=noise(LI,LI,AI);
    
    %calculating Fisher information matrix entries
    FI_zz=squeeze(sum(sum(Dz.^2./noise,1),2))*1e-18; %in 1/nm^2
    FI_xx=squeeze(sum(sum(Dx.^2./noise,1),2))*1e-18; 
    FI_yy=squeeze(sum(sum(Dy.^2./noise,1),2))*1e-18; 
    FI_ss=squeeze(sum(sum(Dsig.^2./noise,1),2)); %in 1/n_photon^2
    FI_bb=squeeze(sum(sum(ones(size(Dsig))*Dbg^2./noise,1),2)); %in 1/n_photon^2
    
    FI_xy=squeeze(sum(sum(Dx.*Dy./noise,1),2))*1e-18;
    FI_xz=squeeze(sum(sum(Dx.*Dz./noise,1),2))*1e-18;
    FI_yz=squeeze(sum(sum(Dy.*Dz./noise,1),2))*1e-18;
    
    FI_xs=squeeze(sum(sum(Dx.*Dsig./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_ys=squeeze(sum(sum(Dy.*Dsig./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_zs=squeeze(sum(sum(Dz.*Dsig./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_xb=squeeze(sum(sum(Dx.*Dbg./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_yb=squeeze(sum(sum(Dy.*Dbg./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_zb=squeeze(sum(sum(Dz.*Dbg./noise,1),2))*1e-9; %in 1/nm*1/n_photon
    FI_bs=squeeze(sum(sum(Dsig.*Dbg./noise,1),2)); %in 1/n_photon^2
             
    FI=zeros(5,5,size(PSF,3)-1); %Fisher matrix
    
    %diagonal elements
    FI(1,1,:)=FI_xx;
    FI(2,2,:)=FI_yy;
    FI(3,3,:)=FI_zz;
    FI(4,4,:)=FI_ss;
    FI(5,5,:)=FI_bb;
    
    %off-diagonal elements
    FI(1,2,:)=FI_xy; FI(2,1,:)=FI_xy;
    FI(1,3,:)=FI_xz; FI(3,1,:)=FI_xz;
    FI(1,4,:)=FI_xs; FI(4,1,:)=FI_xs; 
    FI(1,5,:)=FI_xb; FI(5,1,:)=FI_xb;
    
    FI(2,3,:)=FI_yz; FI(3,2,:)=FI_yz;
    FI(2,4,:)=FI_ys; FI(4,2,:)=FI_ys;
    FI(2,5,:)=FI_yb; FI(5,2,:)=FI_yb;
    
    FI(3,4,:)=FI_zs; FI(4,3,:)=FI_zs;
    FI(3,5,:)=FI_zb; FI(3,5,:)=FI_zb;
    
    FI(4,5,:)=FI_bs; FI(5,4,:)=FI_bs;
    
       
    
    
%     % visualizing fisher matrix entries
    figure(22); colormap gray;
    imagesc(FI(:,:,1),[-.1 .1]); axis equal; axis tight; title('Fisher Matrix');
    
    % calculating CRAMER RAO bounds
    for m=1:(size(PSF,3)-1); %inverting Fisher matrix and taking diagonal values
        FI_inv=inv(FI(:,:,m));
        CR_x(m)=FI_inv(1,1);
        CR_y(m)=FI_inv(2,2);
        CR_z(m)=FI_inv(3,3);
        CR_sig(m)=FI_inv(4,4);
        CR_bg(m)=FI_inv(5,5);
    end
end
    



% figure(4);
% subplot(2,2,1);
% az=((1:length(FI_yy))*dz-length(FI_zz)/2*dz)*1e6;
% plot(az,CR_z,colour); ylim([0 60]);
% title('CRB (z-direction)');
% ylabel('nm^2');
% xlabel('defocus / µm');
% hold on;
% 
% subplot(2,2,2);
% ax=((1:length(FI_xx))*dz-length(FI_xx)/2*dz)*1e6;
% plot(ax,CR_x,colour); ylim([0,30]);
% ylabel('nm^2');
% xlabel('defocus / µm');
% title('x-direction');
% hold on;
% 
% subplot(2,2,3);
% plot(ax,CR_y,colour); ylim([0,30]);
% ylabel('nm^2');
% xlabel('defocus / µm');
% title('y-direction');
% hold on;
% %%



function [x,y,z,N,BG,ZEST] = MLE_Estimator_TOT(PSF_img,ZEst,Int)

global pxSize DAT I PSF_tot

%load('PSF_SAF_5nm.mat') %data exported from "optimal_SAF_PSF.m"; variables: PSF_tot, PSF_UAF, ux, z_vec, RI, NA

uz=5; %z-increment in PSF-stack

%  normalization of energy in each z-slice
energies=(sum(sum(PSF_tot,1),2)); 
PSF_norm=PSF_tot./repmat(energies,[size(PSF_tot,1),size(PSF_tot,2),1]);
% dz_PSF=zeros(size(PSF_tot));
% dz_PSF(:,:,1:end-1)=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
dz_PSF=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
PSF_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF

fw=2; %frame-width; assumed molecule image must be smaller than the PSF-stack images! Ohterwise, NaNs appear when interpolating at different x-y positions
      % frame-width must be larger than any possible x-y-shift; 

      
%create coordinate system 
[nx0 ny0 nz]=size(PSF_tot);
x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
[X,Y]=ndgrid(x,y);

x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
[nx ny]=size(X); %size of molecule image


%% taking out an arbitraty x-y slice of the PSF-stack: this is our "measurement"

photon_no=3000;
BG=30000; %number of total background photons

I=PSF_img;

fun = @(x)estimator3D(x);    
        
        %initial BG estimate
        BGmask1=zeros(size(I));
        BGmask2=zeros(size(I));
        BGmask1(1,:)=1; BGmask1(end,:)=1; BGmask1(:,1)=1; BGmask1(:,end)=1;
        BGmask2(2,2:end-1)=1; BGmask2(end-1,2:end-1)=1; BGmask2(2:end-1,2)=1; BGmask2(2:end-1,end-1)=1;
        BG2=sum(BGmask2(:).*I(:))/sum(BGmask2(:)); %mean signal in inner frame
        BG1=sum(BGmask1(:).*I(:))/sum(BGmask1(:)); %mean signal in outer frame
        BG0=abs(BG1-(BG2-BG1)*2);
        BG0=BG1;
        %BG0=BG/nx/ny; %take true BG as initial estimate: this should be the best
        D=(centerOfMass(abs(I))-floor(10/2+1));
        
        
% % % %         initial estimates of x-y-position,BG,N0 using Gaussfit
% % %         param_ini=[BG0 max(I(:))-BG0 5+D(1) 5+D(2) 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
% % %         lb=[0 0.5*(max(I(:))-BG0) 0 0 0.15]; %lower bounds for fit-parameters
% % %         ub=[2*BG0 1.5*(max(I(:))-BG0) 20 20 3]; %upper bounds
        
        
        param_ini=[BG0 max(I(:))-BG0 0 0 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
        lb=[0 0.5*(max(I(:))-BG0) -nx/3 -nx/3 0.5]; %lower bounds for fit-parameters
        ub=[2*BG0 1.5*(max(I(:))-BG0) nx/3 nx/3 3]; %upper bounds
        
        [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_gauss_and_offset_test,param_ini,x_data,I,lb,ub);
        Gaussfit=fun_gauss_and_offset(Param,x_data);
        N0=Param(2)*2*pi; %energy contained (see e.g. formula in Thunderstorm script)
        x0=Param(3);
        y0=-Param(4);
        BG0=Param(1);
        N0=sum(I(:))-BG0*nx*ny; %initial guess for number of photons
        
        %----MLE estimation----------        
        var0=[x0 y0 round(nz/2) N0 BG0]; %initial estimates
        LB=[-fw/2 -fw/2 1 0.75*N0 0.1*BG0]; %lower bounds
        UB=[+fw/2 +fw/2 nz 1.5*N0 BG0]; %upper bounds     

        Z_EST = 1+abs((Param(5)-.45)/.45*45);
        
        if Z_EST > length(PSF_tot)
            Z_EST=length(PSF_tot)-1;
        elseif Z_EST <1
            Z_EST=1;
        end
        
        if ZEst > length(PSF_tot)
            ZEst=length(PSF_tot)-1;
        elseif ZEst <1
            ZEst=1;
        end
        
        ZEst
        
        options = optimoptions(@fminunc,'DiffMinChange',.000001,'OptimalityTolerance',1e-10);
        optionsFmin = optimset('MaxIter',10000,'MaxFunEvals',10000);
%         tmp=fminunc(fun,[100 N0 BG0 D(1) D(2)],options)
        tmp=fminsearch(fun,[50 N0 BG0 D(1) D(2)],optionsFmin)

        ZEST=tmp(1);
        x=tmp(4);
        y=tmp(5);
        z=tmp(1);
        N=tmp(2);
        BG=tmp(3);
        
% %         x=tmp(1);
% %         y=tmp(2);
% %         z=tmp(3);
% %         N=tmp(4);
% %         BG=tmp(5);


% imagesc(I)
% title(['MLE of ZU-positoin; Z=' num2str(z)]);
% colorbar;
% pause(.25)
end
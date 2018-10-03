%evaluate experimental molecule Data --> 3D localisations

clear all; 

%% user parameters 

%camera parameters
gain=50; 
amp=9.8; %electrons per count 
QE=0.95; %quantum efficiency


%% loading PSF model

%vectashield-PSFs:
%load('PSF_SAF_NA1,7_0-2nm-200nm_RI=1,45_dz=-300nm.mat') %defocused PSF-model; should give better z-estimates
%load('PSF_SAF_0-2nm-200nm_RI=1,45_dz=-300nm_2018-09-26.mat'); 
load('PSF_SAF_5nm.mat'); 

PSF=PSF_tot; 
clear PSF_tot; 

nx0=size(PSF,1); %size of model
ny0=nx0; 

%PSF: normalization of energy in each z-slice
energies=(trapz(trapz(PSF,1),2)); 
PSF_norm=PSF./repmat(energies,[size(PSF,1),size(PSF,2),1]);
dz_PSF=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
PSF_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF

fw=3;


%% reading in data

[Name,Path,~]=uigetfile('*.tif','choose multipage tif data','C:\Users\q004aj\Documents\PROJEKTE\SAF-TIRF_Gerhard Schütz - Immunolog. Synapse\Data');
    image_no=length(imfinfo([Path Name]));
    for m=1:image_no;
      data(:,:,m)=imread([Path Name],m);
    end
disp('done');
    
%% <<<<<<<<<<<<<< METHOD A : click & evaluate
for j=1:1
    % choose image 
    img_no=j;
    figure(1);
    imagesc(data(:,:,img_no)); colormap gray; axis equal; axis tight; 
    title('exp. image');

    % click onto molecules

    figure(1);
    imagesc(data(:,:,img_no)); colormap gray; axis equal; axis tight; 
    title('exp. image');
    [px, py] = getpts(gcf); %click on several molecules and press enter 
    no_mols=length(px); 

    clear tmp;
    for m=1:no_mols;
        tmp(:,:,m)=data(round((py(m)-floor((nx0-1)/2)+fw):(py(m)+ceil((nx0-1)/2))-fw),round((px(m)-floor((nx0-1)/2)+fw):(px(m)+ceil((nx0-1)/2)-fw)),img_no);
    end

    data2=double(tmp/gain*amp/QE); %converting into units of photons; 

    % imagesc(data2(:,:,2)); title('fist exp. image / photons'); axis equal; axis tight; colorbar; 


    %--- Estimating parameters ---

            %create coordinate system 
            [nx0, ny0, nz]=size(PSF_norm);
            x=(fw+1:nx0-fw)-ceil((nx0+1)/2);
            y=(fw+1:ny0-fw)-ceil((nx0+1)/2);
            [X,Y]=ndgrid(x,y);
            clear x_data;
            x_data(:,:,1)=X; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
            x_data(:,:,2)=Y; %coord. data in this form is required for Gaussfits which are used to find initial estimates for the MLE fit
            [nx, ny]=size(X); %size of molecule image
            %defining background mask for initial BG-estimate: 

            BGmask1=zeros(nx,ny);
            BGmask1(1,:)=1; BGmask1(end,:)=1; BGmask1(:,1)=1; BGmask1(:,end)=1;

    clear x_est y_est z_est z_est_nm BG_est N_est
    for m=1:no_mols; 

            I=data2(:,:,m); %pick one experimental image

            fun_LLH=@(v) sum(sum(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))-I.*log(v(5)+v(4)*interpn(PSF_norm,(fw+1:nx+fw)-v(1),(fw+1:ny+fw)'-v(2),v(3))),1),2);


            BG00=sum(BGmask1(:).*I(:))/sum(BGmask1(:)); %initial background est.

            %----Gaussfit for initial estimates of x-y-position,BG,N0-----
            param_ini=[BG00 max(I(:))-BG00 0 0 1]; %[offset amplitude x-shift y-shift width]; initial parameters for Gaussfit
            lb=[0 0.5*(max(I(:))-BG00) -nx/3 -nx/3 0.5]; %lower bounds for fit-parameters
            ub=[2*BG00 1.5*(max(I(:))-BG00) nx/3 nx/3 3]; %upper bounds
            [Param,resnorm,residual,exitflag]=lsqcurvefit(@fun_gauss_and_offset_test,param_ini,x_data,I,lb,ub);
            Gaussfit=fun_gauss_and_offset_test(Param,x_data);
            x0=Param(3);
            y0=Param(4);
            BG0=Param(1); %this serves as initial estimate for the MLE fit
            %N0=Param(2)*2*pi*Param(5)^2; %energy contained (see e.g. formula in Thunderstorm script)
            N0=sum(I(:))-BG0*nx*ny; %initial guess for number of photons - seems to work better 

            %----MLE estimation----------    
            z_ini=80; %initial z-estimate (frame-no.)
            gauss_est(m,:)=[x0 y0 round(nz/2) N0 BG0]; %initial estimates (e.g. from Gaussfit)
            LB=[-nx/3 -nx/3 0 0.5*N0 0.5*BG0]; %lower bounds
            UB=[+nx/3 +nx/3 nz 2*N0 2*BG0]; %upper bounds
            %tmp=fmincon(fun_LLH,x0,[],[],[],[],[1 0.75*N_ini 0.75*BG_ini],[nz 1.5*N_ini 1*BG_ini]);
            tmp=fminsearchbnd(fun_LLH,gauss_est(m,:),LB,UB);
            x_est(m)=tmp(1);
            y_est(m)=tmp(2);
            z_est(m)=tmp(3)
            N_est(m)=tmp(4);
            BG_est(m)=tmp(5);

    %         z_est_nm(m)=interp1(1:length(z_vec),z_vec,tmp(3))*1e9 %estimated z-position in nm
    % 
    %         imagesc(I); axis equal; axis tight; title(['signal=' num2str(N_est(m))]); pause(0);

    end
    Z(j)=tmp(3);
end
% % % figure(2); 
% % % plot3(px,py,z_est_nm,'o'); 
% % % grid on;
% % % hold on; 
% % % zlabel('z / nm');
% % % 
% % % finaldata{img_no}=[px+y_est' py+x_est' BG_est' N_est' z_est_nm'];
% % % 
% % % %for m=1:length(finaldata)
% % % save('loc_data.mat','finaldata','NA','RI','z_vec','ux');


%% data eval
% % % 
% % % xc=90; 
% % % yc=65; 
% % % r_sphere=1e-3; 
% % % r_coord=(1:140)*ux;
% % % z_theory=r_sphere-sqrt(r_sphere^2-r_coord.^2);
% % % figure(4);
% % % 
% % % sig_thresh=4000; %photon threshold; 
% % % for m=1:length(finaldata);
% % %     d=finaldata{m};
% % %     for mm=1:size(d,1); %stepping through all molecules
% % %         mol=d(mm,:);
% % %         if mol(4)>=sig_thresh;
% % %             figure(4); 
% % %             plot(sqrt((mol(1)-(xc)).^2+(mol(2)-(yc)).^2),mol(5),'o'); hold on; 
% % %             ylabel('z /nm');
% % %             %figure(5); 
% % %             %imagesc(mol(1),mol(2),mol(3));
% % %         end
% % %     end
% % % end
% % % plot(r_coord/ux,z_theory*1e9); grid on;
% % % xlabel('radial coord. / pixel');
% % % title('z-localization data');
% % % hold off; 
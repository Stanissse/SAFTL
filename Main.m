clear all
% close all
tic
load('PSF_SAF_5nm.mat')
% Pixelsize in nm
global pxSize DAT middle dist_max PSF_tot
dist_max=600;
pxSize = 115;

Z=[];
X=[];
Y=[];


%% user parameters 

%camera parameters
gain=50; 
amp=9.8; %electrons per count 
QE=0.95; %quantum efficiency

%% Open first file and find point between TOT and UAF img
% Import Fits-Stack
import matlab.io.*
[filenameIMG,PathNameIMG] = uigetfile('*.tif','Select the FITS file');

% Open STACK
info = imfinfo([PathNameIMG,filenameIMG]);
num_images = numel(info);
[x y z]=size(imread([PathNameIMG,filenameIMG], 1, 'Info', info));

imagedata=zeros(x,y,num_images);

for k = 1:num_images
    imagedata(:,:,k) = imread([PathNameIMG,filenameIMG], k, 'Info', info);
end

imagedata = double(imagedata/gain*amp/QE);
[ImagewithoutBG,BG,OrgImg]=EverBG(imagedata);
ImagewithoutBG=ImagewithoutBG-min(ImagewithoutBG(:))+100;


imagedata_=imagedata(:,:,1);
f1=figure('Name','Click in the middle of the two spots');
imagesc(imagedata_)
axis off
[x_coord, y_coord]=ginput(1);
close(f1)
middle=y_coord;

imagedata_=imagedata(:,:,1);
f1=figure('Name','Click in the center of the SAF-Sphere');
imagesc(imagedata_)
axis off
[rx, ry]=ginput(1);
close(f1)

%% SAVE TIFF (If needed)
outputFileName = '-400nm-BG.tif'
for K=1:length(imagedata(1, 1, :))
   imwrite(uint16(ImagewithoutBG(:, :, K)), outputFileName, 'WriteMode', 'append', 'Compression', 'none');
end

%% Ask for single molecule
if AskUser('Do you wnat to look at single molecules? Enter 1 for "yes" or 0 for "no"    ');
    single = 0;
else
    single = 1;
end


%% Select only first element of stack
imgNR=2
if single ==0;
    
    
    imagedata_=imagedata(:,:,imgNR);
    f1=figure('Name','Click at the molecule of interest in the SAF-zone');
    imagesc(imagedata_) 
        % use mouse button to zoom in or out
        % Press Enter to get out of the zoom mode
        zoom on;
        % Wait for the most recent key to become the return/enter key
        waitfor(gcf, 'CurrentCharacter', char(13))
        zoom reset
        zoom off
        axis off
        [PosxSAF, PosySAF]=ginput(1);
        close(f1)
 tic
        R=sqrt((PosxSAF*115-rx*115).^2+(PosySAF*115-ry*115).^2);
        ZestAnal=-sqrt(-(R)^2+(1e6)^2)+1e6;
        posX=int16(PosxSAF);
        posY=int16(PosySAF);
        I=imagedata_(posY-4:posY+5,posX-4:posX+5);
%         figure(22)
%         imagesc(I)
%         title('manual')
        [x,y,z,ii,bg_,ZEst_] = MLE_Estimator_TOT(double(I),25,1600);
toc
end


%% open csv-data (from ThunderStorm)
% Save ThunderStorm-data as: [frame x y sigma intensity]
[filename,PathName] = uigetfile('*.csv','Select the csv file');
DAT = csvread([PathName,filename],1, 0);

%% Calculate ratios
tic
[X, Y, Z_p10, I_, BG,ZEst,Start]=RatioCalculatorMLE(imagedata,rx,ry);
toc


%% delete unphys

% searchs positions of unphysical ratios
Delete = find(ZEst>90);

% sets unphysical ratios to zero
X(Delete) = 0;
Y(Delete) = 0;
ZEst(Delete) = 0;

% delete unphysical ratios
X = X(X~=0);
Y = Y(Y~=0);
ZEst = ZEst(ZEst~=0);

% plot(Start-ZEst,'bx')
%%

% save('Dat_Kugel-500nm_1000.mat','X', 'Y', 'Z_p10', 'I_', 'BG','ZEst')
%%
% Z_p10(Z_p10==0)=[];
% 
% %% delete unphysical values
% [X, Y, Z]=deleteUnPhys(X, Y, Z);
% 
% %%
% Z_p10=Z_p10*uz*1e9;

%% show image
figure(10)
colormap(winter)
% scatter3(X,Y,Z,100,linspace(1,length(Z),length(Z)),'x');
%scatter3(X,Y,Z_p10,100,Z_p10,'x');
%hold on
scatter3(X,Y,ZEst*5,10,linspace(1,length(ZEst),length(ZEst)));
title('MLE-Schätzer (nur TOT)')
toc

%%
R=sqrt((X-(rx+0)*115).^2+(Y-(ry+0)*115).^2);
size(R)
% % select ROI (25000 nm)
Z_=ZEst*5;

Delete = find(R/25000 >= 1);
% sets unphysical ratios to zero
Z_(Delete) = 0;
R(Delete) = 0;


% delete unphysical ratios
R = R(R~=0);
Z_ = Z_(Z_~=0);

figure(99) 
format short g
plot(R/1000,Z_,'bx')
title(['MLE-Schätzer (nur TOT) mean= ' num2str(mean(Z_),3) '(' num2str(std(Z_),2) ')'])
xlabel('Radius [\mum]')
ylabel('z-Abstand [nm]')
hold on
%%#
R_sort=[R/1000,Z_];
R_new=sortrows(R_sort);

% plot(R_new(:,1))

%%
clear NEW

k=2;
i=1;
j=1;

while k<=length(R_new(:,1))
        while (abs(R_new(k-1,1) - R_new(k,1))<.1) && k<length(R_new(:,1))
            NEW(j,i,:)=R_new(k,:);
            k=k+1
            i=i+1;
        end
    k=k+1;
    j=j+1
end


%%
clear STD MEAN

[S_new,~,~]=size(NEW);

for i=1:S_new
    
    STD(i)=std(deleteZeros(NEW(i,:,2)));
    MEAN(i)=mean(deleteZeros(NEW(i,:,2)));
    
end


errorbar(MEAN,STD,'bx')
function m=estimator3D(v)
 x1=v(:,1);
 x2=v(:,2);
 x3=v(:,3);
 x4=v(:,4);
 x5=v(:,5);

 
        global pxSize DAT I PSF_tot
 
%         PSF=PSF_tot; 
% %         clear PSF_tot; 
% 
%         nx0=size(PSF,1); %size of model
%         ny0=nx0; 
% 
%         %PSF: normalization of energy in each z-slice
%         energies=(trapz(trapz(PSF,1),2)); 
%         PSF_norm=PSF./repmat(energies,[size(PSF,1),size(PSF,2),1]);
%         dz_PSF=diff(PSF_norm,1,3); %calculating derivative along z (required for ML estimation)
%         PSF_norm(:,:,end)=[]; %delete last slice to match size of dz_PSF


 

fw=2;

   
% 	 load('PSF_SAF_5nm.mat')
    [~,~,dz]=size(PSF_tot);
    [dx,~,~]=size(I);
    

%     m= sum(sum((x3+x2*interpn(PSF_tot,(fw+1:dx+fw)-x4,(fw+1:dx+fw)'-x5,x1))-I.*log(1e-18+(x3+x2*interpn(PSF_tot,(fw+1:dx+fw)-x4,(fw+1:dx+fw)'-x5,x1)))))+10^5 * any(x2 < 500);
    m = sum(sum(x3+x2*interpn(PSF_tot,(fw+1:dx+fw)-x4,(fw+1:dx+fw)'-x5,x1)-I.*log(x3+x2*interpn(PSF_tot,(fw+1:dx+fw)-x4,(fw+1:dx+fw)'-x5,x1)),1),2);
    
end
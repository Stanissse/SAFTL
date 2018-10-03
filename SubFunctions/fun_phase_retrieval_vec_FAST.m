function [error, I]=fun_phase_retrieval_vec_FAST(Zernikes,Image_stack, kernel, F_kernel, defocus, z_vec, E_BFP, dia_pupil, Zernike_stack, b, display)
%for vectorial phase retrieval:
%function calculates stack of defocused bead images and compares it to a
%measured stack ("Image_stack"); output is an error metric
%Zernikes...vector of Zernike coefs, starting with "tip" (piston is
%omitted)
%Image_stack...measured stack of bead-images at varying values of defocus;
%must be normalized: sum(Image_stack(:))=1
%kernel...convolution kernel for chirped z-trafo
%F_kernel...the Fourier transform of "kernel"
%defocus...phase function for 1m of defocus
%z_vec...contains z-positions of images in Image_stack
%E_BFP...tensor containing back focal plane fields: Ex_Px, Ex_Py, Ex_Pz, Ey_Px, Ey_Py, Ey_Pz... BFP fields
%dia_pupil...size of pupil field in pixel
%Zernike_stack...stack of Zernike modes to be measured
%T_obj...BFP field transmission (objective apodization)
%display...'y' or 'n' if set to 'y', focal plane intensities are displayed
%global Image_stack defocus z_vec ux uk Ex_Px Ex_Py Ex_Pz Ey_Px Ey_Py Ey_Pz dia_pupil Zernike_stack T_obj
%b...3D-FFT of fluorescent bead; used at the very end to low-pass filter
%simulated stack

tic
N_img=size(Image_stack,1); 
no_images=size(Image_stack,3); %number of images in the z-stack
Nx=size(kernel,1)-dia_pupil; %size of simulated focal plane
Nk=dia_pupil; %size of input field (k-space)
%N_pad=2^nextpow2(Nx+Nk-1); %padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
N_pad=Nx+Nk-0;

hugo=zeros(1,1,size(Zernike_stack,3));
hugo(1,1,1:length(Zernikes))=Zernikes; 
pupil_in=exp(1i*sum(Zernike_stack.*repmat(hugo,[dia_pupil,dia_pupil,1]),3));
%pupil_in=exp(1i*sum(ZernikeCalc(1:length(Zernikes),Zernikes,pupil,'Noll')));
% % % imagesc(angle(pupil_in))
conj_kernel=conj(kernel);

for m=1:no_images

    
    FF=pupil_in.*exp(1i*defocus*z_vec(m));
    
    for q=1:size(E_BFP,3) %calculate focal field for every BFP-field component

        %-----CZT 2D--------------------------------
        E_in=E_BFP(:,:,q).*FF;

        f1=embed(E_in,N_pad,0).*conj_kernel;

        %convolving f1 and f2
        hugo=fftshift(ifft2(fft2(ifftshift(f1)).*F_kernel));
        E_out=embed(N_pad*conj_kernel.*hugo,N_img,0); 

        %provisional normalization: intensities are preserved
        E_out=E_out/sum(E_out(:).*conj(E_out(:))).*(sum(E_in(:).*conj(E_in(:))));
        I_out(:,:,q)=E_out.*conj(E_out);
        %---------------------------------------------------
    end
    
    I(:,:,m)=sum(I_out,3); 
        
 
    if strcmp(display,'y');
        subplot(2,1,1);
        imagesc(Image_stack(:,:,m)); title('measured');
        subplot(2,1,2);
        imagesc(I(:,:,m)); title('retrieved'); pause;
    end
end

%----considering effects of bead-size -> low pass filtering -----
I_lp=fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(I))).*b)));

tmp=I_lp/sum(I_lp(:));


% % % weight=[1.2,1.2,1.2,1.2,1.2,1.2,1,1,1,1,1,1.2,1.2,1.2,1.2,1.2,1.2];
% % % 
% % % for i=1:17
% % %    temp(i)=sum(sum(abs(Image_stack(:,:,i)-tmp(:,:,i))))*weight(i);
% % % end


error=sum(abs(Image_stack(:)-tmp(:)))

% 
% error=sum(temp);
toc

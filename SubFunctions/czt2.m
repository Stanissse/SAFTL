function E_out=czt2(E_in,uk,ux_des,Nx)
%2D chirp z-trafo (Tong, Raoqiong, and Robert W. Cox. "Rotation of NMR images using the 2D chirp-z transform." Magnetic resonance in medicine 41.2 (1999): 253-256.
%see also OneNote (AJ)
%E_in...input field
%uk...length unit in E_in
%ux_des...desired object space resolution
%Nx...size of field in pixel
%clear all; 

% %----test parameters--------
%     Nk=64;
%     Nx=64;
%     lambda_0=532e-9;
%     RI=1;
%     NA=0.92;
%     kt_max=2*pi/lambda_0*NA;
%     uk=kt_max/(Nk/2); %k-space unit
% 
%     ux_des=100e-9;
%     [Kx,Ky,Kr,pupil,defocus,~,uk,ux,v]=fun_usualsuspects(Nk,lambda_0,ux_des,NA,RI);
%     E_in=pupil;
% %---------------------------

Nk=size(E_in,1); %size of input field (k-space)
N_pad=2^nextpow2(Nx+Nk-1); %padded size required for convolution is Nx+Nk-1; here we pad more in order to have a gridsize that is a power of 2 (=faster)
N_pad=Nx+Nk-0;

x=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5))';
y=(floor(-N_pad(1)/2+0.5):floor(N_pad(1)/2-0.5));

%E_in=R.^2<100^2;

%-----CZT 2D--------------------------------
    ux_fft=2*pi/(Nk*uk); %fft resolution without padding
    r=ux_des/ux_fft; %required r-factor to meet the desired resolution at the given grid size N
    alpha=r*1/Nk;
  
    %kernel=exp(-1i*alpha*pi*R2); %convolution kernel
    kernel=exp(-1i*alpha*pi*x.^2)*exp(-1i*alpha*pi*y.^2); %create quadratic convolution phase kernel; faster method
    f1=embed(E_in,N_pad,0).*conj(kernel);

    %convolving f1 and f2
    %disp('czt2 time=');
    %tic
    hugo=fftshift(ifft2(fft2(ifftshift(f1)).*fft2(ifftshift(kernel))));
    E_out=embed(N_pad*conj(kernel).*hugo,Nx,0); 
    %toc
    
    
    
%     figure(1);
%     imagesc(abs(E_out)); axis equal; axis tight; title('CZT2'); colorbar;


%----------comparing with FFT2--------------
%     N_pad2=round(2*pi/(ux_des*uk)); %padding required to achieve ux_des
%     disp('FFT2 time=');
%     tmp=ifftshift(embed(E_in,N_pad2,0));
%     tic
%     E=fft2(tmp)/N_pad2; %intensity normalization
%     toc
%     figure(2);
%     E_out2=embed(fftshift(abs(E).^1),Nx,0);
%     imagesc(E_out2); axis equal; axis tight; title('FFT2'); colorbar;
%     ux_fft2=2*pi/((size(E,1))*uk);
% 


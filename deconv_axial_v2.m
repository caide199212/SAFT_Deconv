%function deconv_axial(axial_width)
% deconv_axial(axial_width)
% axial_width e.g. 35e-6   35 um
% 3D deconvolution...step 1: axial
% after 2D SAFT
% set iteration times (10 times)
load SAFT_data
axial_width=36e-6;
%% =====SET PARAMETERS==== [USER]
%axial_width=35e-6;  % 32 um
NUMIT=10;    % iteration times
DAMPAR = 0;     %0.0001 related to SNR
%------PSF-Axial
sgm=axial_width/2.3548;
nz2=50;imz2=(1:nz2)*dimz;
h=exp(-((imz2-round(nz2/2)*dimz)).^2/(2*sgm^2));
h=h/sum(h);  %sum(h)=1
PSF = h;

im3D_dec_axial=zeros(nz,nx,ny);  % data after axial deconvolution
tot=tic;     % total time for calculation
for i=1:nx
    fprintf('------At x=%1.0d------\n',i)
    for j=1:ny
        %fprintf('------At x=%1.0d; y=%1.0d------\n',i,j)
        %==== Lucy deconvolution
        I=squeeze(RF_CF(:,i,j));
        J1 = deconvlucy(I',PSF,NUMIT);  
        im3D_dec_axial(:,i,j)=J1;
        %figure; plot(imz*1e3,I,imz*1e3,J1); %xlim([zc-0.5 zc+0.5])
    end
end
fprintf('Total time is %4.4f min\n',toc(tot)/60)
save dec_axial_data im3D_dec_axial
return


%% -----XZ cross-section and YZ cross-section check
close all
x_mid2=round(nx/2);y_mid2=round(ny/2);
%----Before axial deconvolution
%x_mid=round(Nx/2);y_mid=round(Ny/2);
xz_mid=squeeze(RF_CF(:,y_mid2,:))';
yz_mid=squeeze(RF_CF(x_mid2,:,:))';
%----After axial deconvolution
%xz_mid_RF=squeeze(RF_h(:,y_mid2,:))';
%yz_mid_RF=squeeze(RF_h(x_mid2,:,:))';
xz_mid_axial_deconv=squeeze(im3D_dec_axial(:,y_mid2,:))';
yz_mid_axial_deconv=squeeze(im3D_dec_axial(x_mid2,:,:))';

%-----XZ cross-section
figure(1)
suptitle('XZ cross-section')
subplot(1,2,1)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid));colormap(gray);colorbar;axis image
title('Before axial deconv');xlabel('X(mm)');ylabel('Z(mm)')
subplot(1,2,2)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid_axial_deconv));colormap(gray);colorbar;axis image
title('After axial deconv');xlabel('X(mm)');ylabel('Z(mm)')

%-----YZ cross-section
figure(2)
suptitle('YZ cross-section')
subplot(1,2,1)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid));colormap(gray);colorbar;axis image
title('Before axial deconv');xlabel('Y(mm)');ylabel('Z(mm)')
subplot(1,2,2)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid_axial_deconv));colormap(gray);colorbar;axis image
title('After axial deconv');xlabel('Y(mm)');ylabel('Z(mm)')


%-----XY cross-section  Max Amplituede Projection  section-by-section
%---Original image
% need to adjust section start point manually
% dx=x(2)-x(1);
% dy=y(2)-x(1);
% z_sec=round((op_z(2)-op_z(1))/dz);
no=7; %imz_start=140;
z_id=zeros(2,no);
%z_id(1,:)=[140 240 350 450 560 670 780];
z_id(2,:)=[40 190 340 490 640 790 940];
z_sec2=150;



%----lateral width and SNR(for 7 objects)
lateral_width=zeros(2,no);  % 1: orginal 2:SAFT 3:SAFT+CF
axial_width=zeros(2,no);
SNR_ob=zeros(3,no); %along xy plane
for i=1:no % along object
    %---Before axial deconv
    sig_sec=RF_CF(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max=max(abs(sig_sec),[],3);
    z_max=max(max(xy_max));
    [row col]=find(xy_max==z_max);
    sig_z_max=squeeze(sig_sec(row,col,:));
    axial_width(1,i)=FWHM(sig_z_max)*dimz;  % axial width
    %sig_y=xy_max(:,y_mid);
    sig_y=xy_max(x_mid2,:);                % easy to check one cross-section
    lateral_width(1,i)=FWHM(sig_y)*dimy;    % lateral width
    %---After axial deconv
    sig_sec2=im3D_dec_axial(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max2=max(abs(sig_sec2),[],3);
    z_max2=max(max(xy_max2));
    [row col]=find(xy_max2==z_max2);
    sig_z_max2=squeeze(sig_sec2(row,col,:));
    axial_width(2,i)=FWHM(sig_z_max2)*dimz;  % axial width
    %sig_y2=xy_max2(:,y_mid2);
    sig_y2=xy_max2(x_mid2,:);
    lateral_width(2,i)=FWHM(sig_y2)*dimy;      % lateral width
    
    figure
    suptitle('XY MAP')
    subplot(1,2,1)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max));colormap(gray);colorbar;axis image
    title('Before axial deconv');xlabel('X(mm)');ylabel('Y(mm)')
    subplot(1,2,2)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max2));colormap(gray);colorbar;axis image
    title('After axial deconv');xlabel('X(mm)');ylabel('Y(mm)')
    %     fprintf('Lateral width is %4.3f mm\n',lateral_width*1e3)
    %     figure(3)
    %     plot(x*1e3,y/max(y))
    %     title('Lateral width');xlabel('X(mm)');ylabel('Normalized PA amplitude')
end
%----print axial width
fprintf('Axial width-Before axial deconv \n')
fprintf('%4.3f um \n',axial_width(1,:)*1e6)
fprintf('Axial width-After axial deconv\n')
fprintf('%4.3f um \n',axial_width(2,:)*1e6)
%----print lateral width
fprintf('Lateral width-Before axial deconv \n')
fprintf('%4.3f um \n',lateral_width(1,:)*1e6)
fprintf('Lateral width-After axial deconv\n')
fprintf('%4.3f um \n',lateral_width(2,:)*1e6)

%% 

%save dec_axial_data im3D_dec_axial RF_CF imx imy imz







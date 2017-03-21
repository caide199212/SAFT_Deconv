%function deconv_lateral(lateral_width)
% deconv_lateral(lateral_width)
% lateral_width e.g. 50e-6   50 um
% 3D deconvolution...step 2: lateral
% set iteration times (10 times)

load dec_axial_data
lateral_width=95e-6;
%% =====SET PARAMETERS==== [USER]
%lateral_width=60e-6;  % 50 um
NUMIT=15;    % iteration times
DAMPAR = 0;     %0.0001 related to SNR

%% ------PSF-Lateral
sgm=lateral_width/2.3548;
%hsize=round(sgm/dimx*6/2)*2+1; hsize=min([ny hsize]);
hsize=40;
PSF = fspecial('gaussian',hsize,sgm/dimx);
PSF=PSF/max(max(PSF));
%imagesc(PSF);axis image
%shg;return
im3D_dec_lateral=zeros(nz,nx,ny);  % data after 3D deconvolution both axially and laterally
for i=1:nz
    fprintf('------At z=%1.0d------\n',i)
    %==========II.1 Lucy deconvolution
    I=squeeze(im3D_dec_axial(i,:,:));
   
    J1 = deconvlucy(I,PSF,NUMIT,DAMPAR);  % just choose i1 (or i2)
    im3D_dec_lateral(i,:,:)=J1;
    %plot figure. better not plot to speed up!!!
    %figure
    %subplot(121); imagesc(imx*1e3,imy*1e3,I); axis image; title('Original')
    %subplot(122); imagesc(imx*1e3,imy*1e3,J1); axis image; title('Deconved')
    %h=figure(1);
    %MOV=getframe(fig);
    %aviobj=addframe(aviobj,MOV);
    %F(i)=getframe;
    %pause(0.1)
end
%movie2avi(F,'video_sample.avi','compression','MSVC');
%movie2avi(F,'video_sample.avi','compression','None');
%aviobj=close(aviobj)
save deconv_lateral_data im3D_dec_lateral
return

%% -----XZ cross-section and YZ cross-section check
close all;clc
x_mid2=round(nx/2);y_mid2=round(ny/2);
%----Before axial deconvolution
%x_mid=round(Nx/2);y_mid=round(Ny/2);
xz_mid=squeeze(RF_CF(:,y_mid2,:))';
yz_mid=squeeze(RF_CF(x_mid2,:,:))';
%----After axial deconvolution
xz_mid_axial_deconv=squeeze(im3D_dec_axial(:,y_mid2,:))';
yz_mid_axial_deconv=squeeze(im3D_dec_axial(x_mid2,:,:))';
%----After lateral deconvolution
xz_mid_lateral_deconv=squeeze(im3D_dec_lateral(:,y_mid2,:))';
yz_mid_lateral_deconv=squeeze(im3D_dec_lateral(x_mid2,:,:))';

%-----XZ cross-section
figure(1)
suptitle('XZ cross-section')
subplot(1,3,1)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid));colormap(gray);colorbar;axis image
title('Before axial deconvolution');xlabel('X(mm)');ylabel('Z(mm)')
subplot(1,3,2)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid_axial_deconv));colormap(gray);colorbar;axis image
title('After axial deconvolution');xlabel('X(mm)');ylabel('Z(mm)')
subplot(1,3,3)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid_lateral_deconv));colormap(gray);colorbar;axis image
title('After lateral deconvolution');xlabel('Y(mm)');ylabel('Z(mm)')
%-----YZ cross-section
figure(2)
suptitle('YZ cross-section')
subplot(1,3,1)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid));colormap(gray);colorbar;axis image
title('Before axial deconvolution');xlabel('Y(mm)');ylabel('Z(mm)')
subplot(1,3,2)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid_axial_deconv));colormap(gray);colorbar;axis image
title('After axial deconvolution');xlabel('Y(mm)');ylabel('Z(mm)')
subplot(1,3,3)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid_lateral_deconv));colormap(gray);colorbar;axis image
title('After lateral deconvolution');xlabel('Y(mm)');ylabel('Z(mm)')

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
lateral_width=zeros(3,no);  % 1:Before axial deconvolution 2:After axial deconvolution 3:After lateral deconvolution
axial_width=zeros(3,no);
SNR_ob=zeros(3,no); %along xy plane
for i=1:no % along object
    %---Before axial deconvolution
    sig_sec=RF_CF(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max=max(abs(sig_sec),[],3);
    z_max=max(max(xy_max));
    [row col]=find(xy_max==z_max);
    sig_z_max=squeeze(sig_sec(row,col,:));
    axial_width(1,i)=FWHM(sig_z_max)*dimz;  % axial width
    %sig_y=xy_max(:,y_mid);
    sig_y=xy_max(x_mid2,:);                % easy to check one cross-section
    lateral_width(1,i)=FWHM(sig_y)*dimy;    % lateral width
    %---After axial deconvolution
    sig_sec2=im3D_dec_axial(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max2=max(abs(sig_sec2),[],3);
    z_max2=max(max(xy_max2));
    [row col]=find(xy_max2==z_max2);
    sig_z_max2=squeeze(sig_sec2(row,col,:));
    axial_width(2,i)=FWHM(sig_z_max2)*dimz;  % axial width
    %sig_y2=xy_max2(:,y_mid2);
    sig_y2=xy_max2(x_mid2,:);
    lateral_width(2,i)=FWHM(sig_y2)*dimy;      % lateral width
    %---After lateral deconvolution
    sig_sec3=im3D_dec_lateral(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max3=max(abs(sig_sec3),[],3);
    z_max3=max(max(xy_max3));
    [row col]=find(xy_max3==z_max3);
    sig_z_max3=squeeze(sig_sec3(row(1),col(1),:));
    axial_width(3,i)=FWHM(sig_z_max3)*dimz;  % axial width
    %sig_y3=xy_max3(:,y_mid2);
    sig_y3=xy_max3(x_mid2,:);
    lateral_width(3,i)=FWHM(sig_y3)*dimy;      % lateral width
    
    figure
    suptitle('XY MAP')
    subplot(1,3,1)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max));colormap(gray);colorbar;axis image
    title('Before axial deconvolution');xlabel('X(mm)');ylabel('Y(mm)')
    subplot(1,3,2)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max2));colormap(gray);colorbar;axis image
    title('After axial deconvolution');xlabel('X(mm)');ylabel('Y(mm)')
    subplot(1,3,3)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max3));colormap(gray);colorbar;axis image
    title('After lateral deconvolution');xlabel('X(mm)');ylabel('Y(mm)')
    %     fprintf('Lateral width is %4.3f mm\n',lateral_width*1e3)
    %     figure(3)
    %     plot(x*1e3,y/max(y))
    %     title('Lateral width');xlabel('X(mm)');ylabel('Normalized PA amplitude')
end
%----print axial width
fprintf('Axial width-Before axial deconvolution \n')
fprintf('%4.3f um \n',axial_width(1,:)*1e6)
fprintf('Axial width-After axial deconvolution\n')
fprintf('%4.3f um \n',axial_width(2,:)*1e6)
fprintf('Axial width-After lateral deconvolution\n')
fprintf('%4.3f um \n',axial_width(3,:)*1e6)
%----print lateral width
fprintf('Lateral width-Before axial deconvolution \n')
fprintf('%4.3f um \n',lateral_width(1,:)*1e6)
fprintf('Lateral width-After axial deconvolution\n')
fprintf('%4.3f um \n',lateral_width(2,:)*1e6)
fprintf('Lateral width-After lateral deconvolution\n')
fprintf('%4.3f um \n',lateral_width(3,:)*1e6)

%save deconv_lateral_data im3D_dec_lateral  im3D_dec_axial RF_CF imx imy imz


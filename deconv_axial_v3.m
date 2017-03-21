% 3D deconvolution...step 1: axial
% set iteration times
% after 2D SAFT
% Cai De 2016/01/17
%clc;clear all;close all
%load 'D:\people\Caide\Simulation\SAFT_2D/2d_SIR_circle_SAFT_data2.mat'
close all;clc
nx=length(imx);
ny=41;
%ny=length(imy);
nz=length(imz);
dimx=imx(2)-imx(1); dimy=imy(2)-imy(1);dimz=imz(2)-imz(1);
RF_CF=permute(RF_CF(:,:,20:60),[2 3 1]);
im3D=double(RF_CF);  %RF_CF is single previously
%im3D=RF;
%% =====SET PARAMETERS==== [USER]
c=1480;
axial_width=35e-6;  % 32 um
NUMIT=10;    % iteration times
DAMPAR = 0;     %0.0001 related to SNR
%-----set tolerance
%can tune to smaller value if AMSE not achieved global min,
%but run time could be too long! %Song Liang: AMSE reached local min.
%tol1=1e-2;      % for AMSE
%tol2=1e-2;      % for AER
%% ------1. shift zero level
%------2. reduce nz
%------3. check SNR and apply BPF
%------PSF-Axial
sgm=axial_width/2.3548;
nz2=50;imz2=(1:nz2)*dimz;
h=exp(-((imz2-round(nz2/2)*dimz)).^2/(2*sgm^2));
h=h/sum(h);  %sum(h)=1
PSF = h;

im3D_dec_axial=zeros(nx,ny,nz);  % data after axial deconvolution
tot=tic;     % total time for calculation
for i=1:nx
    fprintf('------At x=%1.0d------\n',i)
    for j=1:ny
        %fprintf('------At x=%1.0d; y=%1.0d------\n',i,j)
        %==========II.1 Lucy deconvolution
        I=squeeze(im3D(i,j,:))';
        %        AMSEall=[]; AERall=[]; AMSE=1;     %AER=1;%AERmin=1;
        %         ttic=tic;
        %         while (AMSE>1e-10)
        %             if NUMIT==1
        %                 J_i_1=I; % Ji-1 the previous one
        %                 luc1_cell=deconvlucy({I},PSF,1,DAMPAR);  %{I} cell array
        %             else
        %                 J_i_1=luc1_cell{2}; % Ji-1 the previous one
        %                 luc1_cell=deconvlucy(luc1_cell,PSF,1,DAMPAR);
        %             end
        %             J_i=luc1_cell{2};   % Ji   the current one
        %             f1=conv2(J_i,PSF,'same');
        %             %------calculate AMSE and AER
        %             AMSE=sum(sum((J_i-J_i_1).^2))/sum(sum(J_i_1.^2));
        %             AER=sum(sum((f1-I).^2))/sum(sum(I.^2));
        %             AMSEall=[AMSEall AMSE]; AERall=[AERall AER];
        %             %plot figure. better not plot to speed up!!!
        %             figure(2); plot(imz*1e3,I,imz*1e3,J_i); %xlim([zc-0.5 zc+0.5])
        %             legend('Original','Deconved')
        %             %------criterion to stop
        %             if NUMIT>=3
        %                 VAR1=abs((AMSE-AMSEall(end-1))/AMSEall(end-1));
        %                 VAR2=abs((AER-AERall(end-1))/AERall(end-1));
        %             end
        %             if NUMIT>20 ...
        %                     && ( (AMSE>min(AMSEall)) || VAR1<tol1 ) ... %STOP at golbal min or variation too small
        %                     && ( (AER>min(AERall)) || VAR2<tol2 )       %STOP at golbal min or variation too small
        %                 [Y,i1]=min(AMSEall);
        %                 [Y,i2]=min(AERall);
        %                 if i1>=20
        %                     fprintf('Total Number of iteration: %d\n',NUMIT)
        %                     fprintf('Number of iteration for min. AMSE: %d\n',i1)
        %                     fprintf('Number of iteration for min. AER: %d\n',i2)
        %                     if VAR1>=tol1
        %                         fprintf('1. AMSE achieved global minimum\n')
        %                     end
        %                     if VAR2>=tol2
        %                         fprintf('2. AER achieved global minimum\n')
        %                     end
        %                     break
        %                 end
        %             end
        %             NUMIT=NUMIT+1;
        %         end
        %         fprintf('Time: %2.2f mins\n\n',toc(ttic)/60)
        J1 = deconvlucy(I,PSF,NUMIT,DAMPAR);  % just choose i1 (or i2)
        im3D_dec_axial(i,j,:)=J1';
        %figure; plot(imz*1e3,I,imz*1e3,J1); %xlim([zc-0.5 zc+0.5])
    end
end
fprintf('Total time is %4.4f min\n',toc(tot)/60)
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







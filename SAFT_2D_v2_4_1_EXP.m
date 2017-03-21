% 2D SAFT---depth-dependent in a circle
% (nt,nx,ny) nt time axis; nx fast axis; ny slow axis
%  nt in first dimension, easy for processing and display
%  loops order, ny first, nx second, nt the last
% SIR correction or sqrt(N)
% focal zone considered separately
% RF use the data before hilbert transform
% RF*CF and then apply hilbert transform
% Cai De 2016/05/23
%clc;clear all;close all

%% load data
load raw_data
%load parameters3
% sigall=permute(sigall,[3 1 2]);     % (nt,Nx,Ny)
% sigall_h=permute(sigall_h,[3 1 2]); % (nt,Nx,Ny)
%load SIR_dream_data sig_DOF imz_SIR   % attention imz_SIR unit [mm]
%close all;clc;clear memory

%% parameters set

Nx=100; % fast axis
Ny=100;  % slow axis
nt=1002;

%step size
dx=10e-6;  % step size of fast axis [m]
dy=10e-6;  % step size of slow axis [m]
dt=2e-9;   % 500 MS/s sampling rate [s]
dz=3e-6;   % step size of time axis 2ns*1.5um/ns=3um; [m]
x=[-round(Nx/2)*dx:dx:(Nx-round(Nx)/2-1)*dx];  % length(x)=Nx;
y=[-round(Ny/2)*dy:dy:(Ny-round(Ny)/2-1)*dy];  % length(y)=Ny;

%% time delay read from scope when laser pulse acts as trigger
% t_delay corresponds to the MIDDLE POINT of Aline,  attention!!!
t_delay=4.4e-6;  % [s] 5.76
t=[t_delay-round(nt/2)*dt:dt:t_delay+(round(nt)/2-1)*dt];
c=1500;  % sound speed [m/s]
z=t*c;


% ---region of interest
% reduced region
Xm=0.2e-3;
Ym=0.2e-3;
imz1=5.7e-3;imz2=7.5e-3;
%step size
dimx=5e-6;  % 20 um
dimy=5e-6;  % 20 um
dimz=3e-6;   % 2 um

imx=-Xm:dimx:Xm;
imy=-Ym:dimy:Ym;
imz=imz1:dimz:imz2;

nx=length(imx);ny=length(imy);nz=length(imz);

% time axis adjustment
% t_shift=50e-6/c; % 100 um shift of focal point in z axis +: upward -: downward
% t=t+t_shift;

%% transducer
c=1500;  % sound speed
f=6.7*1e-3;   %focal length unit: 8.52 mm  8.64 mm 
D=6*1e-3;   %diameter of transducer
Theta=asin(D/2/f);  % polar angle  25бу
Nx=length(x);
Ny=length(y);
[X,Y]=meshgrid(x,y);  % transducer's position matrix
BD=60e-6; %beam diameter in the focal zone for SAFT in the focal zone
nf=40; % focal zone nf*dimz  50

%% 
%----RF and CF
RF=zeros(nz,nx,ny,'single');    % single to save memory
RF_h=zeros(nz,nx,ny,'single'); % RF after hilbert transform
CF=zeros(nz,nx,ny,'single');
RF_CF=zeros(nz,nx,ny,'single');  % RF*CF and then apply hilbert transform
tot=tic;  % total time for calculation
for i=41%1:ny %38 %250:500  % slow axis
    [temp Ny_c]=min(abs(y-imy(i)));    % the closest transducer's x position
    display(num2str(i))
    tic
    parfor j=1:nx % parfor  % fast axis
        %display(num2str(j))
        %[temp Nx_c]=min(abs(x-imx(i)));    % the closest transducer's y position
        [temp Nx_c]=min(abs(x-imx(j)));
        for k=1:nz
            %---2D SAFT
            %---depth-dependent
            [~, Nz_c]=min(abs(t*c-imz(k)));    % the closest z position of signal in Aline
            %[temp Nz_SIR]=min(abs(imz_SIR*1e-3-imz(k))); %
            nu=0;de=0;RF_xy=0;
            dis=sqrt((X-imx(j)).^2+(Y-imy(i)).^2+(f-imz(k))^2);dis=dis';
            dis_z=imz(k)-f;
            if dis_z>0
                dis=f+dis;
            else
                dis=f-dis;
            end
            dis_xy=sqrt((X-imx(j)).^2+(Y-imy(i)).^2);dis_xy=dis_xy';
            hz=abs(imz(k)-f);
            d_xy=hz*sin(Theta);  %  times 2 or 3 doesn't affect the uniform amplitude??
            [row col]=find(dis_xy<=d_xy); %sum in a circle  +2*BD
            N_xy=length(row);
            if abs(dis_z)<=nf*dimz      % in the focal zone
                %[temp t_idx]=min(abs(dis(Nx_c,Ny_c)/c-t));  % choose the nearest one
                %RF_xy=RF_xy+sigall(Nx_c,Ny_c,t_idx);
                %nu=nu+sigall(Nx_c,Ny_c,t_idx);    %numerator
                %de=de+sigall(Nx_c,Ny_c,t_idx)^2;  % denominator of CF
                RF(k,j,i)=sigall(Nz_c,Nx_c,Ny_c);
                CF(k,j,i)=1;%(nu)^2/(1*de);
                %if N_xy==0      % in the focal zone N_xy may be zero
                %[temp t_idx]=min(abs(dis(Nx_c,Ny_c)/c-t));  % choose the nearest one
                %RF_xy=RF_xy+sigall(Nx_c,Ny_c,t_idx);
                %nu=nu+sigall(Nx_c,Ny_c,t_idx);    %numerator
                %de=de+sigall(Nx_c,Ny_c,t_idx)^2;  % denominator of CF
                %RF(i,j,k)=RF_xy*sig_DOF(Nz_SIR);
                %CF(i,j,k)=(nu)^2/(1*de);
            else
                for m=1:N_xy
                    [temp t_idx]=min(abs(dis(row(m),col(m))/c-t));
                    RF_xy=RF_xy+sigall(t_idx,row(m),col(m));
                    %nu=nu+sigall(row(m),col(m),t_idx);    %numerator
                    de=de+sigall(t_idx,row(m),col(m))^2;  % denominator of CF
                end
                nu=RF_xy;
                if N_xy && de
                    RF(k,j,i)=RF_xy/sqrt(N_xy);
                    %RF(i,j,k)=RF_xy*sig_DOF(Nz_SIR);%/sqrt(N_xy);  % divided by N_xy or not?? sqrt(N_xy)???
                    CF(k,j,i)=(nu)^2/(N_xy*de);
                end
            end
        end
%         RF_h(:,j,i)=abs(hilbert(RF(:,j,i)));
%         RF_CF(:,j,i)=abs(hilbert(RF(:,j,i).*CF(:,j,i)));
%         figure
%         temp0=squeeze(sigall_n(38,38,:));
%         plot(temp0)
%         title('Original')
%         figure
%         temp1=squeeze(RF(i,j,:));
%         plot(temp1)
%         title('RF')
%         figure
%         temp2=squeeze(CF(i,j,:));
%         plot(temp2)
%         title('CF')
%         figure
%         plot(abs(hilbert(temp1.*temp2)))
%         title('SAFT+CF')
%         
    end
%     toc
    Bmode_Ori=squeeze(sigall_h(:,:,Ny_c)); % the closest one
    Bmode_RF=squeeze(RF(:,:,i));
    Bmode_CF=squeeze(CF(:,:,i));
    figure
    subplot(1,3,1)
    imagesc(y*1e3,t*1.5*1e6,mat2gray(Bmode_Ori));axis image
    title('Original image');xlabel('Y(mm)');ylabel('Z(mm)')
    subplot(1,3,2)
    imagesc(imx*1e3,imz*1e3,mat2gray(abs(hilbert(Bmode_RF))));axis image
    title('SAFT image');xlabel('Y(mm)');ylabel('Z(mm)')
    subplot(1,3,3)
    imagesc(imx*1e3,imz*1e3,mat2gray(abs(hilbert(Bmode_RF.*Bmode_CF))));axis image
    title('SAFT+CF image');xlabel('Y(mm)');ylabel('Z(mm)')
end
% Hilbert transform
% put outside to avoid problem of parfor
parfor i=1:ny
    for j=1:nx
        RF_h(:,j,i)=abs(hilbert(RF(:,j,i)));
        RF_CF(:,j,i)=abs(hilbert(RF(:,j,i).*CF(:,j,i)));
    end
end


fprintf('Total time is %4.3f min\n',toc(tot)/60)
return
%% -----XZ cross-section and YZ cross-section check
close all
%----Original image
x_mid=round(Nx/2);y_mid=round(Ny/2);
xz_mid=squeeze(sigall_h(:,y_mid,:))';
yz_mid=squeeze(sigall_h(x_mid,:,:))';
%----SAFT image
x_mid2=round(nx/2);y_mid2=round(ny/2);
xz_mid_RF=squeeze(RF_h(:,y_mid2,:))';
yz_mid_RF=squeeze(RF_h(x_mid2,:,:))';
xz_mid_RF_CF=squeeze(RF_CF(:,y_mid2,:))';
yz_mid_RF_CF=squeeze(RF_CF(x_mid2,:,:))';

%-----XZ cross-section
figure(1)
suptitle('XZ cross-section')
subplot(1,3,1)
imagesc(x*1e3,t*1.5*1e6,mat2gray(xz_mid));colormap(gray);colorbar;axis image
title('Original image');xlabel('X(mm)');ylabel('Z(mm)')
subplot(1,3,2)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid_RF));colormap(gray);colorbar;axis image
title('SAFT image');xlabel('X(mm)');ylabel('Z(mm)')
subplot(1,3,3)
imagesc(imx*1e3,imz*1e3,mat2gray(xz_mid_RF_CF));colormap(gray);colorbar;axis image
title('SAFT+CF image');xlabel('Y(mm)');ylabel('Z(mm)')
%-----YZ cross-section
figure(2)
suptitle('YZ cross-section')
subplot(1,3,1)
imagesc(y*1e3,t*1.5*1e6,mat2gray(yz_mid));colormap(gray);colorbar;axis image
title('Original image');xlabel('Y(mm)');ylabel('Z(mm)')
subplot(1,3,2)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid_RF));colormap(gray);colorbar;axis image
title('SAFT image');xlabel('Y(mm)');ylabel('Z(mm)')
subplot(1,3,3)
imagesc(imy*1e3,imz*1e3,mat2gray(yz_mid_RF_CF));colormap(gray);colorbar;axis image
title('SAFT+CF image');xlabel('Y(mm)');ylabel('Z(mm)')

%-----XY cross-section  Max Amplituede Projection  section-by-section
%---Original image
%op_z=op(:,3);
op_z=[5.8 6.1 6.4 6.7 7 7.3 7.6]*1e-3;
dz=c*(t(2)-t(1));
dx=x(2)-x(1);
dy=y(2)-x(1);
z_sec=round((op_z(2)-op_z(1))/dz);
[row col]=find(xz_mid==max(max(xz_mid)));
%z_id=row-round(no/2*z_sec);
% need to adjust section start point manually
z_id=zeros(2,no);
z_id(1,:)=[140 240 350 450 560 670 780];
%---SAFT image
z_id(2,:)=[40 190 340 490 640 790 940];
z_sec2=150;
%z_sec2=round(z_sec*dz/dimz)+1;



%----lateral width and SNR(for 7 objects)
lateral_width=zeros(3,no);  % 1: orginal 2:SAFT 3:SAFT+CF
axial_width=zeros(3,no);
SNR_ob=zeros(3,no); %along xy plane
for i=1:no % along object
    %---Original iamge
    sig_sec=sigall_h(:,:,z_id(1,i):z_id(1,i)+z_sec);
    xy_max=max(abs(sig_sec),[],3);
    z_max=max(max(xy_max));
    [row col]=find(xy_max==z_max);
    sig_z_max=squeeze(sig_sec(row,col,:));
    axial_width(1,i)=FWHM(sig_z_max)*dz;  % axial width
    %sig_y=xy_max(:,y_mid);
    sig_y=xy_max(x_mid,:);                % easy to check one cross-section
    lateral_width(1,i)=FWHM(sig_y)*dy;    % lateral width
    %---SAFT image
    sig_sec2=RF_h(:,:,z_id(2,i):z_id(2,i)+z_sec2);
    xy_max2=max(abs(sig_sec2),[],3);
    z_max2=max(max(xy_max2));
    [row col]=find(xy_max2==z_max2);
    sig_z_max2=squeeze(sig_sec2(row,col,:));
    axial_width(2,i)=FWHM(sig_z_max2)*dimz;  % axial width
    %sig_y2=xy_max2(:,y_mid2);
    sig_y2=xy_max2(x_mid2,:);
    lateral_width(2,i)=FWHM(sig_y2)*dimy;      % lateral width
    %---SAFT+CF image
    sig_sec3=RF_CF(:,:,z_id(2,i):z_id(2,i)+z_sec2);
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
    imagesc(x*1e3,y*1e3,mat2gray(xy_max));colormap(gray);colorbar;axis image
    title('Original image');xlabel('X(mm)');ylabel('Y(mm)')
    subplot(1,3,2)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max2));colormap(gray);colorbar;axis image
    title('SAFT image');xlabel('X(mm)');ylabel('Y(mm)')
    subplot(1,3,3)
    imagesc(imx*1e3,imy*1e3,mat2gray(xy_max3));colormap(gray);colorbar;axis image
    title('SAFT+CF image');xlabel('X(mm)');ylabel('Y(mm)')
    %     fprintf('Lateral width is %4.3f mm\n',lateral_width*1e3)
    %     figure(3)
    %     plot(x*1e3,y/max(y))
    %     title('Lateral width');xlabel('X(mm)');ylabel('Normalized PA amplitude')
end
%----print axial width
fprintf('Axial width-Original \n')
fprintf('%4.3f um \n',axial_width(1,:)*1e6)
fprintf('Axial width-SAFT\n')
fprintf('%4.3f um \n',axial_width(2,:)*1e6)
fprintf('Axial width-SAFT+CF\n')
fprintf('%4.3f um \n',axial_width(3,:)*1e6)
%----print lateral width
fprintf('Lateral width-Original \n')
fprintf('%4.3f um \n',lateral_width(1,:)*1e6)
fprintf('Lateral width-SAFT\n')
fprintf('%4.3f um \n',lateral_width(2,:)*1e6)
fprintf('Lateral width-SAFT+CF\n')
fprintf('%4.3f um \n',lateral_width(3,:)*1e6)

%% save data_0115_02 Bmode_CF Bmode_Ori Bmode_RF Bmode_RFh
%save 2d_circle_SAFT_data2 RF RF_h RF_CF CF sigall sigall_h imx imy imz x y t nx ny Nx Ny lateral_width axial_width SNR_ob
%save 2d_SIR_circle_SAFT_data2 RF RF_h RF_CF CF sigall sigall_h imx imy imz x y t nx ny Nx Ny lateral_width axial_width SNR_ob


function SAFT_2D(option)
% SAFT_2D(option)
% two options 'whole' 'middle'
% middle calculate the middle slice for check
% 2D SAFT---depth-dependent in a circle
% (nt,nx,ny) nt time axis; nx fast axis; ny slow axis
% loops order  ny first, nx second, nt the last
% SIR correction or sqrt(N)
% focal zone considered separately
% RF use the data before hilbert transform
% RF*CF and then apply hilbert transform

%% load data
load raw_data

if strcmp(option,'whole')
    y_range=[1 ny];
elseif strcmp(option,'middle')
    y_range=round((ny+1)/2);
else
    disp('Reset option, ''whole'', or ''middle''')
    return
end
    

%load SIR_dream_data sig_DOF imz_SIR   % attention imz_SIR unit [mm]
%close all;clc;clear memory

% time axis adjustment
% t_shift=50e-6/c; % 100 um shift of focal point in z axis +: upward -: downward
% t=t+t_shift;

%% transducer
c=1500;  % sound speed
f=6.7*1e-3;   %focal length unit: 8.52 mm  8.64 mm 
D=6*1e-3;   %diameter of transducer
Theta=asin(D/2/f);  % half polar angle  
[X,Y]=meshgrid(x,y);  % transducer's position matrix
Z=k_x*X+k_y*Y;          % tilted Z position
% figure
% mesh(Z)
BD=60e-6; %beam diameter in the focal zone for SAFT in the focal zone
nf=40; % focal zone nf*dimz  50

%% 
%----RF and CF
RF=zeros(nz,nx,ny,'single');    % single to save memory
RF_h=zeros(nz,nx,ny,'single'); % RF after hilbert transform
CF=zeros(nz,nx,ny,'single');
RF_CF=zeros(nz,nx,ny,'single');  % RF*CF and then apply hilbert transform
tot=tic;  % total time for calculation
figure
for i=y_range %38 %250:500  % slow axis
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
            %[temp Nz_SIR]=min(abs(imz_SIR*1e-3-imz(k))); 
            nu=0;de=0;RF_xy=0;
            dis=sqrt((X-imx(j)).^2+(Y-imy(i)).^2+((Z+f)-imz(k))^2);dis=dis';
            dis_z=imz(k)-f;  %(Z+f) f is OK for out-of-focus zone, but maybe some difference for focal zone
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
    toc
    Bmode_Ori=squeeze(sigall_h(:,:,Ny_c)); % the closest one
    Bmode_RF=squeeze(RF(:,:,i));
    Bmode_CF=squeeze(CF(:,:,i));
    %figure
    subplot(1,3,1)
    imagesc(x*1e3,t*1.5*1e6,mat2gray(Bmode_Ori));axis image
    title('Original image');xlabel('X (mm)');ylabel('Z (mm)')
    subplot(1,3,2)
    imagesc(imx*1e3,imz*1e3,mat2gray(abs(hilbert(Bmode_RF))));axis image
    title('SAFT image');xlabel('X (mm)');ylabel('Z (mm)')
    subplot(1,3,3)
    imagesc(imx*1e3,imz*1e3,mat2gray(abs(hilbert(Bmode_RF.*Bmode_CF))));axis image
    title('SAFT+CF image');xlabel('X (mm)');ylabel('Z (mm)')
end
% Hilbert transform
% put outside to avoid problem of parfor
for i=1:ny
    for j=1:nx
        RF_h(:,j,i)=abs(hilbert(RF(:,j,i)));
        RF_CF(:,j,i)=abs(hilbert(RF(:,j,i).*CF(:,j,i)));
    end
end

fprintf('Total time is %4.3f min\n',toc(tot)/60)
save SAFT_data RF RF_h CF RF_CF
return

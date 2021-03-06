% 2D SAFT+3D Deconv
% Cai De 2016/07/17
% nz(nt) in first dimension, easy for processing and display
% ====Set Parameters====
% nt:time axis(z axis)  Nx:fast axis  Ny:slow axis
nt=1002; Nx=100; Ny=100;

%----detector position----
dx=10e-6;  % step size of fast axis [m]
dy=10e-6;  % step size of slow axis [m]
x=[-round(Nx/2)*dx:dx:(Nx-round(Nx)/2-1)*dx];  % length(x)=Nx;
y=[-round(Ny/2)*dy:dy:(Ny-round(Ny)/2-1)*dy];  % length(y)=Ny;
k_x=-0.0; k_y=0.0;  % tilted angle with respect to x axis and y axis  z axis is downward
% Start time point of Aline for SAFT 
% time delay read from scope when laser pulse acts as trigger
% t_delay corresponds to the MIDDLE POINT of Aline, attention!!!
dt=2e-9;   % 500 MS/s sampling rate [s]
t_delay=4.34e-6;  % [s]
t=[t_delay-round(nt/2)*dt:dt:t_delay+(round(nt)/2-1)*dt];
c=1500;  % sound speed [m/s]
dz=dt*c;   % step size of time axis 2ns*1.5um/ns=3um; [m]
z=t*c;
% ---region of interest---
% reduced region
Xm=0.2e-3;
Ym=0.2e-3;
imz1=5.7e-3;imz2=7.5e-3;
% pixel size
dimx=5e-6;  % 20 um
dimy=5e-6;  % 20 um
dimz=3e-6;   % 2 um

imx=-Xm:dimx:Xm;
imy=-Ym:dimy:Ym;
imz=imz1:dimz:imz2;

nx=length(imx);ny=length(imy);nz=length(imz);

%% ====data preprocessing====

% re-arrange Bmode data as sigall(nt,Nx,Ny) and sigall_h(nt,Nx,Ny)
% save raw_data sigall sigall_h

% re_arrange(nt,Nx,Ny)
load raw_data
save raw_data % include the parameters above
% Aline bandwidth check

%----Optional processing----
% Remove DC shift
% Matched filter

%return
%% ====2D SAFT==== 
% written as function cannot use parfor
SAFT_2D('middle')
return

%% ====interpolation====

%% ====axial deconvolution====
deconv_axial(35e-6)
return

%% ====lateral deconvolution====
% depth-independent PSF
deconv_lateral(60e-6)









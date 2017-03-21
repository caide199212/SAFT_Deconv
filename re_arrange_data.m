function re_arrange_data(nt,Nx,Ny)
% re-arrange_data(nt,Nx,Ny)
% nt:time axis(z axis)  Nx:fast axis  Ny:slow axis
% save raw_data sigall(nt,Nx,Ny) and sigall_h(nt,Nx,Ny)
% e.g. nt=1002; Nx=100; Ny=100;

sigall=zeros(nt,Nx,Ny);   % original signal
sigall_h=zeros(nt,Nx,Ny); % original signal after hilbert transform
for i=1:Ny
    FileName=['Bmode' num2str(i) '.mat'];
    Bmode=importdata(FileName);
    Bmode_h=abs(hilbert(Bmode));
    sigall(:,:,i)=Bmode;
    sigall_h(:,:,i)=Bmode_h;
end
save raw_data sigall sigall_h
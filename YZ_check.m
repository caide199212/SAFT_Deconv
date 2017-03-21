function YZ_check(data_3D)
% YZ_check(3D_data)
% [nz(nt) nx ny]
[nz nx ny]=size(data_3D);
figure
for i=1:nx
    imagesc(squeeze(data_3D(:,i,:)))
    pause(0.2)
end
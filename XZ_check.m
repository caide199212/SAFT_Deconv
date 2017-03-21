function XZ_check(data_3D)
% % XZ_check(data)
% [nz(nt) nx ny]
% Bmode check
[nz nx ny]=size(data_3D);
figure
for i=1:ny
    imagesc(squeeze(data_3D(:,:,i)))
    pause(0.2)
end

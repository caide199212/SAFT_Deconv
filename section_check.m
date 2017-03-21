function section_check(data_3D,sec1_z,sec2_z,sec3_z,sec4_z,sec5_z)
% section_check(data_3D,sec1_z,sec2_z,sec3_z,sec4_z,sec5_z)
% e.g. sect1_z=102:125
% [nz(nt) nx ny]
% check objects section by section
sec1=data_3D(sec1_z,:,:);
sec2=data_3D(sec2_z,:,:);
sec3=data_3D(sec3_z,:,:);
sec4=data_3D(sec4_z,:,:);
sec5=data_3D(sec5_z,:,:);
xy_max1=squeeze(max(sec1,[],1));
xy_max2=squeeze(max(sec2,[],1));
xy_max3=squeeze(max(sec3,[],1));
xy_max4=squeeze(max(sec4,[],1));
xy_max5=squeeze(max(sec5,[],1));
% y1=xy_max1(:,21);
% y2=xy_max2(:,21);
% y3=xy_max3(:,21);
% y4=xy_max4(:,21);
% y5=xy_max5(:,21);

figure;imagesc(xy_max1')
figure;imagesc(xy_max2')
figure;imagesc(xy_max3')
figure;imagesc(xy_max4')
figure;imagesc(xy_max5')
% figure;plot(y1)
% figure;plot(y2)
% figure;plot(y3)
% figure;plot(y4)
% figure;plot(y5)
save xy_max xy_max1 xy_max2 xy_max3 xy_max4 xy_max5



figure;imagesc(abs(hilbert(Bmode1)));shg
figure;imagesc(abs(hilbert(Bmode2)));shg
figure;imagesc(abs(hilbert(Bmode3)));shg
figure;imagesc(abs(hilbert(Bmode4)));shg
figure;imagesc(abs(hilbert(Bmode5)));shg

Bmode=[Bmode1(1:410,:);Bmode2(411:513,:);Bmode3(514:613,:);Bmode4(614:727,:);Bmode5(728:end,:)];

figure;imagesc(abs(hilbert(Bmode)));shg
Bmode_h=abs(hilbert(Bmode));
sigall=repmat(Bmode,[1 1 100]);
sigall_h=repmat(Bmode_h,[1 1 100]);
save raw_data sigall sigall_h

y1=max(abs(hilbert(Bmode1)));
y2=max(abs(hilbert(Bmode2)));
y3=max(abs(hilbert(Bmode3)));
y4=max(abs(hilbert(Bmode4)));
y5=max(abs(hilbert(Bmode5)));
y=[y1;y2;y3;y4;y5];
figure;plot(1:200,y)

figure;imagesc(Bmode_RF_CF)
y1=max(Bmode_RF_CF(107:137,:));
y2=max(Bmode_RF_CF(208:256,:));
y3=max(Bmode_RF_CF(298:366,:));
y4=max(Bmode_RF_CF(408:457,:));
y5=max(Bmode_RF_CF(517:545,:));
figure;plot(y1);shg
figure;plot(y2);shg
figure;plot(y3);shg
figure;plot(y4);shg
figure;plot(y5);shg


sec1=Bmode_RF_CF(107:137,:);
sec1=Bmode_RF_CF(208:256,:);
sec1=Bmode_RF_CF(107:137,:);
sec2=Bmode_RF_CF(208:256,:);
sec3=Bmode_RF_CF(298:366,:);
sec4=Bmode_RF_CF(408:457,:);
sec5=Bmode_RF_CF(517:545,:);
Aline1=sec1(:,53);
Aline2=sec1(:,52);
Aline3=sec1(:,34);
Aline4=sec1(:,35);
Aline5=sec1(:,52);
figure;plot(Aline1);shg
figure;plot(Aline2);shg
figure;plot(Aline3);shg
figure;plot(Aline4);shg
figure;plot(Aline5);shg

Bmode_lateral=squeeze(im3D_dec_lateral(:,:,41));
sec1=Bmode_lateral(113:137,:);
sec2=Bmode_lateral(208:256,:);
sec3=Bmode_lateral(298:366,:);
sec4=Bmode_lateral(408:457,:);
sec5=Bmode_lateral(517:545,:);
y1=max(sec1);
y2=max(sec1);
y2=max(sec2);
y3=max(sec3);
y4=max(sec4);
y5=max(sec5);
 figure;plot(y1);shg
figure;plot(y2);shg
figure;plot(y3);shg
figure;plot(y4);shg
figure;plot(y5);shg


Aline1=sec1(:,53);
Aline2=sec2(:,52);
Aline3=sec3(:,34);
Aline4=sec4(:,35);
Aline5=sec5(:,52);
figure;plot(Aline1);shg
figure;plot(Aline2);shg
figure;plot(Aline3);shg
figure;plot(Aline4);shg
figure;plot(Aline5);shg






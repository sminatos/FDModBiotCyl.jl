clear all
close all

org=load("~/work/Kiguchi_Nojima_DSI/VSP_Tubewave/Modeling/White_model/matome_Three_Mechanism/Permeable_Layer/THICK/out_frame2B_back2_thick.mat");

now=load("./out_frame2B_back2_thick.mat");



figure;imagesc(org.rec_vz);colorbar
figure;imagesc(now.rec_vz);colorbar

figure;plot(org.rec_vz(:,400),'r-');
hold on;plot(now.rec_vz(:,400),'b.');


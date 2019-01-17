%% An example for the tool.
%   Author:         Xiaochen Qiu from Beihang Univ.
%   Description:    Show how to use the tool by an example of:
%                   Evaluate trajectory position RMSE for a monocular VIO 
%                   on sequence V2_03_difficult from EuRoC dataset.

clc;
clear;
close all;

addpath('DATAS','FUNCTIONS');

%% load ground truth and estimate trajectory
load gt.csv;    % provided by EuRoC
load est.txt;

%% need to adjust the time standard of ground truth
take_off_stamp = 1413394886.455760384;  % timestamp standard in est.txt
gt(:,1) = gt(:,1)*1e-9;
gt(:,1) = gt(:,1)-take_off_stamp;

%% extract timestamp and 3d position from 'gt' and 'est', these are all we need for alignment
time_gt = gt(:,1);
P_gt = gt(:,2:4)';
time_es = est(:,1);
P_es = est(:,9:11)';

%% alignment
[Ids_es, Ids_gt] = findIds (time_es, time_gt, 0.001);
[R_es, t_es, s_es] = sim3DataAlignment (P_es(:,Ids_es), P_gt(:,Ids_gt));

%% do not miss '/s_es' to maintain the estimated scale
P_es_aligned = R_es*P_es + repmat(t_es,1,size(P_es,2))/s_es;

%% draw trajectory
time_matched = time_es(Ids_es);
P_es_matched = P_es_aligned(:,Ids_es);
P_gt_matched = P_gt(:,Ids_gt);
figure;
plot3(P_gt_matched(1,:),P_gt_matched(2,:),P_gt_matched(3,:),'r-');
hold on;
plot3(P_es_matched(1,:),P_es_matched(2,:),P_es_matched(3,:),'b-');
for i = 1:size(P_gt_matched,2)  % draw difference, this block is time consuming
    plot3([P_gt_matched(1,i),P_es_matched(1,i)],...
        [P_gt_matched(2,i),P_es_matched(2,i)],[P_gt_matched(3,i),P_es_matched(3,i)],'y-');
end
axis equal;
grid on;
lgd = legend('ground truth','estimated');
set(lgd,'Fontname','Times New Roman','FontWeight','bold','FontSize',15);
title('V203','FontSize',15);

%% draw error in XYZ axes
errX = P_es_aligned(1,Ids_es)-P_gt(1,Ids_gt);
errY = P_es_aligned(2,Ids_es)-P_gt(2,Ids_gt);
errZ = P_es_aligned(3,Ids_es)-P_gt(3,Ids_gt);
figure;
subplot(311);
plot(errX);
title('position error of axis-X (m)');
subplot(312);
plot(errY);
title('position error of axis-Y (m)');
subplot(313);
plot(errZ);
title('position error of axis-Z (m)');
fprintf('mean error in [X Y Z]: [%fm %fm %fm]\n',mean(errX),mean(errY),mean(errZ));

%% some printing
errVec = P_es_aligned(1:3,Ids_es)-P_gt(1:3,Ids_gt);
N = size(errVec,2);
RMSE_trans = 0;
for i = 1:N
    RMSE_trans = RMSE_trans+norm(errVec(:,i))^2;
end
RMSE_trans = sqrt(RMSE_trans/N);
fprintf('RMSE of translation is %fm\n',RMSE_trans);

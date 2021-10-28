reset_toolbox;
close all

%% ���ɺ�����
% generate HydraMarker

% shape{1,1} = ones(3);
% shape{2,1} = ones(10,1);
% shape{3,1} = ones(6,2);
% [sta,img] = build_marker(10,10,200,shape);

%% ���߶�ȡ�ѱ���ĺ�����
% or load a HydraMar ker
load 10x10_for_3x3_6x2_10x1.mat

%% ��ȡһ�Ű����������ͼƬ
% read an image containing HydraMarker
% img = im2double(rgb2gray(imread('t1.bmp')));
img = im2double(rgb2gray(imread('curve1.jpg')));
% img = imresize(img,720/max(size(img,[1,2])));

%% ʶ�������е�������
% identify the features of HydraMarker
expectN = 2*(size(sta,1)+1)*(size(sta,2)+1);
[ptList,edge] = read_marker(img,sta,7,expectN,3);

%% ��ʾ
% display
figure;
imshow(img);
hold on;
% ���Ʊ� draw edges
Y = ptList(:,1);
X = ptList(:,2);
plot(X(edge'),Y(edge'),'LineWidth',3,'Color','g');
% ���Ƶ� draw dots
scatter(ptList(:,2),ptList(:,1),50,'g','filled','o','LineWidth',1);
% ���Ʋ�ȷ��ID�ĵ� draw unsure IDs
pt_uID = ptList(isnan(ptList(:,3)),:);
scatter(pt_uID(:,2),pt_uID(:,1),100,'r','x','LineWidth',3);
% ����ID draw IDs
pt_ID = ptList(~isnan(ptList(:,3)),:);
text(pt_ID(:,2),pt_ID(:,1),num2str(pt_ID(:,3)),'FontSize',10,'Color','y');

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
% load 10x10_for_3x3_6x2_10x1.mat
load 6x18_for3x3_6x2_10x1.mat

%% ��ȡһ�Ű����������ͼƬ
% read an image containing HydraMarker
% img = im2double(rgb2gray(imread('curve_6.jpg')));
img = im2double((imread('F:\Marker_curved-surface\Experiment\SyntheticCorner\RandomCornerNoise\Noise0.003\4.bmp')));
img_min = min(img(:));
img_max = max(img(:));
img     = (img-img_min)/(img_max-img_min);
% img = imresize(img,2000/max(size(img,[1,2])));

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
scatter(ptList(:,2),ptList(:,1),20,'b','filled','o','LineWidth',1);
% ���Ʋ�ȷ��ID�ĵ� draw unsure IDs
pt_uID = ptList(isnan(ptList(:,3)),:);
scatter(pt_uID(:,2),pt_uID(:,1),100,'r','x','LineWidth',3);
% ����ID draw IDs
pt_ID = ptList(~isnan(ptList(:,3)),:);
text(pt_ID(:,2)+3,pt_ID(:,1)-15,num2str(pt_ID(:,3)),'FontSize',10,'Color','y');

reset_toolbox;
close all

%% 读取已保存的海拉码
% load a HydraMarker
load 10x10_for_3x3_6x2_10x1.mat

%% 读取包含海拉码的视频
% read a video containing HydraMarker
vidObj = VideoReader('marker4.mp4');
vid = VideoWriter('marker_read');
open(vid);
fig = figure;
while hasFrame(vidObj)
    
    img = im2double(rgb2gray(readFrame(vidObj)));
%     img = imresize(img,720/max(size(img,[1,2])));
    img = imnoise(img);
    
    %% 识别海拉码中的特征点
    % identify the features of HydraMarker
    expectN = 2*(size(sta,1)+1)*(size(sta,2)+1);
    [ptList,edge] = read_marker(img,sta,5,expectN*1.5,3);
    
    %% 显示
    % display
    hold off;
    imshow(img);
    hold on;
    % 绘制边 draw edges
    Y = ptList(:,1);
    X = ptList(:,2);
    plot(X(edge'),Y(edge'),'LineWidth',3,'Color','g');
    % 绘制点 draw dots
    scatter(ptList(:,2),ptList(:,1),100,'g','filled','o','LineWidth',1);
    % 绘制不确定ID的点 draw unsure IDs
    pt_uID = ptList(isnan(ptList(:,3)),:);
    scatter(pt_uID(:,2),pt_uID(:,1),100,'r','x','LineWidth',3);
    % 绘制ID draw IDs
    pt_ID = ptList(~isnan(ptList(:,3)),:);
    text(pt_ID(:,2),pt_ID(:,1),num2str(pt_ID(:,3)),'FontSize',15,'Color','y');
    
    pause(0.01);
    
    img_full = frame2im(getframe(fig));
    
    writeVideo(vid,img_full);
end

close(vid);

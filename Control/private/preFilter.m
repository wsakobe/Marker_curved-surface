% 从图像"img"中读取特征点候选, 以列表"ptList"输出.
% Read feature candidates from "img" and output them in "ptList".

% - ptList
%   Nx2的矩阵, 每行对应一个特征点候选在图像中的行、列坐标.
%   Nx2 matrix, each row corresponds to the row-column coordinates of a
%   feature candidate.

% - img
%   输入图像必须是单通道double型.
%   Input image must be single channel, double type.

% - r
%   交叉点的范围(边长为2xr+1的正方形), 本函数基于此排除靠近边界的点.
%   The area of a cross point is a square with side-length "2xr+1".
%   The points near boundaries are rejected based on this value.

% - expectN
%   检出候选点的上限.
%   The upper bound of the number of candidates.

% - sigma
%   检测交叉点的扩展范围, 见文章.
%   the propagating range in cross point detection, see paper.

function ptList = preFilter(img,r,expectN,sigma)

    % "G" 表示在向量和中相互抵消的梯度的总模量, 在对比度高且对称的图案(交叉点)附近较高.
    % "G" is the power of the cancelled out gradients, which is large near
    % high contrast and symmetric patterns (cross point).
    [Gx,Gy] = imgradientxy(img);
    
    Gpow = imgaussfilt((Gx.^2+Gy.^2).^0.5,sigma);
    Gsum = (imgaussfilt(Gx,sigma).^2 + imgaussfilt(Gy,sigma).^2).^0.5;
    G = Gpow-Gsum;
    % 非极大值抑制
    % Non-maximum suppression    
    G(imdilate(G,strel('square',3))~=G)=0;%取3x3邻域中最大值 判断当前值是否为最大值
    G(G(:,:)>0.8)=1;
    figure
    imshow(G)
    % 挑选"expectN"个"G"值最高的点
    % Pick the candidates with top-"expectN" "G" value
    G(1:r,:) = 0; G(end-r+1:end,:) = 0;
    G(:,1:r) = 0; G(:,end-r+1:end) = 0;
    G_sort = G(:); G_sort(G_sort<0.1)=[];
    G_sort = sort(G_sort,'descend');
    [im,in] = ind2sub(size(G),find(G>=G_sort(min(expectN,size(G_sort,1)))));
    ptList = [im,in];
    
end
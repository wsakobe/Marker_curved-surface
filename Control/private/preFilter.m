% ��ͼ��"img"�ж�ȡ�������ѡ, ���б�"ptList"���.
% Read feature candidates from "img" and output them in "ptList".

% - ptList
%   Nx2�ľ���, ÿ�ж�Ӧһ���������ѡ��ͼ���е��С�������.
%   Nx2 matrix, each row corresponds to the row-column coordinates of a
%   feature candidate.

% - img
%   ����ͼ������ǵ�ͨ��double��.
%   Input image must be single channel, double type.

% - r
%   �����ķ�Χ(�߳�Ϊ2xr+1��������), ���������ڴ��ų������߽�ĵ�.
%   The area of a cross point is a square with side-length "2xr+1".
%   The points near boundaries are rejected based on this value.

% - expectN
%   �����ѡ�������.
%   The upper bound of the number of candidates.

% - sigma
%   ��⽻������չ��Χ, ������.
%   the propagating range in cross point detection, see paper.

function ptList = preFilter(img,r,expectN,sigma)

    % "G" ��ʾ�����������໥�������ݶȵ���ģ��, �ڶԱȶȸ��ҶԳƵ�ͼ��(�����)�����ϸ�.
    % "G" is the power of the cancelled out gradients, which is large near
    % high contrast and symmetric patterns (cross point).
    [Gx,Gy] = imgradientxy(img);
    
    Gpow = imgaussfilt((Gx.^2+Gy.^2).^0.5,sigma);
    Gsum = (imgaussfilt(Gx,sigma).^2 + imgaussfilt(Gy,sigma).^2).^0.5;
    G = Gpow-Gsum;
    % �Ǽ���ֵ����
    % Non-maximum suppression    
    G(imdilate(G,strel('square',3))~=G)=0;%ȡ3x3���������ֵ �жϵ�ǰֵ�Ƿ�Ϊ���ֵ
    G(G(:,:)>0.8)=1;
    figure
    imshow(G)
    % ��ѡ"expectN"��"G"ֵ��ߵĵ�
    % Pick the candidates with top-"expectN" "G" value
    G(1:r,:) = 0; G(end-r+1:end,:) = 0;
    G(:,1:r) = 0; G(:,end-r+1:end) = 0;
    G_sort = G(:); G_sort(G_sort<0.1)=[];
    G_sort = sort(G_sort,'descend');
    [im,in] = ind2sub(size(G),find(G>=G_sort(min(expectN,size(G_sort,1)))));
    ptList = [im,in];
    
end
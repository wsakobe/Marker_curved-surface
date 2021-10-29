% 从图像"img"中读取海拉码, 输出点列表"ptList"与连接关系"edge".
% Read HydraMarkers from "img", output feature list "ptList" and adjacency 
% list "edge".

% - ptList
%   Nx3的矩阵, 每行对应一个特征点, 各列分别为特征点在图像中的亚像素行、列坐标,
%   以及ID(特征点在本海拉码中的索引).
%   Nx3 matrix, each row corresponds to a feature, the first two columns
%   represent their subpixel coordinates in "img" (row-column), the last
%   column represents their IDs (indexes in HydraMarker).

% - edge
%   Nx2的矩阵, 每行对应一个连接, 包含被连接两点在"ptList"中的索引.
%   一对特征点被连接, 说明函数认为该对点在海拉码中相邻.
%   Nx2 matrix, each row corresponds to a connection between two features, 
%   containing the indexes of the two connected features in "ptList".
%   A pair of features are connected only if they are identified to be
%   neighbors in HydraMarker.

% - img
%   输入图像必须是单通道double型.
%   Input image must be single channel, double type.

% - sta
%   图"img"中所使用的海拉码的点矩阵.
%   The dot matrix of the HydraMarker used in "img".

% - r, expectN, sigma
%   见相关函数
%   see corresponding functions

function [ptList,edge] = read_marker(img,sta,r,expectN,sigma)
    
    if ~exist('r','var')
        r = 5;
    end
    if ~exist('expectN','var')
        expectN = 100;
    end
    if ~exist('sigma','var')
        sigma = 3;
    end
    
    ptList = preFilter(img,r,expectN,sigma);
    [ptList,ledge] = ptRefine(img,ptList,r);
    [array,edge] = ptStruct(ptList,ledge);
    ptList = ptCurvedSurface(img,ptList,array);
    ptList = ptIdentify(img,sta,array,ptList);
    
end


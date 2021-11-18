% 优化"list", 使之具备亚像素精度, 排除假阳性检测, 并获得脊(交叉点的两条边)的朝向.
% Refine "list" to achieve sub-pixel accuracy, false positive rejection,
% and meanwhile estimate the direction of ledges (the two edges of cross point).

% - ptList
%   Nx2的矩阵, 优化后的"list".
%   Nx2 matrix, refined "list".

% - ledgeList
%   Nx2的矩阵, 每行对应一个特征点, 记录其跳白/黑线的角度(坐标系右半部分, 顺时针由黑转白/白转黑的边缘的角度).
%   Nx2 matrix, each row corresponds to the white/black jump angles of a
%   feature (on the right half of coordinate, the angle that pixel becomes
%   white/black, clock-wise).

% - r
%   交叉点的范围(边长为2xr+1的正方形), 将参考范围内的图像用于优化交叉点信息, 见文章.
%   The area of a cross point is a square with side-length "2xr+1".
%   The cross point is refined based on the image content in this area, see
%   paper.

function [ptList,ledgeList] = ptRefine(img,list,r)
   
    % 基于梯度优化交叉点位置
    % refine cross point locations based on gradient
%     [Gx,Gy] = imgradientxy(img);
%     for iter = 1 : 2
%         for it = 1 : size(list,1)
%             im = round(list(it,1));
%             in = round(list(it,2));
% 
%             [M,N] = ndgrid(im-r:im+r,in-r:in+r);
%             Gm = imgaussfilt(Gy(im-r:im+r,in-r:in+r),1);
%             Gn = imgaussfilt(Gx(im-r:im+r,in-r:in+r),1);
%             G = [Gm(:),Gn(:)];
%             p = sum([M(:),N(:)].*G,2);
%             list(it,:) = (G\p)'; %最小二乘利用梯度得亚像素坐标
%         end
%         % 清除靠近边界的点
%         % remove the detected points near boundaries
        illegal = (any(list<r+2,2) | list(:,1)>size(img,1)-r-2 | list(:,2)>size(img,2)-r-2);
        list(illegal,:)=[];
%     end

    % 基于超平面模型（二元二次方程）计算脊的朝向
    % 基于脊的朝向生成标准模板, 计算相关度以排除假阳性检测
    % calculate the ledge angles based on hyperplane model (binary
    % quadratic equation)
    % build standard template based on the ledge angles, and calculate
    % correlation score to reject false positive
    ledge = zeros(size(list,1),2);
    corr = zeros(size(list,1),1);
    angBias = zeros(size(list,1),1);
    
    [u,v] = meshgrid(-r:r);
    ut = u(:);  vt = v(:);
    A = [ut.^2,ut.*vt,vt.^2,ones(size(ut,1),1)];
    for it = 1 : size(list,1)
        iy = round(list(it,1));
        ix = round(list(it,2));
        [X,Y] = meshgrid(ix-r-1:ix+r+1,iy-r-1:iy+r+1);
        b = interp2(X,Y,img(iy-r-1:iy+r+1,ix-r-1:ix+r+1),...
                    list(it,2)+u,list(it,1)+v);
        c = A\b(:);
        theta = rad2deg(atan(roots([c(3),c(2),c(1)])));        
        
        if ~isreal(theta)
            ledge(it,:) = nan(1,2);
            continue;
        end
        
        ledge(it,:) = theta;
        template = sign(A(:,1:3)*c(1:3));
        corr(it) = corr2(template,b(:));
        angBias(it) = rad2deg(abs(angdiff(deg2rad([theta(1),theta(2)]))));
        angBias(it) = abs(angBias(it)-90);
        
        % 区分跳白和跳黑的ledge, 这取决于"k"的符号.
        % Identify the white/black jump angle.
        % It depends on the signum of "k".
        % "k"为正时, "theta"中较大的值为跳白线, 反之则为跳黑线. 
        % When "k" is positive, the larger one in "theta" is the white jump
        % angle, otherwise, the smaller one is.
        sign_k = sign(theta(1)*theta(2)*c(1));
        if sign_k == 1
            ledge(it,:) = [max(theta),min(theta)];
        else
            ledge(it,:) = [min(theta),max(theta)];
        end
        
    end
    % 删除不足2条脊、低相关度和脊夹角过低的点.
    % Remove the detected points with ledges less than 2, low correlation
    % score or sharp included angle.
    idx = isnan(ledge(:,1)) ...
        | corr < max(corr)-0.38 ...
        | angBias > 30;
    list(idx,:) = [];
    ledge(idx,:) = [];

    % 拼合相互距离较小的点
    % Merge the points closed to each other
    if size(list,1)<2
        ptList = list;
        ledgeList = ledge;
        return;
    end
    ptTree = linkage(list,'average','chebychev');
    ptIdx = cluster(ptTree,'Cutoff',2,'Criterion','distance');
    
    ptList = zeros(length(unique(ptIdx)),2);
    ledgeList = zeros(length(unique(ptIdx)),2);
    for it = 1 : length(unique(ptIdx))
       ptList(it,:) = mean(list(ptIdx==it,:),1);
       ledgeList(it,:) = mean(ledge(ptIdx==it,:),1);
    end
    figure
    imshow(img);
    hold on
    scatter(ptList(:,2),ptList(:,1),20,'g','filled');
end
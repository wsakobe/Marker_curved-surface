% 基于给定的参数生成海拉码并给出其点矩阵.
% Generate a HydraMarker with its dot matrix based on parameters.

% - sta
%   海拉码的点矩阵, 指示各方块内反色点的有无.
%   The dot matrix of HydraMaker, indicating the existence of inverse dot
%   in each square.

% - img
%   海拉码的图像.
%   The image of HydraMarker.

% - M,N,L
%   海拉码点矩阵的预期行、列, 以及图像中每个方块的像素边长.
%   Expected row and column numbers of the dot matrix, and the pixel length
%   of the side of the square in image.

% - shape
%   见相关函数
%   see corresponding functions

function [sta,img] = build_marker(M,N,L,shape)

    % 初始化
    % Initialization
    sta = nan(M,N);
    sta(1,1) = 0;
    
    init_pool = build_pool(shape);
    
    % 记录坍缩过程的二叉树
    % The tree recording collapsing process
    ind = 1;
    tree = [zeros(M*N,1),ones(M*N,2)];
    tree(1,:) = [ind,0,1];
    tree_ptr = 1;
        
    while 1
        % 显示当前图案
        % display the current pattern
        img = sta;
        img(isnan(sta)) = 0.5;
        imshow(imresize(img,30,'nearest'));
        pause(0.1);
        
        %% 处理上一次坍缩
        % handle the effect of last collapsing
        % 更新"cloud"
        % update "cloud"
        [cloud,~,need_trace_back] = update_cloud(sta,init_pool);
        if need_trace_back
            trace_condition = 'repetitive_key';
            sta(tree(tree_ptr,1)) = NaN;
            tree_ptr = tree_ptr - 1;
            continue;
        end
        
        % 如果"cloud.size"中存在0, 说明上次坍缩导致了某些元素无状态可选, 需要回溯.
        % 如果"cloud.size"全为NaN, 说明"sta"完全确定, 坍缩完毕.
        % if "cloud.size" has 0, some elements can be neither 0 nor 1 due
        % to last collapsing, roll back.
        % if "cloud.size" are all NaN, "sta" is fully collapsed, end the
        % loop.
        if min(cloud.size(:)) == 0
            trace_condition = 'unable_to_collapse';
            sta(tree(tree_ptr,1)) = NaN;
            tree_ptr = tree_ptr - 1;
            continue;
        elseif isnan(min(cloud.size(:)))
            break;
        end
        
        %% 开始下一次坍缩
        % next collapsing
        tree_ptr = tree_ptr + 1;
        % 优先坍缩只有一种选择的元素
        % 导致只有一种选择的情况有两种: "pool"仅支持一种取值; 已经尝试了另一种取值.
        % collapse the elements who has only one choice
        % It might be caused by: "pool" only support one; another one has been
        % tried.
        fix0 = find(cloud.pos==0);
        fix1 = find(cloud.pos==1);
        if tree(tree_ptr,1)~=0
            ind = tree(tree_ptr,1);
            if tree(tree_ptr,2) == 1 && cloud.pos(ind)<1
                value = 0;
            elseif tree(tree_ptr,3) == 1 && cloud.pos(ind)>0
                value = 1;
            else
                sta(tree(tree_ptr,1)) = NaN;
                sta(tree(tree_ptr-1,1)) = NaN;
                tree(tree_ptr,:) = [0,1,1];
                tree_ptr = tree_ptr - 2;
                continue;
            end
        elseif ~isempty(fix0)
            ind = fix0(1);
            value = 0;
        elseif ~isempty(fix1)
            ind = fix1(1);
            value = 1;
        else
            % 找到最危险的元素
            %（需要优先坍缩的元素，即坍缩其它位置可能导致其无状态可选）
            % find the elements with least choices 
            minSize = min(cloud.size(:));
            ind = find(cloud.size==minSize);
            % 在最危险的像素之中选择一个失衡最严重的元素
            % find the one with most biased possibilities
            imba = nan(M,N);
            imba(ind) = abs(cloud.pos(ind)-0.5);
            ind = find(imba==max(imba(:)));
            ind = ind(1);
            % 坍缩至概率较大的取值
            % collapse to the value with larger possibility in the first
            % time
            if cloud.pos(ind) <= 0.5
                value = 0;
            else
                value = 1;
            end
        end
        
        % 更新"tree"
        % update "tree"
        sta(ind) = value;
        tree(tree_ptr,1) = ind;
        tree(tree_ptr,value+2) = 0;
 
        %% 显示进度
        % display
        settle_num = sum(~isnan(sta(:)));
        fprintf("\n%0.2f%%",100*settle_num/(M*N));
        
    end
        
    % 生成图像
    % generate image
    img = dot2img(sta,L,0);
    
    % 不需要, 但可用于检查自识别单元的唯一性.
    % It is not necessary to generate the keys, but it can be used to
    % examine the uniqueness of self-identifying entires.
    % key = dot2key(sta,init_pool);
    
end
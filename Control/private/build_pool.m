% 基于给定的自识别单元形状, 生成所有可行的自识别单元实例.
% Generate self-identifying entries based on the preset shapes.

% - shape
%   "shape"是一个cell, 每个cell单元描述了一种自识别单元的形状.
%   每个cell单元包含一个矩阵, 尺寸与自识别单元的最小外接矩形相同.
%   其通过1和0取值描述自识别单元的形状: 1-属于自识别单元; 0-不属于自识别单元.
%   例如, [0 1 0;1 1 1;0 1 0]描述了一种十字形的自识别单元.
%   "shape" is a cell. Each cell describes the shape of a self-identifying
%   entry. It is represented as a matrix with 0 and 1, which mean "does not
%   belong to the entry" and "belong to the entry" respectively.
%   For example, [0 1 0;1 1 1;0 1 0] describes a cross shape entry.

% - pool
%   "pool"是一个cell, 其单元行数与"shape"相同, 但包含5列.
%   每一行对应一种形状的自识别单元.
%   第5列指示此种形状的旋转重复性质, 例如,
%       [1,2,3,4]表示第1-4列的实例拥有相同的形状, 存在解读角度的混淆.
%           (发生于"shape"90°旋转对称时)
%       [1,3;2,4]表示第1、3列的实例和第2、4列的实例分别拥有相同的形状, 内部存在解读角度的混淆, 互相之间则不存在.
%           (发生于"shape"180°旋转对称时)
%       [1;2;3;4]表示第1-4列的实例的形状均不同, 不存在解读角度的混淆.
%       第1-4列记录了此种自识别单元的所有可行实例, 其中第2-4列分别是第1列单元90°、180°、270°旋转后的结果.
%       不属于自识别单元的元素在"pool"中为NaN.
%   "pool" is a cell, which has the same row with "shape", but has 5
%   columns. Each one corresponds to a shape of self-identifying entry.
%   The 5-th column indicates how the shapes repeat themselves.
%       [1,2,3,4] means the entries in 1-4 columns have the same shape,
%       which might confuse the reading algorithm (it happens when the
%       entry repeats itself after 90° rotation).
%       [1,3;2,4] means the entries in 1 and 3 columns, 2 and 4 columns
%       have the same shape (it happens when the entry repeats itself after
%       at least 180° rotation).
%       [1;2;3;4] means the entries in the four columns all have different
%       shapes.
%       The 1-4 columns record all the possible entries of corresponding
%       shape, where the entries in the 2-4 columns are the
%       90°/180°/270°-rotated versions of the ones in the first column.
%       The elements that do not belong to the shape is NaN.

function pool = build_pool(shape)
    
    pool = cell(size(shape,1),5);
    
    %% 找出所有可行实例
    % Find out all the possible entries
    for t = 1 : size(shape,1)

        [M,N] = size(shape{t});
        valid_ind = find(shape{t}==1);
        valid_length = size(valid_ind(:),1);
            % 有效元素（不为0）的个数
            % the number of valid elements
        valid_code = (double(dec2bin(0:2^valid_length-1))-48)';
            % char与double的0和1在转换时需要加减48（见ASCII码表）
        code = -ones(M*N,size(valid_code,2));
        code(valid_ind,:) = valid_code;
        
        pool{t,1} = reshape(code,[M,N,size(valid_code,2)]);
    end
    
    %% 判断旋转重复的性质
    % identify how the shapes repeat themselves
    for t = 1 : size(shape,1)
        [shape_M,shape_N] = size(shape{t});
        if (shape_M == shape_N) && all(all(shape{1}==rot90(shape{1})))
            pool{t,5} = [1,2,3,4];
        elseif all(all(shape{1}==rot90(shape{1},2)))
            pool{t,5} = [1,3;2,4];
        else
            pool{t,5} = [1;2;3;4];
        end   
    end
    
    %% 移除具有旋转歧义和旋转混淆的实例
    % remove the entries that confuse with themselves/others
    for t = 1 : size(shape,1)
        [inst_M,inst_N,inst_num] = size(pool{t,1});
        inst = reshape(pool{t,1},[inst_M*inst_N,inst_num])';
        confuse_ind = zeros(inst_num,3);
        
        if size(pool{t,5},1) == 1
            k_iter = 1 : 3;
        elseif  size(pool{t,5},1) == 2
            k_iter = 2;
        else
            k_iter = [];
        end
        
        for k = k_iter
            inst_r = reshape(rot90(pool{t,1},k),[inst_M*inst_N,inst_num])';
            [~,confuse_ind(:,k)] = ismember(inst,inst_r,'rows');
        end
        
        confuse_ind(confuse_ind < repmat((1:inst_num)',[1,3])) = NaN;
        
        confuse_ind = confuse_ind(:);
        confuse_ind(isnan(confuse_ind)) = [];
        pool{t,1}(:,:,confuse_ind) = [];
        
        pool{t,1}(pool{t,1}==-1) = NaN;
            % 为了比对自识别单元的形状, 在移除歧义和混淆之前, 使用-1替代NaN
            % 因为"ismember([1 NaN],[1 NaN],'rows')=0".
            % NaN is replaced by -1 to compare the shapes, because
            % "ismember([1 NaN],[1 NaN],'rows')=0".
    end
    
    %% 填充第2-4列
    % fill the 2-4 columns
    for t = 1 : size(shape,1)
        pool{t,2} = rot90(pool{t,1});
        pool{t,3} = rot90(pool{t,2});
        pool{t,4} = rot90(pool{t,3});
    end
end


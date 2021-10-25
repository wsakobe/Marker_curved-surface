% ���ڸ�������ʶ��Ԫ��״, �������п��е���ʶ��Ԫʵ��.
% Generate self-identifying entries based on the preset shapes.

% - shape
%   "shape"��һ��cell, ÿ��cell��Ԫ������һ����ʶ��Ԫ����״.
%   ÿ��cell��Ԫ����һ������, �ߴ�����ʶ��Ԫ����С��Ӿ�����ͬ.
%   ��ͨ��1��0ȡֵ������ʶ��Ԫ����״: 1-������ʶ��Ԫ; 0-��������ʶ��Ԫ.
%   ����, [0 1 0;1 1 1;0 1 0]������һ��ʮ���ε���ʶ��Ԫ.
%   "shape" is a cell. Each cell describes the shape of a self-identifying
%   entry. It is represented as a matrix with 0 and 1, which mean "does not
%   belong to the entry" and "belong to the entry" respectively.
%   For example, [0 1 0;1 1 1;0 1 0] describes a cross shape entry.

% - pool
%   "pool"��һ��cell, �䵥Ԫ������"shape"��ͬ, ������5��.
%   ÿһ�ж�Ӧһ����״����ʶ��Ԫ.
%   ��5��ָʾ������״����ת�ظ�����, ����,
%       [1,2,3,4]��ʾ��1-4�е�ʵ��ӵ����ͬ����״, ���ڽ���ǶȵĻ���.
%           (������"shape"90����ת�Գ�ʱ)
%       [1,3;2,4]��ʾ��1��3�е�ʵ���͵�2��4�е�ʵ���ֱ�ӵ����ͬ����״, �ڲ����ڽ���ǶȵĻ���, ����֮���򲻴���.
%           (������"shape"180����ת�Գ�ʱ)
%       [1;2;3;4]��ʾ��1-4�е�ʵ������״����ͬ, �����ڽ���ǶȵĻ���.
%       ��1-4�м�¼�˴�����ʶ��Ԫ�����п���ʵ��, ���е�2-4�зֱ��ǵ�1�е�Ԫ90�㡢180�㡢270����ת��Ľ��.
%       ��������ʶ��Ԫ��Ԫ����"pool"��ΪNaN.
%   "pool" is a cell, which has the same row with "shape", but has 5
%   columns. Each one corresponds to a shape of self-identifying entry.
%   The 5-th column indicates how the shapes repeat themselves.
%       [1,2,3,4] means the entries in 1-4 columns have the same shape,
%       which might confuse the reading algorithm (it happens when the
%       entry repeats itself after 90�� rotation).
%       [1,3;2,4] means the entries in 1 and 3 columns, 2 and 4 columns
%       have the same shape (it happens when the entry repeats itself after
%       at least 180�� rotation).
%       [1;2;3;4] means the entries in the four columns all have different
%       shapes.
%       The 1-4 columns record all the possible entries of corresponding
%       shape, where the entries in the 2-4 columns are the
%       90��/180��/270��-rotated versions of the ones in the first column.
%       The elements that do not belong to the shape is NaN.

function pool = build_pool(shape)
    
    pool = cell(size(shape,1),5);
    
    %% �ҳ����п���ʵ��
    % Find out all the possible entries
    for t = 1 : size(shape,1)

        [M,N] = size(shape{t});
        valid_ind = find(shape{t}==1);
        valid_length = size(valid_ind(:),1);
            % ��ЧԪ�أ���Ϊ0���ĸ���
            % the number of valid elements
        valid_code = (double(dec2bin(0:2^valid_length-1))-48)';
            % char��double��0��1��ת��ʱ��Ҫ�Ӽ�48����ASCII���
        code = -ones(M*N,size(valid_code,2));
        code(valid_ind,:) = valid_code;
        
        pool{t,1} = reshape(code,[M,N,size(valid_code,2)]);
    end
    
    %% �ж���ת�ظ�������
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
    
    %% �Ƴ�������ת�������ת������ʵ��
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
            % Ϊ�˱ȶ���ʶ��Ԫ����״, ���Ƴ�����ͻ���֮ǰ, ʹ��-1���NaN
            % ��Ϊ"ismember([1 NaN],[1 NaN],'rows')=0".
            % NaN is replaced by -1 to compare the shapes, because
            % "ismember([1 NaN],[1 NaN],'rows')=0".
    end
    
    %% ����2-4��
    % fill the 2-4 columns
    for t = 1 : size(shape,1)
        pool{t,2} = rot90(pool{t,1});
        pool{t,3} = rot90(pool{t,2});
        pool{t,4} = rot90(pool{t,3});
    end
end


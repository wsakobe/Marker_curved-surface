% ���ڸ����Ĳ������ɺ����벢����������.
% Generate a HydraMarker with its dot matrix based on parameters.

% - sta
%   ������ĵ����, ָʾ�������ڷ�ɫ�������.
%   The dot matrix of HydraMaker, indicating the existence of inverse dot
%   in each square.

% - img
%   �������ͼ��.
%   The image of HydraMarker.

% - M,N,L
%   �����������Ԥ���С���, �Լ�ͼ����ÿ����������ر߳�.
%   Expected row and column numbers of the dot matrix, and the pixel length
%   of the side of the square in image.

% - shape
%   ����غ���
%   see corresponding functions

function [sta,img] = build_marker(M,N,L,shape)

    % ��ʼ��
    % Initialization
    sta = nan(M,N);
    sta(1,1) = 0;
    
    init_pool = build_pool(shape);
    
    % ��¼̮�����̵Ķ�����
    % The tree recording collapsing process
    ind = 1;
    tree = [zeros(M*N,1),ones(M*N,2)];
    tree(1,:) = [ind,0,1];
    tree_ptr = 1;
        
    while 1
        % ��ʾ��ǰͼ��
        % display the current pattern
        img = sta;
        img(isnan(sta)) = 0.5;
        imshow(imresize(img,30,'nearest'));
        pause(0.1);
        
        %% ������һ��̮��
        % handle the effect of last collapsing
        % ����"cloud"
        % update "cloud"
        [cloud,~,need_trace_back] = update_cloud(sta,init_pool);
        if need_trace_back
            trace_condition = 'repetitive_key';
            sta(tree(tree_ptr,1)) = NaN;
            tree_ptr = tree_ptr - 1;
            continue;
        end
        
        % ���"cloud.size"�д���0, ˵���ϴ�̮��������ĳЩԪ����״̬��ѡ, ��Ҫ����.
        % ���"cloud.size"ȫΪNaN, ˵��"sta"��ȫȷ��, ̮�����.
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
        
        %% ��ʼ��һ��̮��
        % next collapsing
        tree_ptr = tree_ptr + 1;
        % ����̮��ֻ��һ��ѡ���Ԫ��
        % ����ֻ��һ��ѡ������������: "pool"��֧��һ��ȡֵ; �Ѿ���������һ��ȡֵ.
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
            % �ҵ���Σ�յ�Ԫ��
            %����Ҫ����̮����Ԫ�أ���̮������λ�ÿ��ܵ�������״̬��ѡ��
            % find the elements with least choices 
            minSize = min(cloud.size(:));
            ind = find(cloud.size==minSize);
            % ����Σ�յ�����֮��ѡ��һ��ʧ�������ص�Ԫ��
            % find the one with most biased possibilities
            imba = nan(M,N);
            imba(ind) = abs(cloud.pos(ind)-0.5);
            ind = find(imba==max(imba(:)));
            ind = ind(1);
            % ̮�������ʽϴ��ȡֵ
            % collapse to the value with larger possibility in the first
            % time
            if cloud.pos(ind) <= 0.5
                value = 0;
            else
                value = 1;
            end
        end
        
        % ����"tree"
        % update "tree"
        sta(ind) = value;
        tree(tree_ptr,1) = ind;
        tree(tree_ptr,value+2) = 0;
 
        %% ��ʾ����
        % display
        settle_num = sum(~isnan(sta(:)));
        fprintf("\n%0.2f%%",100*settle_num/(M*N));
        
    end
        
    % ����ͼ��
    % generate image
    img = dot2img(sta,L,0);
    
    % ����Ҫ, �������ڼ����ʶ��Ԫ��Ψһ��.
    % It is not necessary to generate the keys, but it can be used to
    % examine the uniqueness of self-identifying entires.
    % key = dot2key(sta,init_pool);
    
end
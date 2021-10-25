% ���ڵ�ǰͼ��"sta"����"cloud"��״̬.
% Update "cloud" based on the current "sta"

% - sta
%   "sta"��һ������, ����Ԫ�ش�������״̬: 1-��; 0-��; NaN-����.
%   "sta" is a matrix with 3 states: 1-has a dot; 0-has not; NaN-unset.

% - cloud
%   "cloud"���"sta"��δ��Ԫ�صĿ�����, ����2������:
%   "cloud.pos"��ʾ��δ��Ԫ��ȡֵΪ1�Ŀ�����(���������͡�λ�õ���ʶ��Ԫ֧��ȡ1�ı���);
%   "cloud.size"��ʾ��δ��Ԫ�ر���ʶ��Ԫ֧��ȡ0��1������(��ͬ���͡�λ�õ���ʶ��Ԫ֮��ȡ��Сֵ);
%   �����Ѷ�Ԫ��, "cloud"�����ȡֵ��ΪNaN;
%   �����޷�ȡֵ��Ԫ��, "cloud.pos"��ȡֵΪNaN, "cloud.size"��ȡֵΪ0.
%   "cloud" describes the possibilities of the unset states in "sta":
%   "cloud.pos" means the possibilities that the states being 1 (supported by
%   the entries from different poses);
%   "cloud.size" means the support number of 0 and 1;
%   For settled states, "cloud.pos" and "cloud.size" are all NaN;
%   For the states that be neither 0 nor 1, "cloud.pos" is NaN,
%   "cloud.size" is 0.

% - need_trace_back
%   boolֵ, ĳ����ʶ��Ԫ���������������λ��ʱΪtrue, ˵����һ��̮���ƻ���Ψһ��
%   bool value, it becomes true when an entry hits more than one position,
%   meaning that the uniqueness of self-identifying entires are broken.

function [cloud,pool,need_trace_back] = update_cloud(sta,init_pool)
    
    [M,N] = size(sta);
    need_trace_back = false;
    
    %% ���ѳ��ֵ�ʵ����"init_pool"���Ƴ�, �γ�"pool"
    % remove appeared entries from "init_pool" to form "pool"
    pool = init_pool;
    hit_map = hit(sta,init_pool,false);
    for t = 1 : size(hit_map,1)
        hit_time = zeros(size(hit_map{t,1},3),4);
        hit_time(:,1) = permute(sum(hit_map{t,1},[1,2],'omitnan'),[3,1,2]);
        hit_time(:,2) = permute(sum(hit_map{t,2},[1,2],'omitnan'),[3,1,2]);
        hit_time(:,3) = permute(sum(hit_map{t,3},[1,2],'omitnan'),[3,1,2]);
        hit_time(:,4) = permute(sum(hit_map{t,4},[1,2],'omitnan'),[3,1,2]);
        hit_time = sum(hit_time,2);
        
        if any(hit_time(:)>1)
            % ĳ����ʶ��Ԫ���������������λ��, ˵����һ��̮���ƻ���Ψһ��, ��Ҫ����
            % An entry hits more than one position, meaning that the
            % uniqueness of self-identifying entires are broken.

            cloud = [];
            pool = [];
            need_trace_back = true;
            return;
        end
        % ��������Ŀ���ʵ����"init_pool"���Ƴ�, ͬ�е�ʵ��Ҳ��Ҫ���Ƴ�.
        % Removed the entries that fully hit "sta".
        % The entries in the same rows are removed too.
        pool{t,1}(:,:,hit_time>0) = [];
        pool{t,2}(:,:,hit_time>0) = [];
        pool{t,3}(:,:,hit_time>0) = [];
        pool{t,4}(:,:,hit_time>0) = [];
    end
    
    %% ʹ��ʣ���"pool"����"cloud"
    % update "cloud" by the remaining entries in "pool"
    % "acc_p0"��"acc_p1"�Ǹ�Ԫ�ص��ۻ�ȡֵ����
    %   ��ĳһԪ�ؿ��ܱ���ͬλ�á���ͬ���͵���ʶ��Ԫģ�����, ����ȡֵ���������໥����, ͨ���۳˼��㣩
    % "acc0"��"acc1"�Ǹ�Ԫ��ȡֵΪ0��1ʱ, �ܿ���ʵ��֧�ֵ���Ŀ
    %   ��ĳһԪ�ؿ��ܱ���ͬλ�á���ͬ���͵���ʶ��Ԫģ�����, ����Ϊ����ʵ��֧���������ֵ��
    % "acc_p0"/"acc_p1" are the accumulated possibilities of the states
    % being 0/1 (an element might be included by multiple entries from
    % different poses. The possibilities are independent with each other,
    % thus they are multipled.)
    % "acc0"/"acc1" are the supported numbers of the states being 0/1 (an
    % element might be included by multiple entries. The supported number
    % equals the minimum).
    
    acc_p0 = ones(M,N);    
    acc_p1 = ones(M,N);
    acc0 = realmax*ones(M,N);
    acc1 = realmax*ones(M,N);
    
    hit_map = hit(sta,pool,true);
    
    for t = 1 : size(pool,1)
        % ���ݻ���λ����һ�۳˴�����ʶ��Ԫģ�������ȡֵ����.
        % δ���е�λ�������ֿ�����:
        % ����������Ԫ�ؾ���ȷ��; ������δȷ����Ԫ�ؼȲ���ȡ0Ҳ����ȡ1.
        % Accumulate the possibilities based on the hit results.
        % There are two situations when a region is not hit by any entry:
        % all the elements in this region is settled;
        % some elements in this region can be neither 0 nor 1.
        for shape_m = 1 : size(pool{t,5},1)
            cur_shape = pool{t,5}(shape_m,1);
            [inst_M,inst_N,~] = size(pool{t,cur_shape});
            for acc_m = 1 : size(hit_map{t,cur_shape},1)
                for acc_n = 1 : size(hit_map{t,cur_shape},2)
                    count0 = 0;
                    count1 = 0;
                    for shape_n = 1 : size(pool{t,5},2)
                        cur_dir = pool{t,5}(shape_m,shape_n);
                        hit_inst = pool{t,cur_dir}(:,:,logical(hit_map{t,cur_dir}(acc_m,acc_n,:)));
                        count0 = count0 + sum((hit_inst)~=1,3);
                        count1 = count1 + sum((hit_inst)~=0,3);
                    end
                    
                    temp_p0 = count0./(count0+count1);
                    temp_p1 = count1./(count0+count1);
                    
                    unconcerned = isnan(pool{t,cur_shape}(:,:,1));
                    count0(unconcerned) = realmax;
                    count1(unconcerned) = realmax;
                    temp_p0(unconcerned) = 1;
                    temp_p1(unconcerned) = 1;
                    
                    m_range = acc_m:acc_m+inst_M-1;
                    n_range = acc_n:acc_n+inst_N-1;
                    acc_p0(m_range,n_range) = acc_p0(m_range,n_range) .* temp_p0;
                    acc_p1(m_range,n_range) = acc_p1(m_range,n_range) .* temp_p1;
                    acc0(m_range,n_range) = min(acc0(m_range,n_range),count0);
                    acc1(m_range,n_range) = min(acc1(m_range,n_range),count1);
                end
            end
        end
    end % for t = 1 : size(pool,1)

    acc_p0(~isnan(sta)) = NaN;
    acc_p1(~isnan(sta)) = NaN;
    acc0(~isnan(sta)) = NaN;
    acc1(~isnan(sta)) = NaN;
    
    cloud.pos = acc_p1./(acc_p0+acc_p1);
    cloud.size = acc0 + acc1;

end
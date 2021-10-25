% ʹ��"pool"�е�ʵ��hitͼ��"sta", ��û���ͼ"hit_map".
% Hit "sta" by the entries in "pool", output hit map "hit_map".

% - nanflag
%   ָʾ"sta"�е�NaN�Ƿ�ͬʱ����0��1.
%   ��Ϊ"true", ������"sta"�ж�ӦԪ�ص�ȡֵ���, ����Ϊ����.
%   ��Ϊ"false", ��ֻ��"pool"�ж�ӦԪ�ص�ȡֵͬΪNaN, ����Ϊ����.
%   Indicates whether the NaN in "sta" are treated as both 0 and 1.
%   If "true", these elements are always hit.
%   If "false", these elements are hit only if the corresponding elements
%   in the entries are NaN too.

% - hit_map
%   "hit_map"��һ��cell, ӵ�е�cell��Ԫ������"pool"��ͬ, ��û�е�5��.
%   ÿ��cell��Ԫ��һ������, �ߴ���"sta"��ͬ, ��������1��ά��, �ֱ��ʾ���ĸ�ʵ��������Ԫ��.
%   0-δ����; 1-����; NaN-������(���Դ�Ԫ��Ϊ���Ͻ�, ���޷�������ʶ��Ԫ).
%   "hit_map" is a cell with the same row number of "pool", but without the
%   5-th column.
%   Each cell is a matrix with the same size of "sta", but it has multiple
%   slices to indicate which entry hits the element.
%   0-miss; 1-hit; NaN-meaningless (entries are not fitted in this position).

function hit_map = hit(sta,pool,nanflag)
    
    [M,N] = size(sta);
    hit_map = cell(size(pool,1),4);
    for t = 1 : size(pool,1)
        
        for k = 1 : 4
            [inst_M,inst_N,inst_num] = size(pool{t,k});
            hit_map{t,k} = zeros(M,N,inst_num);
            % "hit_map{t,2}(m,n,i)=1"����"pool{t,2}"�ĵ�i��ʵ��������"sta"����Ԫ��"sta(m,n)"Ϊ���Ͻ�Ԫ�ص�ͬ�ߴ�����.
            % "hit_map{t,2}(m,n,i)=1" means the i-th entry of "pool{t,2}"
            % hits the region of "sta" that start from the element "sta(m,n)" (top-left).
            
            center_m = floor((size(pool{t,k},1)+1)/2);
            center_n = floor((size(pool{t,k},2)+1)/2);
            for i = 1 : inst_num
                SE = pool{t,k}(:,:,i);
                hit_map{t,k}(:,:,i) = (imerode(sta==1|(isnan(sta)&nanflag),SE==1) ...
                                     & imerode(sta==0|(isnan(sta)&nanflag),SE==0));
            end  
            % �Ƴ���������Ϊ��ʶ��Ԫ���ĵı߿�����
            % remove the regions near boundaries
            hit_map{t,k}(1:center_m-1,:,:) = [];  hit_map{t,k}(end-inst_M+center_m+1:end,:,:) = [];
            hit_map{t,k}(:,1:center_n-1,:) = [];  hit_map{t,k}(:,end-inst_N+center_n+1:end,:) = [];
        end
        
    end % end for t = 1 : size(pool,1)
    
end


% 使用"pool"中的实例hit图案"sta", 获得击中图"hit_map".
% Hit "sta" by the entries in "pool", output hit map "hit_map".

% - nanflag
%   指示"sta"中的NaN是否同时视作0和1.
%   若为"true", 则无论"sta"中对应元素的取值如何, 都视为击中.
%   若为"false", 则只有"pool"中对应元素的取值同为NaN, 才视为击中.
%   Indicates whether the NaN in "sta" are treated as both 0 and 1.
%   If "true", these elements are always hit.
%   If "false", these elements are hit only if the corresponding elements
%   in the entries are NaN too.

% - hit_map
%   "hit_map"是一个cell, 拥有的cell单元行数与"pool"相同, 但没有第5列.
%   每个cell单元是一个矩阵, 尺寸与"sta"相同, 但额外有1个维度, 分别表示是哪个实例击中了元素.
%   0-未击中; 1-击中; NaN-无意义(若以此元素为左上角, 则无法容纳自识别单元).
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
            % "hit_map{t,2}(m,n,i)=1"表明"pool{t,2}"的第i个实例击中了"sta"中以元素"sta(m,n)"为左上角元素的同尺寸区域.
            % "hit_map{t,2}(m,n,i)=1" means the i-th entry of "pool{t,2}"
            % hits the region of "sta" that start from the element "sta(m,n)" (top-left).
            
            center_m = floor((size(pool{t,k},1)+1)/2);
            center_n = floor((size(pool{t,k},2)+1)/2);
            for i = 1 : inst_num
                SE = pool{t,k}(:,:,i);
                hit_map{t,k}(:,:,i) = (imerode(sta==1|(isnan(sta)&nanflag),SE==1) ...
                                     & imerode(sta==0|(isnan(sta)&nanflag),SE==0));
            end  
            % 移除不可能作为自识别单元中心的边框区域
            % remove the regions near boundaries
            hit_map{t,k}(1:center_m-1,:,:) = [];  hit_map{t,k}(end-inst_M+center_m+1:end,:,:) = [];
            hit_map{t,k}(:,1:center_n-1,:) = [];  hit_map{t,k}(:,end-inst_N+center_n+1:end,:) = [];
        end
        
    end % end for t = 1 : size(pool,1)
    
end


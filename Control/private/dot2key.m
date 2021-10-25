% 根据"sta"的当前状态, 以及自识别单元的形状"pool", 生成ID矩阵.
% Build the ID matrix based on "sta" and the shape of entries.

function key = dot2key(sta,pool)
    
    [M,N] = size(sta);
    key = cell(size(pool,1),4);
    
    for t = 1 : size(pool,1)
        for shape_m = 1 : size(pool{t,5},1)
            for shape_n = 1 : size(pool{t,5},2)
                cur_shape = pool{t,5}(shape_m,shape_n);
                [inst_M,inst_N,~] = size(pool{t,cur_shape});
                key{t,cur_shape} = nan(M-inst_M+1,N-inst_N+1);
                for im = 1 : M-inst_M+1
                    for in = 1 : N-inst_N+1
                        patch = sta(im:im+inst_M-1,in:in+inst_N-1);
                        if any(any(isnan(patch)))
                            continue;
                        end
                        patch_rotated = rot90(patch,1-cur_shape);
                        patch_seq = patch_rotated(:)';
                        shape_seq = reshape(~isnan(pool{t,1}(:,:,1)),[1,inst_M*inst_N])';
                        patch_seq(~shape_seq) = [];
                        key{t,cur_shape}(im,in) = bin2dec(char(patch_seq+48));
                    end
                end
            end
        end
    end % for t = 1 : size(shape,1)
    
end


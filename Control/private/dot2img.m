% ����"sta"�ĵ�ǰ״̬, �Լ�������������ر߳�, ���ƺ�����ͼ��"img".
% Draw the image of HydraMarker in "img" based on "sta" and the pixel
% length of squares "L".

% - flag
%   0-���̸�; 1-����
%   two styles are supported: 0-chessboard; 1-grid.
function img = dot2img(sta,L,flag)
    
    [M,N] = size(sta);
    M = M+2;
    N = N+2;
    %% ����, ���̸������
    % draw the background
    if flag == 0
        img = checkerboard(L,M,N);
        img = img(1:L*M,1:L*N);
    else
        img = ones(L*M,L*N);
        line_matrix = [(0:L:L*N)',zeros(N+1,1),(0:L:L*N)',L*M*ones(N+1,1);
                       (0:L:L*N)'+1,zeros(N+1,1),(0:L:L*N)'+1,L*M*ones(N+1,1);
                       zeros(M+1,1),(0:L:L*M)',L*N*ones(M+1,1),(0:L:L*M)';
                       zeros(M+1,1),(0:L:L*M)'+1,L*N*ones(M+1,1),(0:L:L*M)'+1;];
        img = insertShape(img,'Line',line_matrix,'LineWidth',ceil(L/20),'Color','black');
    end
    
    %% ǰ��, ����"sta"���
    % draw the foreground
    circle_x = (1+L)/2:L:size(img,2);
    circle_x = repmat(circle_x,[M,1]);
    circle_x = circle_x(:);
    circle_y = (1+L)/2:L:size(img,1);
    circle_y = repmat(circle_y,[1,N]);
    circle_y = circle_y(:);
    
    circle_matrix = [circle_x,circle_y,repmat(ceil(L/6),[size(circle_x,1),1])];
    sta_pad = padarray(sta,[1,1],0,'both');
    img = insertShape(img,'FilledCircle',circle_matrix(sta_pad(:)==1,:),'Color','black','Opacity',1);
    
    % ��������̸񱳾�, ����Ҫ�ںڵ׷��鴦��׵�.
    % if the background style is chessboard, some dots should be white.
    if flag == 0
        b_square = ~xor(rem(1:N,2),rem((1:M)',2));
        img = insertShape(img,'FilledCircle',circle_matrix(sta_pad(:)==1 & b_square(:),:),'Color','white','Opacity',1);
    end
    
    %% ����һ��ı߿�
    % reduce the width of pad
    halfL = ceil(L/2);
    img = img(halfL:end-halfL+1,halfL:end-halfL+1,:);
    
end


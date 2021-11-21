%sub-pixel extraction on curved surface
function ptList = ptCurvedSurface(img, ptList, array)
    figure;
    imshow(img);
    neighbors = zeros(size(ptList,1),2);
    sample_width = 20;
    sample_freq = 11;
    NVF_x=[];          %normalVectorField
    NVF_y=[];
    global coeff_1
    global coeff_2
    for it = 1 : size(array,1)
        samplePoints_dirx_x=zeros(size(array{it},1),size(array{it},2)*sample_freq);
        samplePoints_diry_x=zeros(size(array{it},1)*sample_freq,size(array{it},2));
        samplePoints_dirx_y=zeros(size(array{it},1),size(array{it},2)*sample_freq);
        samplePoints_diry_y=zeros(size(array{it},1)*sample_freq,size(array{it},2));
        count_x = zeros(size(array{it},1),1);
        count_y = zeros(size(array{it},2),1);
        boundary_x=zeros(size(array{it},2),20);
        for ix = 1 : size(array{it},1)
            for iy = 1 : size(array{it},2)
                if isnan(array{it}(ix,iy))
                    continue;
                end
                nowPoint = array{it}(ix,iy);
                if ix < size(array{it},1) && ~isnan(array{it}(ix+1,iy))
                    neighbors(nowPoint,1) = array{it}(ix+1,iy);
                end
                if iy < size(array{it},2) && ~isnan(array{it}(ix,iy+1))
                    neighbors(nowPoint,2) = array{it}(ix,iy+1);
                end
                juncCur1_x = [];
                juncCur1_y = [];
                juncCur2_x = [];
                juncCur2_y = [];
                for j = 1 : 2
                    if neighbors(nowPoint,j) == 0 
                        continue;
                    end
                    bias_x=0;
                    bias_y=0;
                    x_divide = round(linspace(ptList(nowPoint,1),ptList(neighbors(nowPoint,j),1),sample_freq));
                    y_divide = round(linspace(ptList(nowPoint,2),ptList(neighbors(nowPoint,j),2),sample_freq));
                    for k = 2 : size(x_divide,2) - 1
                        sub_x = sample_width;
                        sub_y = sample_width;
                        if abs(ptList(nowPoint,1)-ptList(neighbors(nowPoint,j),1)) > abs(ptList(nowPoint,2)-ptList(neighbors(nowPoint,j),2))
                            [sub_y] = sigmoidFit(img(x_divide(k),y_divide(k)-sample_width+bias_y:y_divide(k)+sample_width+bias_y));
                            if sub_y <= 0 || sub_y >= sample_width * 2
                                continue;
                            end
                        else
                            [sub_x] = sigmoidFit(img(x_divide(k)-sample_width+bias_x:x_divide(k)+sample_width+bias_x,y_divide(k))');
                            if sub_x == 0 || sub_x >= sample_width * 2
                                continue;
                            end
                        end
                        if j==1
                            juncCur1_x = [juncCur1_x; x_divide(k) - sample_width + sub_x + bias_x - 1];
                            juncCur1_y = [juncCur1_y; y_divide(k) - sample_width + sub_y + bias_y - 1];
                        else
                            juncCur2_x = [juncCur2_x; x_divide(k) - sample_width + sub_x + bias_x - 1];
                            juncCur2_y = [juncCur2_y; y_divide(k) - sample_width + sub_y + bias_y - 1];
                        end
                        bias_x=bias_x+round(sub_x)-sample_width;
                        bias_y=bias_y+round(sub_y)-sample_width;
                    end   
                end
                
                if ~isempty(juncCur2_x)
                    samplePoints_dirx_x(ix, count_x(ix)+1:count_x(ix)+size(juncCur2_x,1)+1)=[ptList(nowPoint,1);juncCur2_x];
                    samplePoints_dirx_y(ix, count_x(ix)+1:count_x(ix)+size(juncCur2_x,1)+1)=[ptList(nowPoint,2);juncCur2_y];
                    count_x(ix)=count_x(ix)+size(juncCur2_x,1)+1;
                end
                if ~isempty(juncCur1_x)
                    samplePoints_diry_x(count_y(iy)+1:count_y(iy)+size(juncCur1_x,1)+1, iy)=[ptList(nowPoint,1);juncCur1_x];
                    samplePoints_diry_y(count_y(iy)+1:count_y(iy)+size(juncCur1_x,1)+1, iy)=[ptList(nowPoint,2);juncCur1_y];
                    count_y(iy)=count_y(iy)+size(juncCur1_x,1)+1;
                end
%                 hold on
%                 scatter(samplePoints_diry_y,samplePoints_diry_x,10,'g','filled');
%                 scatter(samplePoints_dirx_y,samplePoints_dirx_x,10,'b','filled');
                
                %quadratic curve fitting
%                 if exist('juncCur1_x','var') && size(juncCur1_x,1) > 4 && size(juncCur2_x,1) > 4
%                     coeff_1=fit_ellipse(juncCur1_y, juncCur1_x);
%                     fimplicit(@mfun1,[min(juncCur1_y)-100 max(juncCur1_y)+100 min(juncCur1_x)-100 max(juncCur1_x)+100],'-r','LineWidth',1.5);
%                     
%                     coeff_2=fit_ellipse(juncCur2_y, juncCur2_x);
%                     fimplicit(@mfun2,[min(juncCur2_y)-100 max(juncCur2_y)+100 min(juncCur2_x)-100 max(juncCur2_x)+100],'-y','LineWidth',1.5);
                    
                    %solve the intersection points of the two quadratic curves
%                     syms x y
%                     eqns = [coeff_1(1)*x^2+coeff_1(2)*x*y+coeff_1(3)*y^2+coeff_1(4)*x+coeff_1(5)*y+1==0, coeff_2(1)*x^2+coeff_2(2)*x*y+coeff_2(3)*y^2+coeff_2(4)*x+coeff_2(5)*y+1==0];
%                     [solx, soly] = solve(eqns, [x y]);
%                     solx = real(double(solx));
%                     soly = real(double(soly));
%                     err=(solx-ptList(nowPoint,2)).^2+(soly-ptList(nowPoint,1)).^2;
%                     [~,pos]=min(err);
%                     x_sub=solx(pos);
%                     y_sub=soly(pos);
%                     ptList_refined(nowPoint,:) = [y_sub(1) x_sub(1)];
%                     scatter(x_sub(1), y_sub(1), 30,'r','filled','o','LineWidth',1);
%                 end
            end
        end
        for ib=1:size(samplePoints_diry_x,2)
            for ia=3:size(samplePoints_diry_x,1)-2
                if samplePoints_diry_x(ia+2,ib)==0 || samplePoints_diry_x(ia-2,ib)==0 
                    continue;
                end
                slope=(samplePoints_diry_x(ia+2,ib)-samplePoints_diry_x(ia-2,ib))/(samplePoints_diry_y(ia+2,ib)-samplePoints_diry_y(ia-2,ib));
                theta=atan(slope);
                len=((samplePoints_diry_y(ia,ib)-samplePoints_diry_y(ia-2,ib))/(samplePoints_diry_y(ia+2,ib)-samplePoints_diry_y(ia-2,ib)))*(samplePoints_diry_x(ia+2,ib)-samplePoints_diry_x(ia-2,ib))+samplePoints_diry_x(ia-2,ib)-samplePoints_diry_x(ia,ib);
                u=len*sin(theta)*cos(theta)*50;
                v=len*cos(theta)*cos(theta)*50;
%                 quiver(samplePoints_diry_y(ia,ib),samplePoints_diry_x(ia,ib),-u,v,'-r','LineWidth',1.5);
                NVF_x(ia,ib).dir=atan2(v,-u)/pi*180;
                NVF_x(ia,ib).norm=sqrt(u^2+v^2);
            end
        end
        for ia=1:size(samplePoints_dirx_x,1)
            for ib=3:size(samplePoints_dirx_x,2)-2
                if samplePoints_dirx_x(ia,ib+2)==0 || samplePoints_dirx_x(ia,ib-2)==0
                    continue;
                end
                slope=(samplePoints_dirx_x(ia,ib+2)-samplePoints_dirx_x(ia,ib-2))/(samplePoints_dirx_y(ia,ib+2)-samplePoints_dirx_y(ia,ib-2));
                theta=atan(slope);
                len=((samplePoints_dirx_y(ia,ib)-samplePoints_dirx_y(ia,ib-2))/(samplePoints_dirx_y(ia,ib+2)-samplePoints_dirx_y(ia,ib-2)))*(samplePoints_dirx_x(ia,ib+2)-samplePoints_dirx_x(ia,ib-2))+samplePoints_dirx_x(ia,ib-2)-samplePoints_dirx_x(ia,ib);
                u=len*sin(theta)*cos(theta)*50;
                v=len*cos(theta)*cos(theta)*50;
%                 quiver(samplePoints_dirx_y(ia,ib),samplePoints_dirx_x(ia,ib),-u,v,'-y','LineWidth',1.5);
                NVF_y(ia,ib).dir=atan2(v,-u)/pi*180;
                NVF_y(ia,ib).norm=sqrt(u^2+v^2);
            end
        end
        for ib = 1:size(NVF_x,2)
            boundary_x(1,ib)=1;
            cnt=2;
            for ia = 6:size(NVF_x,1)-5
                if ~isempty(NVF_x(ia,ib).dir) && NVF_x(ia,ib).norm < 35 && (abs(NVF_x(ia,ib).dir-NVF_x(ia+1,ib).dir)>150 || abs(NVF_x(ia,ib).dir-NVF_x(ia-1,ib).dir)>150) && judgeBoundary(NVF_x(ia-5:ia+5,ib),5)
                    boundary_x(cnt,ib)=ia;
                    cnt=cnt+1;
                end
            end
            boundary_x(cnt,ib)=size(NVF_x,1)+2;
            hold on
            scatter(samplePoints_diry_y(boundary_x(boundary_x(:,ib)~=0,ib),ib),samplePoints_diry_x(boundary_x(boundary_x(:,ib)~=0,ib),ib),80,'filled','g');
        end    
    end
end

function [subpixel, peak, width] = sigmoidFit(Input)
    fun = @(x,xdata)x(1)./(1+exp(-1*(xdata-x(2))/x(3)))+x(4);
    x0=[0.8 size(Input,2)/2 2 0.2];
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt','Display','none');
    lb = [];    
    ub = [];
    [coe,~,~,exitFlag,~] = lsqcurvefit(fun, x0, 1:size(Input,2), Input, lb, ub, options);
    if exitFlag == 0
        subpixel = 0;
    else
        subpixel = coe(2);
        peak = coe(1);
        width = coe(3);
    end
    
%     figure(3)
%     hold off
%     plot(1:size(Input,2),Input,'*');
%     hold on
%     y=fun(coe,1:0.01:size(Input,2));
%     plot(1:0.01:size(Input,2),y,'LineWidth',2);
end

function z=mfun2(x,y)
    global coeff_2
    z = coeff_2(1).*x.^2 + coeff_2(2).*x.*y + coeff_2(3).*y.^2 + coeff_2(4).*x + coeff_2(5).*y + coeff_2(6);
end

function z=mfun1(x,y)
    global coeff_1
    z = coeff_1(1).*x.^2 + coeff_1(2).*x.*y + coeff_1(3).*y.^2 + coeff_1(4).*x + coeff_1(5).*y + coeff_1(6);
end

function a = fit_ellipse(x, y)
    D1 = [x.^2, x.*y, y.^2]; % quadratic part of the design matrix
    D2 = [x, y, ones(size(x))]; % linear part of the design matrix
    S1 = D1'* D1; % quadratic part of the scatter matrix
    S2 = D1'* D2; % combined part of the scatter matrix
    S3 = D2'* D2; % linear part of the scatter matrix
    T = - inv(S3) * S2'; % for getting a2 from a1

    M = S1 + S2 * T; % reduced scatter matrix
    M = [M(3, :) ./ 2; - M(2, :); M(1, :) ./ 2]; % premultiply by inv(C1)
    [evec, eval] = eig(M); % solve eigensystem

    cond = 4 * evec(1, :) .* evec(3, :) - evec(2, :).^2; % evaluate a’Ca
    a1 = evec(:, cond > 0); % eigenvector for min. pos. eigenvalue
    a = [a1; T * a1]; % ellipse coefficients
    if ~isreal(a)
        a=zeros(6,1);
        coeff_line = polyfit(x, y, 1);
        a(4)=-coeff_line(1);
        a(5)=1;
        a(6)=-coeff_line(2);
    end
end

function res=judgeBoundary(input,span)
    count=0;
    for i = 1:span
        if ~isempty(input(i).dir) && ~isempty(input(2*span+2-i).dir) && abs(abs(input(i).dir-input(2*span+2-i).dir)-180) < 30 && abs(abs(input(i).dir)+abs(input(2*span+2-i).dir)-180) < 30
            count = count+1;
        end
    end
    if count > 2
        res=1;
    else
        res=0;
    end
end
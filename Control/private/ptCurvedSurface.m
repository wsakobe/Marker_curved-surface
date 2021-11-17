%sub-pixel extraction on curved surface
function ptList = ptCurvedSurface(img, ptList, array)
    figure;
    neighbors = zeros(size(ptList,1),4);
    sample_width = 20;
    global coeff_1
    global coeff_2
    for it = 1 : size(array,1)
        for ix = 1 : size(array{it},1)
            for iy = 1 : size(array{it},2)
                if isnan(array{it}(ix,iy))
                    continue;
                end
                nowPoint = array{it}(ix,iy);
                if ix > 1 && ~isnan(array{it}(ix-1,iy))
                    neighbors(nowPoint,1) = array{it}(ix-1,iy);
                end
                if ix < size(array{it},1) && ~isnan(array{it}(ix+1,iy))
                    neighbors(nowPoint,2) = array{it}(ix+1,iy);
                end
                if iy > 1 && ~isnan(array{it}(ix,iy-1))
                    neighbors(nowPoint,3) = array{it}(ix,iy-1);
                end
                if iy < size(array{it},2) && ~isnan(array{it}(ix,iy+1))
                    neighbors(nowPoint,4) = array{it}(ix,iy+1);
                end
                juncCur1_x = [];
                juncCur1_y = [];
                juncCur2_x = [];
                juncCur2_y = [];
                for j = 1 : 4
                    if neighbors(nowPoint,j) == 0 
                        continue;
                    end
                    bias_x=0;
                    bias_y=0;
                    x_divide = round(linspace(ptList(neighbors(nowPoint,j),1),ptList(nowPoint,1),11));
                    y_divide = round(linspace(ptList(neighbors(nowPoint,j),2),ptList(nowPoint,2),11));
                    for k = 2 : size(x_divide,2) - 1
                        sub_x = sample_width;
                        sub_y = sample_width;
                        if abs(ptList(nowPoint,1)-ptList(neighbors(nowPoint,j),1)) > abs(ptList(nowPoint,2)-ptList(neighbors(nowPoint,j),2))
                            [sub_y,score] = sigmoidFit(img(x_divide(k),y_divide(k)-sample_width+bias_y:y_divide(k)+sample_width+bias_y));
                            if sub_y > sample_width * 2
                                continue;
                            end
                        else
                            [sub_x,score] = sigmoidFit(img(x_divide(k)-sample_width+bias_x:x_divide(k)+sample_width+bias_x,y_divide(k))');
                            if sub_x > sample_width * 2
                                continue;
                            end
                        end
                        if j<3
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
                hold off
                imshow(img);
                hold on
                scatter(juncCur2_y,juncCur2_x,10,'g','filled');
                scatter(juncCur1_y,juncCur1_x,10,'b','filled');
                
                %quadratic curve fitting
                if exist('juncCur1_x','var') && size(juncCur1_x,1) > 4 && size(juncCur2_x,1) > 4
                    coeff_1=fit_ellipse(juncCur1_y, juncCur1_x);
                    fimplicit(@mfun1,[min(juncCur1_y)-10 max(juncCur1_y)+10 min(juncCur1_x)-10 max(juncCur1_x)+10],'-r','LineWidth',1.5);
                    
                    coeff_2=fit_ellipse(juncCur2_y, juncCur2_x);
                    fimplicit(@mfun2,[min(juncCur2_y)-10 max(juncCur2_y)+10 min(juncCur2_x)-10 max(juncCur2_x)+10],'-y','LineWidth',1.5);

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
                end
            end
        end
    end
end

function [subpixel, peak, width] = sigmoidFit(Input)
    min_peak = 0.2; %折扣因子
    max_width = 5;
    fun = @(x,xdata)x(1)./(1+exp(-1*(xdata-x(2))/x(3)))+x(4);
    x0=[0.8 size(Input,2)/2 2 0.2];
    coe = lsqcurvefit(fun, x0, 1:size(Input,2), Input);
    subpixel = coe(2);
    peak = coe(1);
    width = coe(3);
    
    figure(3)
    hold off
    plot(1:size(Input,2),Input,'*');
    hold on
    y=fun(coe,1:0.01:size(Input,2));
    plot(1:0.01:size(Input,2),y,'LineWidth',2);
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
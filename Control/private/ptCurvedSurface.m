%sub-pixel extraction on curved surface
function ptList = ptCurvedSurface(img, ptList, array)
    figure;
    neighbors = zeros(size(ptList,1),4);
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
                    if ptList(neighbors(nowPoint,j),1) < ptList(nowPoint,1)
                        x_dir = 1;
                    else
                        x_dir = -1;
                    end
                    if ptList(neighbors(nowPoint,j),2) < ptList(nowPoint,2)
                        y_dir = 1;
                    else
                        y_dir = -1;
                    end
                    x_divide = round(linspace(ptList(neighbors(nowPoint,j),1),ptList(nowPoint,1),7));
                    y_divide = round(linspace(ptList(neighbors(nowPoint,j),2),ptList(nowPoint,2),7));
                    for k = 2 : size(x_divide,2) - 2
                        stdPixel_x = img(x_divide(k),y_divide(k));
                        stdPixel_y = img(x_divide(k),y_divide(k));
                        changeCnt_x = 0;
                        changeCnt_y = 0;
                        cnt = 1;
                        while changeCnt_x < 3 && changeCnt_y < 3 && cnt < 50
                            nowPixel_x = img(x_divide(k) + x_dir * cnt, y_divide(k));
                            nowPixel_y = img(x_divide(k), y_divide(k) + y_dir * cnt);
                            if abs(nowPixel_x - stdPixel_x) > 0.15
                                changeCnt_x = changeCnt_x + 1;
                            else
                                changeCnt_x = 0;
                            end
                            if abs(nowPixel_y - stdPixel_y) > 0.15
                                changeCnt_y = changeCnt_y + 1;
                            else
                                changeCnt_y = 0;
                            end
                            cnt = cnt + 1;  
                        end
                        if changeCnt_x == 3 
                            if j < 3
                                juncCur1_x = [juncCur1_x; x_divide(k) + x_dir * (cnt - 5) - 10 + sigmoidFit(img(x_divide(k)+x_dir*(cnt-5)-10:x_divide(k)+x_dir*(cnt-5)+10,y_divide(k))')];
                                juncCur1_y = [juncCur1_y; y_divide(k)];
                            else
                                juncCur2_x = [juncCur2_x; x_divide(k) + x_dir * (cnt - 5) - 10 + sigmoidFit(img(x_divide(k)+x_dir*(cnt-5)-10:x_divide(k)+x_dir*(cnt-5)+10,y_divide(k))')];
                                juncCur2_y = [juncCur2_y; y_divide(k)];
                            end
                        else
                            if j < 3
                                juncCur1_x = [juncCur1_x; x_divide(k)];
                                juncCur1_y = [juncCur1_y; y_divide(k) + y_dir * (cnt - 5) - 10 + sigmoidFit(img(x_divide(k),y_divide(k)+y_dir*(cnt-5)-10:y_divide(k)+y_dir*(cnt-5)+10))];
                            else
                                juncCur2_x = [juncCur2_x; x_divide(k)];
                                juncCur2_y = [juncCur2_y; y_divide(k) + y_dir * (cnt - 5) - 10 + sigmoidFit(img(x_divide(k),y_divide(k)+y_dir*(cnt-5)-10:y_divide(k)+y_dir*(cnt-5)+10))];
                            end
                        end
                    end   
                end
                hold off
                imshow(img);
                hold on
                scatter(juncCur2_y,juncCur2_x,10,'g','filled');
                scatter(juncCur1_y,juncCur1_x,10,'b','filled');
                
                %quadratic curve fitting
                if exist('juncCur1_x','var') && size(juncCur1_x,1) > 4 && size(juncCur2_x,1) > 4
                    A_1=[juncCur1_y.^2 juncCur1_x.*juncCur1_y juncCur1_x.^2 juncCur1_y juncCur1_x];
                    A_1=reshape(A_1,size(juncCur1_x,1),5);
                    B_1=-1*ones(size(juncCur1_x,1),1);
                    coeff_1=A_1\B_1;
                    fimplicit(@mfun1,[min(juncCur1_y) max(juncCur1_y) min(juncCur1_x) max(juncCur1_x)],'-r','LineWidth',1.5);

                    A_2=[juncCur2_y.^2 juncCur2_x.*juncCur2_y juncCur2_x.^2 juncCur2_y juncCur2_x];
                    A_2=reshape(A_2,size(juncCur2_x,1),5);
                    B_2=-1*ones(size(juncCur2_x,1),1);
                    coeff_2=A_2\B_2;
                    fimplicit(@mfun2,[min(juncCur2_y) max(juncCur2_y) min(juncCur2_x) max(juncCur2_x)],'-y','LineWidth',1.5);

                    %solve the intersection points of the two quadratic curves
                    syms x y
                    eqns = [coeff_1(1)*x^2+coeff_1(2)*x*y+coeff_1(3)*y^2+coeff_1(4)*x+coeff_1(5)*y+1==0, coeff_2(1)*x^2+coeff_2(2)*x*y+coeff_2(3)*y^2+coeff_2(4)*x+coeff_2(5)*y+1==0];
                    [solx, soly] = solve(eqns, [x y]);
                    solx = real(double(solx));
                    soly = real(double(soly));
                    x_sub = solx(abs(solx-ptList(nowPoint,2))==min(abs(solx-ptList(nowPoint,2))));
                    y_sub = soly(abs(soly-ptList(nowPoint,1))==min(abs(soly-ptList(nowPoint,1))));
                    ptList(nowPoint,:) = [y_sub(1) x_sub(1)];
                    scatter(x_sub(1), y_sub(1), 30,'r','filled','o','LineWidth',1);
                end
            end
        end
    end
end

function subpixel = sigmoidFit(Input)
    fun = @(x,xdata)x(1)./(1+exp(-1*(xdata-x(2))/x(3)))+x(4);
    x0=[0.8 size(Input,2)/2 2 0.2];
    coe = lsqcurvefit(fun, x0, 1:size(Input,2), Input);
    subpixel = coe(2);
%     
%     figure
%     plot(1:size(Input,2),Input);
%     hold on
%     y=fun(coe,1:0.01:size(Input,2));
%     plot(1:0.01:size(Input,2),y);
end

function z=mfun2(x,y)
    global coeff_2
    z = coeff_2(1).*x.^2 + coeff_2(2).*x.*y + coeff_2(3).*y.^2 + coeff_2(4).*x + coeff_2(5).*y + 1;
end

function z=mfun1(x,y)
    global coeff_1
    z = coeff_1(1).*x.^2 + coeff_1(2).*x.*y + coeff_1(3).*y.^2 + coeff_1(4).*x + coeff_1(5).*y + 1;
end
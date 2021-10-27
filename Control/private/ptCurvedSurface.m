%曲面交叉点亚像素�?�?
%sub-pixel extraction on curved surface
function ptList = ptCurvedSurface(img, ptList, edge)
    figure;
    imshow(img);
    hold on
    for it = 1 : size(ptList,1)
        find(edge(:,1)==ptList(it,1));
        surroundPoint(it,:) = ptList(edge(:,1)==ptList(it,1));
        %quadratic curve fitting
        if exist('juncCur1_x','var') && size(juncCur1_x,2) > 4 && size(juncCur2_x,2) > 4
            xmin_1=min(juncCur1_x);
            xmax_1=max(juncCur1_x);
            A_1=[juncCur1_y.^2 juncCur1_x.*juncCur1_y juncCur1_x.^2 juncCur1_y juncCur1_x];
            A_1=reshape(A_1,size(juncCur1_x,2),5);
%             if (rank(A_1)<5)
%                 break;
%             end
            B_1=-1*ones(size(juncCur1_x,2),1);
            coeff_1=A_1\B_1;

            xmin_2=min(juncCur2_y);
            xmax_2=max(juncCur2_y);
            A_2=[juncCur2_y.^2 juncCur2_x.*juncCur2_y juncCur2_x.^2 juncCur2_y juncCur2_x];
            A_2=reshape(A_2,size(juncCur1_x,2),5);
%             if (rank(A_2)<5)
%                 break;
%             end
            B_2=-1*ones(size(juncCur1_x,2),1);
            coeff_2=A_2\B_2;

            syms x y;
            h1=ezplot(coeff_1(1)*y^2+coeff_1(2)*x*y+coeff_1(3)*x^2+coeff_1(4)*y+coeff_1(5)*x+1==0,[10 50]);
            h2=ezplot(coeff_2(1)*y^2+coeff_2(2)*x*y+coeff_2(3)*x^2+coeff_2(4)*y+coeff_2(5)*x+1==0,[10 50]);
            set(h1,'color',[0,1,0],'LineWidth',2)
            set(h2,'color',[1,0,0],'LineWidth',2)
            scatter(juncCur1_x, juncCur1_y,30,'g','filled','o','LineWidth',1); 
            scatter(juncCur2_x, juncCur2_y,30,'b','filled','o','LineWidth',1); 
                        
            %solve the intersection points of the two quadratic curves
            s=solve(coeff_1(1)*y^2+coeff_1(2)*x*y+coeff_1(3)*x^2+coeff_1(4)*y+coeff_1(5)*x+1==0,coeff_2(1)*y^2+coeff_2(2)*x*y+coeff_2(3)*x^2+coeff_2(4)*y+coeff_2(5)*x+1==0,x,y);
            X=double(s.x);Y=double(s.y);
            [x_sub,y_sub] = min(abs(X-im)+abs(Y-in));
            scatter(x_sub, y_sub, 50,'r','filled','o','LineWidth',1);
        end
    end
end
%曲面交叉点亚像素检测
%sub-pixel extraction on curved surface
function ptList = ptCurvedSurface(img, ptList, edge)
    figure;
    imshow(img);
    hold on
    for it = 1 : size(list,1)
        cnt1=1;
        cnt2=1;
        for maskR = 2:8
            clear frameOff_x frameOff_y frameSgn frameVal;
            frameOff_x(1:maskR*2) = -maskR:maskR-1;
            frameOff_x(maskR*2+1:4*maskR) = maskR*ones(1,maskR*2);
            frameOff_x(maskR*4+1:6*maskR) = maskR:-1:-maskR+1;
            frameOff_x(maskR*6+1:8*maskR) = -maskR*ones(1,maskR*2);
            frameOff_y(1:maskR*2) = -maskR*ones(1,maskR*2);
            frameOff_y(maskR*2+1:4*maskR) = -maskR:maskR-1;
            frameOff_y(maskR*4+1:6*maskR) = maskR*ones(1,maskR*2);
            frameOff_y(maskR*6+1:8*maskR) = maskR:-1:-maskR+1;

            im = list(it,2);
            in = list(it,1);
            for i = 1:size(frameOff_x,2)
                frameVal(i) = img(frameOff_y(i)+in, frameOff_x(i)+im);
            end
            frameMean = sum(frameVal)/size(frameOff_x,2);
            for j = 1:size(frameOff_x,2)
                if img(frameOff_y(j)+in, frameOff_x(j)+im) > frameMean
                    frameSgn(j) = 1;
                else
                    frameSgn(j) = 0;
                end
            end
            jmpCnt=zeros(2,1);
            for j = 1:size(frameSgn,2)
                if j==1 
                    ib=size(frameSgn,2);
                else
                    ib=j-1;
                end
                if (frameSgn(ib)~=frameSgn(j))
                    if jmpCnt(frameSgn(j)+1)<2
                        jmpIdx(frameSgn(j)+1,jmpCnt(frameSgn(j)+1)+1) = j;
                    end
                    jmpCnt(frameSgn(j)+1)=jmpCnt(frameSgn(j)+1)+1;
                end
            end

            if (jmpCnt(1)==2) && (jmpCnt(2)==2)
                juncCur1_x(cnt1)=frameOff_x(1,jmpIdx(1,1))+im;
                juncCur1_x(cnt1+1)=frameOff_x(1,jmpIdx(1,2))+im;
                juncCur1_y(cnt1)=frameOff_y(1,jmpIdx(1,1))+in;
                juncCur1_y(cnt1+1)=frameOff_y(1,jmpIdx(1,2))+in;
                cnt1=cnt1+2;

                juncCur2_x(cnt2)=frameOff_x(1,jmpIdx(2,1))+im;
                juncCur2_x(cnt2+1)=frameOff_x(1,jmpIdx(2,2))+im;
                juncCur2_y(cnt2)=frameOff_y(1,jmpIdx(2,1))+in;
                juncCur2_y(cnt2+1)=frameOff_y(1,jmpIdx(2,2))+in;
                cnt2=cnt2+2;
            else
                break;
            end
        end
        
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
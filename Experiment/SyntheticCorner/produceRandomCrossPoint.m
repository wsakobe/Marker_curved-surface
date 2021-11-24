% produce a random cross point with projective transform, blur and noise
% output size 64x64
% input:
%       blurDegree , the sigma used in function imgaussfilt
%       noiseDegree, the sigma used in function imnoise
% output:
%       img, random cross point
%       angle, the angle between optical axis and cross point surface
function [img] = produceRandomCrossPoint(blurDegree,noiseDegree,cnt)
    if ~exist('blurDegree','var')
        blurDegree = 1;
    end
    if ~exist('noiseDegree','var')
        noiseDegree = 0.001;
    end
        
    halfM = 320;
    halfN = 320;
    M = 2*halfM;
    N = 2*halfN;
    cnt=1
    
    img = zeros(M, N);
    coeff_x(1) = rand/2000;
    coeff_y(1) = rand/2000;
    coeff_x(2) = rand/100;
    coeff_y(2) = rand/100;
    coeff_x(3) = halfM - halfM^2 * coeff_x(1) - halfM * coeff_x(2);
    coeff_y(3) = halfN - halfN^2 * coeff_y(1) - halfN * coeff_y(2);
    for i = 64:576
        for j = 64:576
            if ((j < coeff_x(1)*i^2 + coeff_x(2)*i + coeff_x(3)) && (i > coeff_y(1)*j^2 + coeff_y(2)*j + coeff_y(3))) || ((j > coeff_x(1)*i^2 + coeff_x(2)*i + coeff_x(3)) && (i < coeff_y(1)*j^2 + coeff_y(2)*j + coeff_y(3))) 
                img(j,i) = 1;
            end
        end
    end
    img(577:640,577:640)=1;
    img(1:64,1:64)=1;

    for i = 64:576
        for j = 1:64
            if i < coeff_y(1)*j^2 + coeff_y(2)*j + coeff_y(3)
                img(j,i)=1;
            end
        end
    end
    for i = 64:576
        for j = 1:64
            if i < coeff_x(1)*j^2 + coeff_x(2)*j + coeff_x(3)
                img(i,j)=1;
            end
        end
    end
    for i = 64:576
        for j = 577:640
            if i > coeff_y(1)*j^2 + coeff_y(2)*j + coeff_y(3)
                img(j,i)=1;
            end
        end
    end
    for i = 64:576
        for j = 577:640
            if i > coeff_x(1)*j^2 + coeff_x(2)*j + coeff_x(3)
                img(i,j)=1;
            end
        end
    end
    imshow(img);
    
%     img(1:halfM,1:halfN) = 1;
%     img(halfM+1:M,halfN+1:N) = 1;

    % random viewpoint
%     angle = pi/3*rand;
%     R = angle2dcm(2*pi*rand,angle,0);
%     R(3,:) = [0,0,1];
%     img = imwarp(img,imref2d([M,N],[-1, 1],[-1, 1]),projective2d(R),...
%         'OutputView',imref2d([M,N],[-2, 2],[-2, 2]),'FillValues',0.5);
    
    % blur and noise
    img = imresize(img,[160 160]);
    if (blurDegree)
        img = imgaussfilt(img,blurDegree);
    end
    img = imnoise(img,'gaussian',0,noiseDegree);
    imshow(img);
    filename=['.\RandomCornerNoise\Noise',num2str(noiseDegree),'\',num2str(cnt),'.bmp'];
    imwrite(img,filename);
end


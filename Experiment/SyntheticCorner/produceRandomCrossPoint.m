% produce a random cross point with projective transform, blur and noise
% output size 64x64
% input:
%       blurDegree , the sigma used in function imgaussfilt
%       noiseDegree, the sigma used in function imnoise
% output:
%       img, random cross point
%       angle, the angle between optical axis and cross point surface
function img = produceRandomCrossPoint(blurDegree,noiseDegree,cnt)
    if ~exist('blurDegree','var')
        blurDegree = 3;
    end
    if ~exist('noiseDegree','var')
        noiseDegree = 0.003;
    end
        
    halfM = 320;
    halfN = 320;
    M = 2*halfM;
    N = 2*halfN;
    cnt=4;
    
    img = zeros(M, N);
    radius1 = 600+rand*600;
    radius2 = 600+rand*600;
    circle1 = [radius1+halfM halfM];
    circle2 = [halfN radius2+halfN];
    coeff_x(1) = 1e-4 + rand/4000
    coeff_y(1) = 1e-4 + rand/4000
    coeff_x(2) = rand/100;
    coeff_y(2) = rand/100;
    coeff_x(3) = halfM - halfM^2 * coeff_x(1) - halfM * coeff_x(2);
    coeff_y(3) = halfN - halfN^2 * coeff_y(1) - halfN * coeff_y(2);
    for i = 1:M
        for j = 1:M
            point = [i j];
            if (((norm(point - circle1) <= radius1) && (norm(point - circle2) <= radius2)) || ((norm(point - circle1) > radius1) && (norm(point - circle2) > radius2))) 
                img(j,i) = 1;
            end
        end
    end

    for i = M/10:M*9/10
        for j = M/10:M*9/10
           img(j,i) = ~img(j,i);
        end
    end

    % random viewpoint
%     angle = pi/6*rand;
%     R = angle2dcm(2*pi*rand,angle,0);
%     R(3,:) = [0,0,1];
%     img = imwarp(img,imref2d([M,N],[-1, 1],[-1, 1]),projective2d(R),...
%         'OutputView',imref2d([M,N],[-1, 1],[-1.5, 1.5]),'FillValues',0.5);
    
    % blur and noise
%     img = imresize(img,[160 160]);
    if (blurDegree)
        img = imgaussfilt(img,blurDegree);
    end
    img = imnoise(img,'gaussian',0,0.07);
    img = img / 2;
    imshow(img);
    filename=['.\RandomCornerNoise\Noise',num2str(noiseDegree),'\',num2str(cnt),'.bmp'];
    imwrite(img,filename);
end


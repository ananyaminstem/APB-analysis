close all;
clear all;

[FileName,PathName,FilterIndex]=uigetfile('*.tif','Select the pre-bleach image');
fnam = strcat(PathName,FileName);
FileName=FileName(1:end-4);
BW4 = imread(fnam);
pre_bleach=BW4(:,:,1);
figure;imshow(pre_bleach,[]);

[FileName1,PathName1,FilterIndex1]=uigetfile('*.tif','Select the post-bleach image');
fnam1 = strcat(PathName1,FileName1);
FileName1=FileName1(1:end-4);
BW5 = imread(fnam1);
post_bleach=BW5(:,:,1);
figure;imshow(post_bleach,[]);

diff_img = post_bleach - pre_bleach;
figure;imshow(diff_img,[]);


T = graythresh(pre_bleach);

Im_bw = im2bw(post_bleach,0.006);

figure;imshow(Im_bw,[])

 

fh = imfill(Im_bw,'holes');

Im_L = bwlabeln(fh, 8);

S = regionprops(Im_L, 'Area');

Im_regProp = ismember(Im_L, find([S.Area] >= 2500));

figure;imshow(Im_regProp,[])

 BW1 = edge(Im_regProp,'canny');

figure;imshow(BW1,[]);

S1 = regionprops(Im_regProp, 'Centroid','MajorAxisLength','PixelIdxList','PixelList');

centroid = S1.Centroid;

% pixels =S1.PixelList;
% figure;
[C,h]=contour(Im_regProp,1);
hold on;plot(centroid(1), centroid(2),'g+');
tic
figure;imshow(post_bleach,[]);
for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(post_bleach,xi,yi,200);

    int(:,i) = c;
    if mod(i,50)==0
        hold on;plot(xi,yi,'r')
    end
end

for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(pre_bleach,xi,yi,200);

    int1(:,i) = c;

end

for i =2:size(C,2)

    xi = [C(1,i) round(centroid(1))];

    yi = [C(2,i) round(centroid(2))];

    [cx,cy,c,xi,yi] = improfile(diff_img,xi,yi,200);

    int2(:,i) = c;

end

toc

mean_int_postbleach = mean(int,2);
mean_int_prebleach = mean(int1,2);
mean_int_diffimg = mean_int_postbleach - mean_int_prebleach;
% mean_int_diffimg = mean(int2,2);


figure; plot(mean_int_prebleach);
hold on; plot(mean_int_postbleach,'r');
hold on; plot(mean_int_diffimg,'g');

legend('pre bleach','post bleach','difference');

% figure;
for i=2:size(C,2)
    [max_int(i),loc(i)] = max(int(:,i));

%   plot(int(:,1:100)); hold on; 

end
% figure; plot(loc(1:10));


 

%%

% 
% 
% 
%  
% 
% S2 = regionprops(BW1, 'PixelIdxList','PixelList');
% 
%  
% 
% [BW1,thresh,gv,gh] = edge(Im_regProp,'sobel');
% 
% edgeDir = atan2(gv, gh);
% 
% figure;imshow(BW1,[])
% 
% figure;imshow(edgeDir,[])
% 
%  
% 
% [Nx,Ny]=surfnorm(Im_regProp);
% 
% figure;imshow(Nx,[])
% ujvbl000000000000000000000000000000000000000000000000

% quiver(Nx);
% 
%  
% 
% im_f=imfuse(BW,BW1);
% 
% figure;imshow(im_f,[]);
% 
%  
% 
% Im_circ_cell=zeros(512,512);
% 
% for i = 1:size(BW,1)
% 
%     for j=1:size(BW,1)
% 
%         if (abs(i-size(BW,1)/2))^2 + (abs(j-size(BW,1)/2))^2 < (size(BW,1)/2)^2
% 
%             Im_circ_cell(i,j)=1; 
% 
%         end
% 
%     end
% 
% end
% 
% figure;imshow(Im_circ_cell,[])
% 
% Im_edge_circ = edge(Im_circ_cell,'canny');
% 
% figure;imshow(Im_edge_circ,[]);colormap(hot);title('canny');title('Edge detection Dapi')
% 
%  
% 
% im_f1=imfuse(BW,Im_edge_circ);
% 
% figure;imshow(im_f1,[]);
% 
%  
% 
%  
% 
% xi = [ S1.Centroid(1), ]
% 
% [cx,cy,c] = improfile(BW,xi,yi,n)
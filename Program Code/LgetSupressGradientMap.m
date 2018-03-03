function [px,py]=LgetSupressGradientMap(IM,SupressFactor,shown)
f = double(IM)/255;
% f0 = gaussianBlur(f,1);
f0 = wiener2(f,[10 10]);
%     disp(' Comute the traditional external force ...');
[px,py] = gradient(f0);
%% make the border gradient to be zeros
px([1:2, end-1:end],:)=0;py([1:2, end-1:end],:)=0;
px(:,[1:2, end-1:end])=0;py(:,[1:2, end-1:end])=0;

%% supress the weak gradient
G_mag=sqrt(px.^2+ py.^2);
% TGMag=mean(G_mag(:));   %max(G_mag(:));%figure(23);hist(G_mag,50);
%       allMag=G_mag(:);allMag(find(allMag>150/255))=[];
%       allMag(find(allMag<10/255))=[];
%% !!!! important thrshold for  supressing the weak gradient
TGMag=mean(G_mag(:));
ValidGM=G_mag>=TGMag*SupressFactor;
px(~ValidGM)=0;   py(~ValidGM)=0;

%% display the results
if shown
    figure(2);imshow(IM,'InitialMagnification','fit');hold on;
    quiver(px,py,5,'y');
    hold off;
end
end
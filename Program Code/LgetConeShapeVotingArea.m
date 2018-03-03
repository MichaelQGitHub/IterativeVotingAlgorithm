% cal the cone shape voting area describe in 
% X.Qi,F.Xing,D.J.Foran,andL.Yang,
% ??Robustsegmentationofoverlappingcellsinhistopathologyspecimensusingparallelseed detectionandrepulsivelevelset,??
% IEEETransactionsonBiomedicalEngineering,vol.59,no.3,pp.754?C765,2012
% x, y represent the current piont
% ux, uy represent the shift Gaussian center piont
% r_max, r_min represent the maximum and minimum range of the cone
% theta is the one side anlge allowed for the cone
% imsize is the size of image
% this function returns the pixel list that are located in the valie cone
% shape area for the voting
function [ptsIdx_ValidArea,bw_valid]=LgetConeShapeVotingArea(x,y,ux,uy,r_max,r_min,theta,imsize)

% %% test parameters
% x=100;
% y=100;
% ux=150;
% uy=150;
% r_max=500;
% r_min=300;
% theta=pi/6;
% imsize=[1000 1000];

%% get ring first
% bw=zeros(imsize);
% bw=(x-ux).^2+(y-uy).^2<r_max;

[px,py] = meshgrid(1:imsize(2),1:imsize(1));

bw1=((px-x).^2+(py-y).^2)<r_max^2;
% show(bw1);
bw2=((px-x).^2+(py-y).^2)<r_min^2;
% show(bw2);

bw_ring=bw1&~bw2;
% show(bw_ring);

%% get beam
alpha=atan((uy-y)/(ux-x));
% current Graident voting  Angle is denoted as alpha
% alpha=atan((uy-y)/(ux-x));
% cal k for two lines
k1=tan(alpha-theta);
k2=tan(alpha+theta);
b1=y-k1*x;
b2=y-k2*x;

bw3=k1*px-py+b1<0;
% if isnan(uy)
%     pause();
% end
if bw3(floor(uy),floor(ux))~=true
    bw3=k1*px-py+b1>0;
end
%show(bw3);

bw4=k2*px-py+b2<0;
if bw4(floor(uy),floor(ux))~=true
    bw4=k2*px-py+b2>0;
end
%show(bw4);


bw_beam=bw3& bw4;
% show(bw_beam);
%% get cone shape valid area
% LshowMaskCountouronIM(bw_beam,bw_ring);%hold on;
  
bw_valid=bw_beam&bw_ring;
% show(bw_valid);

% LshowBWonIM(bw_ring,bw_valid,1);hold on;
% plot(x,y,'r*');
% plot(ux,uy,'rs');
% plot([1:imsize(2)],[k1.*[1:imsize(2)]+b1],'-b');
% plot([1:imsize(2)],[k2.*[1:imsize(2)]+b2],'-b');
% hold off
%% 
ptsIdx_ValidArea=find(bw_valid(:)>0);
end
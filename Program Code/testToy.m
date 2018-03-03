clear all

IM=imread('/5DiffCircles.jpg');
IM=IM(:,:,1);

%% show gradient map
[px,py]=LgetSupressGradientMap( IM(:,:,1),0,0);

temp_edge=edge(IM(:,:,1), 'canny',.5);
% show(temp_edge,3);

bw_edge=temp_edge;
px(~bw_edge)=0;py(~bw_edge)=0;
show(bw_edge,12);hold on;
quiver(px,py,5,'y');
hold off;

%% test the Symetric Voting method

Para.VotingGap=2;
Para.rmin=1;  % determine the voting range
Para.rmax=66;
% the angle of cone-shape
% delta = pi/6;
Para.theta=pi/6;

% Gaussian Variance
Para.Sigma=4;
Para.debug=0;
Para.ConeshapeRestrict=0;
Para.N=4;
theta_min=pi/30;
Para.thetaSet=[theta_min:(Para.theta-theta_min)/(Para.N-1):Para.theta];
% Para.ObjColor='Black';
Para.ObjColor='white';

[Gx,Gy] = gradient(double(IM));

[im_Vote,Allim_Vote]= LIterativeVoting2007(bw_edge,Gx,Gy,Para);
%% show result

LshowBWonIM(bw_edge,Allim_Vote(:,:,1));
for i=1:Para.N
    LshowBWonIM(bw_edge,Allim_Vote(:,:,i));
end
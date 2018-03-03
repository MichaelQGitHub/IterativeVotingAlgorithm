function [im_Vote,Allim_Vote]= LIterativeVoting2007(bw_edge,Gx,Gy,Para)
% Input
%     -bw_i the binary edge image
%     -img the original image for observation
% Output
%     -RC_i the mask image with markers
% Program written by Cheng LU
% Shaanxi Normal University
% 2014 June 24th

% This function implement the method proposed in
% Parvin, B., Yang, Q., Han, J., Chang, H., Rydberg, B., & Barcellos-Hoff,
% M. H. (2007). Iterative voting for inference of structural saliency and
% characterization of subcellular events. IEEE Transactions on Image
% Processing : A Publication of the IEEE Signal Processing Society, 16(3),
% 615?23. Retrieved from http://www.ncbi.nlm.nih.gov/pubmed/17357723

% if ~exist('debug','var')
%     Para.debug=0;
% end
%% parameters setting
VotingGap=Para.VotingGap; % we don't need all pts for the voting, this parameter defineds the sampling gap
rmin = Para.rmin;
rmax = Para.rmax;
IterationN=Para.N;
% the angle of cone-shape
% delta = pi/6;
% theta=Para.theta;
% Gaussian Variance
Sigma = Para.Sigma;

bw_edge_d=bw_edge;
% get gradient magnitude
Gx(~bw_edge_d)=0;
Gy(~bw_edge_d)=0;
Gmag = sqrt(Gx.^2+Gy.^2);
%% change all the gradient direction to its opposit if object is darker than background
if strcmp(Para.ObjColor,'Black')
    Gx=-Gx;
    Gy=-Gy;
end
map_Acc=zeros(size(bw_edge));

%% compute the number of connected components
[m,n] = size(bw_edge_d);
c1 = bwconncomp(bw_edge_d);

for iteration=IterationN:-1:1 % control the iteration
    % accumulated map for recording the voting at each iteration, the V in the paper
    map_Acc=zeros(size(bw_edge));
    numVotingPts=1; %ind for the pts involved in voting
    theta=Para.thetaSet(iteration);
    
    for i = 1:c1.NumObjects
        fprintf('On the %d/%d edge\n',i,c1.NumObjects);
        curPtList=c1.PixelIdxList{i};
         for j=1:VotingGap:length(curPtList)
             %fprintf('j=%d\n',j);
            curPtIdx=curPtList(j);
            [r,v]=ind2sub([m,n],curPtIdx);
            % at first, the voting direction is based on gradient direction
            if Gmag(curPtIdx)==0
                continue;
            else
                if iteration==IterationN
                    cos_theta =  Gx(curPtIdx)/Gmag(curPtIdx);
                    sin_theta =  Gy(curPtIdx)/Gmag(curPtIdx);
                    % compute shift Gaussian kernel center
                    ux = v+(rmax+rmin)*(cos_theta)/2;
                    uy = r+(rmax+rmin)*(sin_theta)/2;
                    if ~(ux<1 || uy<1 || ux>n || uy>m)
                        Allux_old(numVotingPts)=ux;
                        Alluy_old(numVotingPts)=uy;
                        numVotingPts=numVotingPts+1;
                    end
                else % update the voting direction based on the maximun pixel in
                    % pre-voting area of pre-voting map
                    cos_theta =  Gx(curPtIdx)/Gmag(curPtIdx);
                    sin_theta =  Gy(curPtIdx)/Gmag(curPtIdx);
                    % compute shift Gaussian kernel center
                    ux = v+(rmax+rmin)*(cos_theta)/2;
                    uy = r+(rmax+rmin)*(sin_theta)/2;
                    if ~(ux<1 || uy<1 || ux>n || uy>m)
                        [ptsIdx_ValidArea,bw_valid]=...
                            LgetConeShapeVotingArea(v,r,Allux_old(numVotingPts),Alluy_old(numVotingPts),...
                            rmax,rmin,Para.thetaSet(iteration+1),[m n]);
                        bw_valid_old=bw_valid;
                        % find the location of maximum pixel
                        map_Acc_old=Allim_Vote(:,:,iteration+1);
                        [maxV,maxIdx]=max(map_Acc_old(bw_valid));
                        maxIdxinIM=ptsIdx_ValidArea(maxIdx);
                        [max_r,max_c]=ind2sub([m,n],maxIdxinIM);
                        if Para.debug
                            LshowBWonIM(bw_valid,map_Acc_old,1);hold on;
                            plot(max_c,max_r,'*r');
                            hold off;
                        end
                        
                        dx=max_c-v;dy=max_r-r;
                        cos_theta=dx/sqrt(dx^2+dy^2);
                        sin_theta=dy/sqrt(dx^2+dy^2);
                        % compute shift Gaussian kernel center
                   
                        ux = v+(rmax+rmin)*(cos_theta)/2;
                        uy = r+(rmax+rmin)*(sin_theta)/2;
                        if ~(ux<1 || uy<1 || ux>n || uy>m)
                            Allux_old(numVotingPts)=ux;
                            Alluy_old(numVotingPts)=uy;
                            numVotingPts=numVotingPts+1;
                        else % keep previous result
                            numVotingPts=numVotingPts+1;
                        end
                    end
                end
            end
            % check if the shift kernel center is outside the image,
            % if so we ignore this pt
            if ux<1 || uy<1 || ux>n || uy>m
                continue;
            else % voting in valid region, this region is defined as a cone shape
                
                [ptsIdx_ValidArea,bw_valid]=LgetConeShapeVotingArea(v,r,ux,uy,rmax,rmin,theta,[m n]);
                
                if Para.debug
                    if iteration~=IterationN
                        LshowTwoKindofCountouronIM(bw_valid_old,bw_valid,map_Acc_old,1);
                        hold on;
                        plot(max_c,max_r,'*r');
                        quiver(Gx,Gy,5,'y');
                        hold off;
                        %                         pause(.5);
                    else
                        LshowBWonIM(bw_valid,map_Acc,1);
                        hold on;
                        quiver(Gx,Gy,5,'y');
                        hold off;
                    end
                end
                
                if isempty(ptsIdx_ValidArea)
                    continue;
                else
                    %%% cal the Gaussian kernel
                    [px,py] = meshgrid(1:size(bw_edge,2),1:size(bw_edge,1));
                    gau = 1/(sqrt(2*pi)*Sigma).*exp(-(sum(([px(:),py(:)]-[repmat(ux,size(px(:)),1),repmat(uy,size(px(:)),1)]).^2,2))./(2*Sigma^2));
                    map_gau=reshape(gau,size(bw_edge,1),size(bw_edge,2));
                    map_gau(~bw_valid)=0;
                    map_Acc=map_Acc+Gmag(curPtIdx).*map_gau;
                    %                 map_Acc=map_Acc+map_gau;
                    if Para.debug
                        %                         show(bw_valid,1);
                        %                         show(map_Acc,3);
                        %                         show(map_gau,2);
                        %                         LshowMaskCountouronIM(bw_valid,map_gau,4);
                    end
                end
            end
        end
    end
    Allim_Vote(:,:,iteration)=map_Acc;
end
im_Vote=map_Acc;



function [Sm,Sb,Sf] = SCL(nodesPositionArray,elementNodesArray, sideWest,elementStressAtNodes)

% stress clasification line
tSCL =  abs(nodesPositionArray(sideWest(1),2)-nodesPositionArray(sideWest(end),2));
stressAtWestNode(size(sideWest,1),3) = 0; %ordenados de manera acendente en y
stressAtWestNode2(size(sideWest,1),1) = 0; %ordenados de manera acendente en y
y(size(sideWest,1),1) = 0;
for i = 1:size(sideWest,1)    
    [fil,col] = find(elementNodesArray == sideWest(i));    
    if length(fil)>1
%         ele1 = [fil(1) col(1)];
%         ele2 = [fil(2) col(2)];
        stressAtWestNode(i,1) = (elementStressAtNodes(fil(1),col(1),1)+elementStressAtNodes(fil(2),col(2),1))/2;
        stressAtWestNode(i,2) = (elementStressAtNodes(fil(1),col(1),2)+elementStressAtNodes(fil(2),col(2),2))/2;
        stressAtWestNode(i,3) = (elementStressAtNodes(fil(1),col(1),3)+elementStressAtNodes(fil(2),col(2),3))/2;
    else
        stressAtWestNode(i,1) = elementStressAtNodes(fil,col,1);
        stressAtWestNode(i,2) = elementStressAtNodes(fil,col,2);
        stressAtWestNode(i,3) = elementStressAtNodes(fil,col,3);        
    end    
    y(i) = nodesPositionArray(sideWest(i),2);
    stressBending(i,1) = stressAtWestNode(i,1)*(0.5*tSCL-y(i));
    stressBending(i,2) = stressAtWestNode(i,2)*(0.5*tSCL-y(i));
    stressBending(i,3) = stressAtWestNode(i,3)*(0.5*tSCL-y(i));  
    
end


Sxxm = trapz(y,stressAtWestNode(:,1))/tSCL;
Syym = trapz(y,stressAtWestNode(:,2))/tSCL;
Sxym = trapz(y,stressAtWestNode(:,3))/tSCL;
Sxxb = trapz(y,stressAtWestNode(:,1).*(tSCL/2-y))*6/tSCL^2;
% Syyb = trapz(y,stressBending(:,2))*6/tSCL^2; no segun norma
% Sxyb = trapz(y,stressBending(:,3))*6/tSCL^2; no segun norma toma pa vo  pontelli
SxxF = stressAtWestNode(end,1)-(Sxxm-Sxxb);
SyyF = stressAtWestNode(end,2)-(Syym);
SxyF = stressAtWestNode(end,3)-(Sxym);

SvmFEA =sqrt(stressAtWestNode(end,1)^2+stressAtWestNode(end,2)^2+3*stressAtWestNode(end,3)^2); 
Svmm = sqrt(Sxxm^2+Syym^2+3*Sxym^2);
Svmb = Sxxb;
SvmF = sqrt(SxxF^2+SyyF^2+3*SxyF^2);
SvmF2 = SvmFEA -(Svmm-Svmb); %esta da mas chica que la SvmF

stressLine = abs(Sxxb).*(2/tSCL.*y-1)+Sxxm;
figure
plot(y,stressAtWestNode(:,1))
grid on; hold on
plot(y,stressLine)
title('stress xx: FEA and Sm+Sb')
legend('FEA', '\sigma_m+\sigma_b')

Sm = [Sxxm Syym Sxym Svmm];
Sb = [Sxxb Svmb];
Sf = [SxxF SyyF SxyF SvmF];

end

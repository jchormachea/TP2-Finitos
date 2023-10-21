% stress clasification line

%Side west esta descendente en y
%N1 = (L/2-x)*(L-x)*2/L^2
%N2 = x*(L-x)*4/L^2
%N3 = x*(x-L/2)*2/L^2
% int(N1) = L/6 int(N2) = 2L/3 int(N3) = L/6
% int2(N1) = (L*t)/12; int2(N2) = -(L*(L - t))/3; int2(N3) = (L*(2*L - t))/12
clear stressAtWestNode stressAtWestNode2

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
Sxxb = trapz(y,stressBending(:,1))*6/tSCL^2;
Syyb = trapz(y,stressBending(:,2))*6/tSCL^2;
Sxyb = trapz(y,stressBending(:,3))*6/tSCL^2;
x = linspace(0,tSCL,fix(size(sideWest,1)));

stressLine = Sxxm+x.*Sxxb;
figure
plot(x,stressAtWestNode(:,1))
grid on; hold on
plot(x,stressLine)
title('stress FEA and Sm+Sb')
legend('FEA', 'Sm+Sb')



% %integral sigma_m
% sumLs = 0;
% for i = 1:fix(size(sideWest,1)/2)
%     iNode = 1+2*(i-1);
%     iNode2 = iNode+1;
%     iNode3 = iNode+2;
%     node1 = sideWest(iNode);
%     node2 = sideWest(iNode2);
%     node3 = sideWest(iNode3);
%     elementL = abs(nodesPositionArray(node1,2)-nodesPositionArray(node3,2));
%     sumLs = sumLs+elementL;
%     intN1 = elementL/6; %integrales de N para sigma_m
%     intN2 = 2*elementL/3;
%     intN3 = elementL/6;
%     int2N1 = (L*t)/12; %integrales de N*(0.5*t-x) para sigma_b
%     int2N2 = -(L*(L - t))/3; 
%     int2N3 = (L*(2*L - t))/12;
%     weightStressWest(i,1) = weightStressWest(i,1) + (intN1*stressAtWestNode(iNode,1)+intN2*stressAtWestNode(iNode2,1)+intN3*stressAtWestNode(iNode3,1))*elementL;
%     weightStressWest(i,2) = weightStressWest(i,2) + (intN1*stressAtWestNode(iNode,2)+intN2*stressAtWestNode(iNode2,2)+intN3*stressAtWestNode(iNode3,2))*elementL;
%     weightStressWest(i,3) = weightStressWest(i,3) + (intN1*stressAtWestNode(iNode,3)+intN2*stressAtWestNode(iNode2,3)+intN3*stressAtWestNode(iNode3,3))*elementL;
%     weightStressWest2(i,1) = weightStressWest2(i,1) + (int2N1*stressAtWestNode(iNode,1)+int2N2*stressAtWestNode(iNode2,1)+int2N3*stressAtWestNode(iNode3,1))*elementL;
%     weightStressWest2(i,2) = weightStressWest2(i,2) + (int2N1*stressAtWestNode(iNode,2)+int2N2*stressAtWestNode(iNode2,2)+int2N3*stressAtWestNode(iNode3,2))*elementL;
%     weightStressWest2(i,3) = weightStressWest2(i,3) + (int2N1*stressAtWestNode(iNode,3)+int2N2*stressAtWestNode(iNode2,3)+int2N3*stressAtWestNode(iNode3,3))*elementL;
%     
% end
%pruebo integracion numerica por nodo, sin func de forma
function [Sm,Sb,Sf,SmNx,SbNx,SfNx] = SCL(Sy, nodesPositionArray,elementNodesArray, sideWest,elementStressAtNodes,stressNx)

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
end
%escalado
fscale = Sy/max(stressAtWestNode(:,1));
stressAtWestNode = fscale*stressAtWestNode(:,:); %llevo Sx a fluencia.
fscaleNx = Sy/max(stressNx(:,2)); %llevo Sx a fluencia.
stressNx(:,2) = fscaleNx*stressNx(:,2);
stressNx(:,3) = fscaleNx*stressNx(:,3);

%matlab SCL
Sxxm = trapz(y,stressAtWestNode(:,1))/tSCL;
Syym = trapz(y,stressAtWestNode(:,2))/tSCL;
Sxym = trapz(y,stressAtWestNode(:,3))/tSCL;
Sxxb = trapz(y,stressAtWestNode(:,1).*(tSCL/2-y))*6/tSCL^2;
% Syyb = trapz(y,stressBending(:,2))*6/tSCL^2; no segun norma
% Sxyb = trapz(y,stressBending(:,3))*6/tSCL^2; no segun norma toma pa vo  pontelli
SxxF = stressAtWestNode(end,1)-(Sxxm-Sxxb);
SyyF = stressAtWestNode(end,2)-(Syym);
SxyF = stressAtWestNode(end,3)-(Sxym);

SvmFEA = (stressAtWestNode(:,1).^2+stressAtWestNode(:,2).^2+3*stressAtWestNode(:,3).^2).^0.5;
Svmm = sqrt(Sxxm^2+Syym^2+3*Sxym^2);
Svmb = abs(Sxxb);
SvmF = sqrt(SxxF^2+SyyF^2+3*SxyF^2);

stressLine = abs(Sxxb).*(2/tSCL.*y-1)+Sxxm;

%Nx SCL (el corte es nulo por eso no se exportó)
SxxmNx = trapz(stressNx(:,1),stressNx(:,2))/tSCL;
SyymNx = trapz(stressNx(:,1),stressNx(:,3))/tSCL;
SxxbNx = trapz(stressNx(:,1),stressNx(:,2).*(tSCL/2-stressNx(:,1)))*6/tSCL^2;
SxxFNx = stressNx(end,2)-(SxxmNx-SxxbNx);
SyyFNx = stressNx(end,3)-(SyymNx);

SvmFEANx = (stressNx(:,2).^2+stressNx(:,3).^2).^0.5;
SvmmNx = sqrt(SxxmNx^2+SyymNx^2);
SvmbNx = abs(SxxbNx);
SvmFNx = sqrt(SxxFNx^2+SyyFNx^2);

stressLineNx = abs(SxxbNx).*(2/tSCL.*stressNx(:,1)-1)+SxxmNx;


figure
plot(y,stressAtWestNode(:,1),'b')
% plot(y,SvmFEA,'b')
grid on; hold on
plot(stressNx(:,1),stressNx(:,2),'k-.') %tension en x de Nx
plot(y,stressLine,'r')
plot(stressNx(:,1),stressLineNx,'m*')
title('stress xx: FEA and Sm+Sb')
legend('MATLAB FEA', 'Nx FEA', '\sigma_m+\sigma_b MATLAB', '\sigma_m+\sigma_b Nx')
ylim([Sxxb 260]) 

Sm = [Sxxm Syym Sxym Svmm];
Sb = [Sxxb Svmb];
Sf = [SxxF SyyF SxyF SvmF];

SmNx = [SxxmNx SyymNx SvmmNx];
SbNx = [SxxbNx SvmbNx];
SfNx = [SxxFNx SyyFNx SvmFNx];

end

% Case launcher
clc; clear; close all; 
election = 5;
%% Preprocess
tic
elementType='Q8';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
% Mesh generation
t = 2; %mm
w = 50; %mm
d = w/1.1;

% [elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator(elementType,shape,mLength,mHeight,nElementsInLength,nElementsInHeight,ang,distortion);
% nodesData = 'Nodos_prueba.csv';
% conecData = 'Conectividades_prueba.txt';
nodesfiles = {'NodosRD005.csv' 'NodosRD01.csv' 'NodosRD015.csv' 'NodosRD022.csv' 'NodosRD03.csv' 'Nodos_Pieza2.csv' 'Nodos_Pieza2.csv'};
conecfiles = {'ConectividadesRD005.txt' 'ConectividadesRD01.txt' 'ConectividadesRD015.txt' 'ConectividadesRD022.txt' 'ConectividadesRD03.txt' 'Conectividades_Pieza2.txt' 'Conectividades_Pieza2.txt'};
forcesfiles  = {'FuerzaRD005.csv' 'FuerzasRD01.csv' 'FuerzasRD015.csv' 'FuerzaRd022.csv' 'FuerzaRD03.csv' 'Fuerzas_pieza1_ej2'};
rcases = [0.05*d 0.1*d 0.15*d 0.22*d 0.3*d]; 

if election <6
    r = rcases(election);
end
nodesData = nodesfiles{election};
conecData = conecfiles{election};
forceData = forcesfiles{election};
[nodesPositionArray,elementNodesArray] = DataRead(nodesData,conecData);

%ubicacion de lados
%oeste
val = min(nodesPositionArray(:,1));
nodesWest = find(nodesPositionArray(:,1)==val); 
nodesWest = [nodesWest nodesPositionArray(nodesWest,:)];
sideWest = sortrows(nodesWest,3,'ascend'); %ordeno ascendente en y
sideWest = sideWest(:,1); %nodos lado 4(oeste)
%este
val = max(nodesPositionArray(:,1));
nodesEast = find(nodesPositionArray(:,1)==val); 
nodesEast = [nodesEast nodesPositionArray(nodesEast,:)];
sideEast = sortrows(nodesEast,3,'descend'); %ordeno descendente en y
sideEast = sideEast(:,1); %nodos lado 2(este)
%sur(Ojota)
val = min(nodesPositionArray(:,2));
nodesSouth = find(nodesPositionArray(:,2)==val); 
nodesSouth = [nodesSouth nodesPositionArray(nodesSouth,:)];
sideSouth = sortrows(nodesSouth,2,'ascend'); %ordeno ascendente en x
sideSouth = sideSouth(:,1); %nodos lado 1(sur)


nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs


% Material properties
E = 200000;%Mpa
nu = 0.25;
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

if election < 6
        % Boundary conditions
        boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
        boundaryConditionsArray(sideWest,1) = true; %apoyo movil restringe x
        boundaryConditionsArray(sideSouth,1) = true; %apoyo movil restringe x
        boundaryConditionsArray(sideSouth(1),:) = true; %apoyo fijo
        % Load definition
        distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
        pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
        fuerzas = readmatrix(forceData);    
        pointLoadsArray(fuerzas(:,1),1) = fuerzas(:,2);

        q = 1000/(25^2/2); %*y (distribuida lineal con y)
%         for i = 1:fix(size(sideEast,1)/2)%cantidad de elementos en el lado
%             iNode = 1+2*(i-1);
%             iNodeA = sideEast(iNode);
%             iNodeB = sideEast(iNode+1);
%             iNodeC = sideEast(iNode+2);
%             L = abs(nodesPositionArray(iNodeA,2)-nodesPositionArray(iNodeC,2));
%             qA = q*nodesPositionArray(iNodeA,2);
%             qB = q*nodesPositionArray(iNodeB,2);
%             qC = q*nodesPositionArray(iNodeC,2);    
%             directMatrix = L/30*[4 2 -1; 2 16 2; -1 2 4]; 
%             nodalForces = directMatrix*[qA qB qC]';
%             pointLoadsArray([iNodeA iNodeB iNodeC],1) = pointLoadsArray([iNodeA iNodeB iNodeC],1) + nodalForces; %N	
%         end
        
elseif election == 6
    
    % Boundary conditions
        boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
        boundaryConditionsArray(sideWest,1) = true; %apoyo movil restringe x        
        boundaryConditionsArray(sideSouth(1),:) = true; %apoyo fijo
        % Load definition
        distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
        pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
        fuerzas = readmatrix(forceData);    
        pointLoadsArray(fuerzas(:,1),1) = fuerzas(:,2);
        
%         q = 2000/(25^2/2); %*y (distribuida lineal con y)
%         
%         for i = 1:fix(size(sideEast,1)/2)%cantidad de elementos en el lado
%             iNode = 1+2*(i-1);
%             iNodeA = sideEast(iNode);
%             iNodeB = sideEast(iNode+1);
%             iNodeC = sideEast(iNode+2);
%             L = abs(nodesPositionArray(iNodeA,2)-nodesPositionArray(iNodeC,2));
%             qA = q*nodesPositionArray(iNodeA,2);
%             qB = q*nodesPositionArray(iNodeB,2);
%             qC = q*nodesPositionArray(iNodeC,2);    
%             directMatrix = L/30*[4 2 -1; 2 16 2; -1 2 4]; 
%             nodalForces = directMatrix*[qA qB qC]';
%             pointLoadsArray([iNodeA iNodeB iNodeC],1) = pointLoadsArray([iNodeA iNodeB iNodeC],1) + nodalForces; %N	
%         end    
    
else %election == 7 es la misma pieza que el caso 6 pero con otras condiciones de borde por la simetría
    
    % Boundary conditions
    boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
    boundaryConditionsArray(sideWest,1) = true; %apoyo movil restringe x
    boundaryConditionsArray(sideSouth,1) = true; %apoyo movil restringe x
    boundaryConditionsArray(sideSouth(1),:) = true; %apoyo fijo
    % Load definition
    distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
    pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction   
    
    fuerzas = readmatrix(forceData);    
    pointLoadsArray(fuerzas(:,1),1) = fuerzas(:,2);
    
end
    


%% Solver

% Stiffness calculation and assembly
[stiffnessMatrix]=assembleStiffnessMatrix(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,t);

% Matrix reduction
isFixed = reshape(boundaryConditionsArray',1,[])';
isFree = ~isFixed;

% Loads Vector rearrangement
loadsVector = reshape(pointLoadsArray',1,[])';

% Equation solving
displacementsReducedVector = stiffnessMatrix(isFree,isFree)\loadsVector(isFree);

% Reconstruction
displacementsVector = zeros(nTotalDof,1);
displacementsVector(isFree) = displacementsVector(isFree) + displacementsReducedVector;
%% Postprocess
deformations = reshape(displacementsVector,nDimensions,nNodes)';

%Stress recovery
% [elementStressAtNodes]=GaussStressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
[elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% Deformed plot
meshPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)','b','No');

% Stresses plot
bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));
% bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',reshape(displacementsVector,nDimensions,nNodes)');
toc


%% calculo de concentrador (si corresponde)
fprintf('Caso %d\n',election)
smaxNx = [117.41, 93.221, 82.76, 75.229, 70.022]; %Mpa se sacaron de los resultados de Nx
if election < 6
    M = q/3*(25^3+25^3);
    I = t*d^3/12;
    c = 0.5*d;
    So = M*c/I;
    Smaxfea = max(max(elementStressAtNodes(:,:,1)));
    SmaxNx = smaxNx(election);
    Kf = Smaxfea/So;
    KfNx = SmaxNx/So;
    fprintf('w/d = 1.1 r/d = %.3f\n',r/d)
    fprintf('S_o = %.3f S_max matlab = %.3f S_max Nx = %.3f MPa\n', So, Smaxfea, SmaxNx)
    fprintf('Factor de concentración de tensiones FEA Matlab: %.3f\n',Kf)
    fprintf('Factor de concentración de tensiones FEA Nx: %.3f\n',KfNx)
%     fprintf('Factor de concentración de tensiones teórico: %.3f\n',Kfteo) 
    load Wd11PlacaMomento.mat
    load Wd11PlacaMomentoEdM.mat
    figure
    plot(Wd11PlacaMomento(:,1), Wd11PlacaMomento(:,2),'m-.')
    grid on; hold on
    plot(Wd11PlacaMomentoEdM(:,1), Wd11PlacaMomentoEdM(:,2),'c-.')
    plot(r/d, Kf, '*')
end

if election > 5 %casos para SCL

    [Sm,Sb,Sf] = SCL(nodesPositionArray,elementNodesArray, sideWest,elementStressAtNodes);
    
end

%% Grafico de todos los concentradores obtenidos
KFEAMatlab = [0.05 2.454; 0.1 1.937; 0.15 1.710; 0.22 1.574; 0.3 1.467];
KFEANx = [0.05 2.426; 0.1 1.926; 0.15 1.710; 0.22 1.554; 0.3 1.447];
figure
plot(Wd11PlacaMomento(:,1),Wd11PlacaMomento(:,2),'m-.')
grid on; hold on
plot(Wd11PlacaMomentoEdM(:,1), Wd11PlacaMomentoEdM(:,2),'b-.')
plot(KFEAMatlab(:,1),KFEAMatlab(:,2),'r*')
plot(KFEANx(:,1),KFEANx(:,2),'ko')
title('Concentradores')
legend('libro Concentradores', 'libro Shigley', 'Matlab','Nx')

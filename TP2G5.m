% Case launcher
clc; clear; close all; 

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

nodesData = 'NodosRD022.csv';
conecData = 'ConectividadesRD022.txt';
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


meshPlot(elementNodesArray,nodesPositionArray,'b','No');
nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs


% Material properties
E = 200000;%Mpa
nu = 0.25;
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(sideWest,1) = true; %apoyo movil restringe x
boundaryConditionsArray(sideSouth,1) = true; %apoyo movil restringe x
boundaryConditionsArray(sideSouth(1),:) = true; %apoyo fijo
% Load definition
distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction

q = 2000/(25^2/2); %*y (distribuida lineal con y)
% q = 20/50;
for i = 1:fix(size(sideEast,1)/2)%cantidad de elementos en el lado
    iNode = 1+2*(i-1);
    iNodeA = sideEast(iNode);
    iNodeB = sideEast(iNode+1);
    iNodeC = sideEast(iNode+2);
    L = abs(nodesPositionArray(iNodeA,2)-nodesPositionArray(iNodeC,2));
    qA = q*nodesPositionArray(iNodeA,2);
    qB = q*nodesPositionArray(iNodeB,2);
    qC = q*nodesPositionArray(iNodeC,2);    
    directMatrix = L/30*[4 2 -1; 2 16 2; -1 2 4]; 
    nodalForces = directMatrix*[qA qB qC]';
    pointLoadsArray([iNodeA iNodeB iNodeC],1) = pointLoadsArray([iNodeA iNodeB iNodeC],1) + nodalForces; %N	
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

M = q/3*(25^3+25^3);
I = t*d^3/12;
c = 0.5*d;
So = M*c/I;
Smaxfea = max(max(elementStressAtNodes(:,:,1)));
Kf = Smaxfea/So;
fprintf('Factor de concentración de tensiones: %.3f\n',Kf)



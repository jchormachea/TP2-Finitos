% Case launcher
clc; clear; close all; 

%% Preprocess
tic
elementType='Q8';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
% Mesh generation
t = 2; %mm

% [elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator(elementType,shape,mLength,mHeight,nElementsInLength,nElementsInHeight,ang,distortion);
% name = 'Malla1.csv';
% [elementNodesArray,nodesPositionArray,sideWest,sideEast] = ReadData(name);
elementNodesArray = load('conectividades.mat');
elementNodesArray = elementNodesArray.element;
nodesPositionArray = load('nodos.mat');
nodesPositionArray = nodesPositionArray.nodos;
meshPlot(elementNodesArray,nodesPositionArray,'b','No');

%ubicacion de lados
val = min(nodesPositionArray(:,1));
sideWest = find(nodesPositionArray(:,1)==val); %nodos lado 4(oeste)
val = max(nodesPositionArray(:,1));
sideEast = find(nodesPositionArray(:,1)==val); %nodos lado 4(este)

nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

% Material properties
E = 200000;%Mpa
nu = 0.3;
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
boundaryConditionsArray(sideWest,:) = true; %empotrado
% Load definition
distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction

% % pointLoadsArray (sideNodes(2,:),1) = (-1:2/(size(sideNodes(2,:),2)-1):1)/3;
% % pointLoadsArray (sideNodes(2,2:2:end),1) = 2*pointLoadsArray (sideNodes(2,2:2:end),1);
% % pointLoadsArray (vertexNodes(2),1) = -1/6;
% % pointLoadsArray (vertexNodes(4),1) =  1/6;
pointLoadsArray(sideEast,1) = 1000; %N
pointLoadsArray(sideEast,end) = -1000; %N
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
%Stress recovery
[elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);

% Deformed plot
% meshPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)','b','Yes');

% Stresses plot
bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));
toc
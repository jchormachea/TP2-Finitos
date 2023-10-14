% Case launcher
clc; clear; close all; 

%% Preprocess

elementType='Q4';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
loadCase='Custom2';         %'Uniform' 'Constant Bending' 'Variable Bending' 'Custom'

% Mesh generation
t = 1; %mm
shape = 'Straight'; %'Curved' 'Straight'
mLength = 3; %in
mHeight = 2; %in
nElementsInLength =1;
nElementsInHeight = 2;
ang = 90;
distortion = 0;

% [elementNodesArray,nodesPositionArray,vertexNodes,sideNodes]=quadrilateralDomainMeshGenerator(elementType,shape,mLength,mHeight,nElementsInLength,nElementsInHeight,ang,distortion);

nodesPositionArray = [  -1 -1
                        0 -2
                        1 -1
                        -1 0
                        0 0
                        1 0
                        -1 1
                        0 2
                        1 1];
elementNodesArray = [1 2 5 4
                    2 3 6 5
                    4 5 8 7
                    5 6 9 8];

vertexNodes = [1 2 3 7 8 9];
sideNodes = [1 2 3
            3 6 9
            7 8 9
            1 4 7];

% elementNodesArray

nElements=size(elementNodesArray,1);    %Number of elements
nNodes=size(nodesPositionArray,1);      %Number of nodes
nTotalDof=nNodes*nDimensions;           %Number of total dofs

% Material properties
% E = 200;%Mpa
% nu = 1/3;
% E = 30e6; %psi
% nu = 0.25;
E = 1;
nu = 1/3;
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);

% Boundary conditions
boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
% boundaryConditionsArray(vertexNodes(1),1:2) = true;
% boundaryConditionsArray(sideNodes(end,find(sideNodes(end,:))),1) = true;
% boundaryConditionsArray(sideNodes(4,[1 2 3]),:) = true; %Custom Q4
% boundaryConditionsArray(sideNodes(4,[1 2 3 4 5]),:) = true; %Custom Q8,9
% boundaryConditionsArray(vertexNodes([1,3]),:) = true; %Custom 3
% boundaryConditionsArray(vertexNodes(2),2) = true; %Custom 3
boundaryConditionsArray([1 3 4 6 7 9],:) = true; %Custom 2
% Load definition
distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction

% Load case load imposition
switch loadCase
    case 'Uniform'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),1) = 1;
                pointLoadsArray (vertexNodes(2),1) = 1/2;
                pointLoadsArray (vertexNodes(4),1) = 1/2;
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,1:end),1) = 1/3;
                pointLoadsArray (sideNodes(2,2:2:end),1) = 2/3;
                pointLoadsArray (vertexNodes(2),1) = 1/6;
                pointLoadsArray (vertexNodes(4),1) = 1/6;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/abs(sum(pointLoadsArray(:,1)));
    case 'Constant Bending'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),1) = -1:2/(size(sideNodes(2,:),2)-1):1;
                pointLoadsArray (vertexNodes(2),1) = -1/2;
                pointLoadsArray (vertexNodes(4),1) =  1/2;
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,:),1) = (-1:2/(size(sideNodes(2,:),2)-1):1)/3;
                pointLoadsArray (sideNodes(2,2:2:end),1) = 2*pointLoadsArray (sideNodes(2,2:2:end),1);
                pointLoadsArray (vertexNodes(2),1) = -1/6;
                pointLoadsArray (vertexNodes(4),1) =  1/6;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/sum(pointLoadsArray(:,1).*nodesPositionArray(:,2));
    case 'Variable Bending'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray (sideNodes(2,:),2) = -((-1:2/(size(sideNodes(2,:),2)-1):1).*(-1:2/(size(sideNodes(2,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(4,:),2) = ((-1:2/(size(sideNodes(4,:),2)-1):1).*(-1:2/(size(sideNodes(4,:),2)-1):1)+1);
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                % Right side shear
                pointLoadsArray (sideNodes(2,:),2) = -((-1:2/(size(sideNodes(2,:),2)-1):1).*(-1:2/(size(sideNodes(2,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(2,2:2:end),2) = 2*pointLoadsArray (sideNodes(2,2:2:end),1);
                pointLoadsArray (vertexNodes(2),1) = 0;
                pointLoadsArray (vertexNodes(4),1) =  0;
                % Left side shear
                pointLoadsArray (sideNodes(4,:),2) = ((-1:2/(size(sideNodes(4,:),2)-1):1).*(-1:2/(size(sideNodes(4,:),2)-1):1)+1);
                pointLoadsArray (sideNodes(4,2:2:end),2) = 2*pointLoadsArray (sideNodes(4,2:2:end),1);
                pointLoadsArray (vertexNodes(1),1) = 0;
                pointLoadsArray (vertexNodes(3),1) =  0;
        end
        % Scaling
        pointLoadsArray=pointLoadsArray/abs(sum(pointLoadsArray(sideNodes(4,:),2)))/10;
    case 'Custom'
        switch elementType
            case {'CST' 'Q4'}
%                 pointLoadsArray (sideNodes(2,2),1) = -1;
                pointLoadsArray (sideNodes(3,:),2) = -300*t*mLength/2;%300psi¨*t*L (0,5in)
                pointLoadsArray (vertexNodes(4),2) = pointLoadsArray (vertexNodes(4),2) -1000; %lb
%                 pointLoadsArray (sideNodes(2,:),1) = 1;
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,3),1) = -1;
        end
        
        
    case 'Custom2'
        switch elementType
            case {'CST' 'Q4'}
                pointLoadsArray ([2,8],2) = -1;
                
        end
        switch elementType
            case {'LST' 'Q8' 'Q9'}
                pointLoadsArray (sideNodes(2,3),1) = -1;
        end       
        
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
displacementsVector
%% Postprocess
%Stress recovery
[elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);

% Deformed plot
meshPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)','b','Yes');

% Stresses plot
% bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));

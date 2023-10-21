function [elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,nodalDisplacements)
% Stress recovery from displacements
% 
% [elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacements)
%
% elementStressAtNodes:      Element stresses at element nodes
%
% elementType:          Type of element 'CST' LST' 'Q4' 'Q8' 'Q9'
% nodesPositionArray:   Nodal position in cartesian coordinates
% elementNodesArray:    Element conectivity matrix
% constitutiveMatrix:   Constitutive Matrix
%

%% Definitions
switch elementType
    case 'CST'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates=[0 0
                                 1 0
                                 0 1];
    case 'LST'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates=[0   0
                                 1   0
                                 0   1
                                 0.5 0
                                 0.5 0.5
                                 0   1];
    case 'Q4'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1];
    case 'Q8'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1
                                    0 -1
                                    1  0
                                    0  1
                                   -1  0];
    case 'Q9'
        nDimensions = 2;                                %Problem dimension
        naturalNodalCoordinates = [-1 -1
                                    1 -1
                                    1  1
                                   -1  1
                                    0 -1
                                    1  0
                                    0  1
                                   -1  0
                                    0  0];
end

nElements=size(elementNodesArray,1);            %Number of elements

nNodes=size(nodesPositionArray,1);              %Number of nodes
nElementalNodes = size(elementNodesArray,2);    %Number of node in each element
nElementalDof = nDimensions*nElementalNodes;    %Number of elemental Dofs

nTotalDof = nDimensions*nNodes;                 %Number of node in each element

nPoints = size(naturalNodalCoordinates,1);      %Number of points to get stress


%% Stress recovery

% Shape functions derivatives at every point in natural coordinates
shapeFunctionsDerivatives = getShapeFunctionsDerivatives(naturalNodalCoordinates,elementType);

for iElement = 1:nElements
    
    elementalDofs = convertNode2Dof(elementNodesArray(iElement,:),nDimensions);
    elementalNodesPosition = nodesPositionArray(elementNodesArray(iElement,:),:);
    for iPoint = 1:nPoints
        
        % Jacobian calculation at nodal positions
        jacobian = shapeFunctionsDerivatives(:,:,iPoint)*elementalNodesPosition;
        jacobianDeterminant=det(jacobian);
        
        % Shape functions derivatives in structural coordinates
        structuralShapeFunctionsDerivatives = jacobian\shapeFunctionsDerivatives(:,:,iPoint);
        
        B = zeros(size(constitutiveMatrix,1),nElementalDof);
        B(1,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(1,:);
        B(2,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(2,:);
        B(3,1:2:nElementalDof-1) = structuralShapeFunctionsDerivatives(2,:);
        B(3,2:2:nElementalDof) = structuralShapeFunctionsDerivatives(1,:);
        
        elementStressAtNodes(iElement,iPoint,:) = constitutiveMatrix*B*nodalDisplacements(elementalDofs);
    end
end


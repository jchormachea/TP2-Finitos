function ShapeFunctionCreatorAndViewer(shapeFunctionToShow)
% Symbolic Shape Funtions Generator and Viewer
%
% shapeFunctionToShow: Node to show shape function
%
% ShapeFunctionCreatorAndViewer(1)
% 

close all

%% Symbolic Shape Function Generator
syms xi eta real

% Functionality
functionality=[1 xi eta xi^2 xi*eta eta^2 xi^2*eta xi*eta^2];

% Natural coordinates nodal location
naturalNodes=[-1 -1
               1 -1
               1  1
              -1  1
               0 -1
               1  0
               0  1
              -1  0];
          
          %Control
          if size(functionality,2)==size(naturalNodes,1)
              
              nNodes=size(functionality,2);
          else
              disp('Inconsistent functionality and natural coordinates nodes locations')
              return
          end

for iNode=1:nNodes;    
    A(iNode,:)=subs(functionality,{xi,eta},{naturalNodes(iNode,1),naturalNodes(iNode,2)});
end

N=functionality*inv(A);

disp(N')

% Testing

subs(N,{xi,eta},{naturalNodes(:,1),naturalNodes(:,2)})

% Shape Function Plot
figure
ezsurf(N(shapeFunctionToShow),[-1 1 -1 1]); hold on

%% Isoparametric structural element viewer

[xiNaturalCoordinate,etaNaturalCoordinate] = meshgrid(-1:.1:1, -1:.1:1);
nPointsX=size(xiNaturalCoordinate,1); nPointsY=size(etaNaturalCoordinate,1);

structuralNodesCoordinates=[-10 -10
                             10 -10
                             10  10
                            -10  10
                              0 -10
                             10   0   %Use this
                              0  10
                            -10   0];

% Function creation for evaluation                        
                        
shapeFunctionsDerivatives=[diff(N,xi)
                           diff(N,eta)];

handleShapeFunctionDerivatives=matlabFunction(shapeFunctionsDerivatives);

handleShapeFunction=matlabFunction(N);
                       
   for iPointsX=1:nPointsX
       for iPointsY=1:nPointsY;

           % Numerical evaluation
           evaluatedShapeFunction=handleShapeFunction(xiNaturalCoordinate(iPointsX,iPointsY),etaNaturalCoordinate(iPointsX,iPointsY));
           
           evaluatedShapeFunctionDerivatives=handleShapeFunctionDerivatives(xiNaturalCoordinate(iPointsX,iPointsY),etaNaturalCoordinate(iPointsX,iPointsY));
           
           % Points Location
           
           Xiso(iPointsX,iPointsY)=evaluatedShapeFunction*structuralNodesCoordinates(:,1);
           Yiso(iPointsX,iPointsY)=evaluatedShapeFunction*structuralNodesCoordinates(:,2);
           
           % Jacobian
 
           J(iPointsX,iPointsY)=det(evaluatedShapeFunctionDerivatives*structuralNodesCoordinates);
           
       end
   end

figure;
surfc(Xiso,Yiso,ones(nPointsX,nPointsY),J); xlabel('x'); ylabel('y'); colorbar;

view([0,0,1])


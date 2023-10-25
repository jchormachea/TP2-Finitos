%% Código TP2 FEA Hormachea Jose Conrado 61439 - Nieto Franco 61459
% este código se hizo con matlab R2020a. Puede que no corra si se utiliza otra versión

%% Inciliazar código
clc; clear; close all; 
election = 7;
% % valores de election según cada caso
% % 1: concentrador con r/d = 0.05
% % 2: concentrador con r/d = 0.1
% % 3: concentrador con r/d = 0.15
% % 4: concentrador con r/d = 0.22
% % 5: concentrador con r/d = 0.3
% % 6: SCL para la placa 1 con w/d = 2 y r/d = 0.15
% % 7: SCL para la placa 2 con w/d = 2 y r/d = 0.15
%% Preprocess

elementType='Q8';          %'CST' 'LST' 'Q4' 'Q8' 'Q9'
problemType='Stress';       %'Stress' 'Strain' 'Axisymmetric'
nDimensions=2;              %Problem dimension
% Mesh generation
t = 2; %mm
w = 50; %mm
d = w/1.1; % solo se usa para election de 1 a 5

%nombres de archivos
nodesfiles = {'NodosRD005.csv' 'NodosRD01.csv' 'NodosRD015.csv' 'NodosRD022.csv' 'NodosRD03.csv' 'Nodos_pieza1_ej2.csv' 'Nodos_pieza2_ej2.csv'};
conecfiles = {'ConectividadesRD005.txt' 'ConectividadesRD01.txt' 'ConectividadesRD015.txt' 'ConectividadesRD022.txt' 'ConectividadesRD03.txt' 'ConectividadesPieza1EJ2.txt' 'ConectividadesPieza2EJ2.txt'};
forcesfiles  = {'FuerzaRD005.csv' 'FuerzasRD01.csv' 'FuerzasRD015.csv' 'FuerzaRd022.csv' 'FuerzaRD03.csv' 'Fuerza_pieza1_ej2.csv' 'Fuerza_pieza2_ej2.csv'};
rcases = [0.05*d 0.1*d 0.15*d 0.22*d 0.3*d]; 

if election <6
    r = rcases(election); %para calcular Kf
end

%lectura y eleccion de archivos
nodesData = nodesfiles{election};
conecData = conecfiles{election};
forceData = forcesfiles{election};
[nodesPositionArray,elementNodesArray] = DataRead(nodesData,conecData);
%grafico 
% meshPlot(elementNodesArray,nodesPositionArray,'b','No');

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
%sur
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
Sy = 250; %Mpa
[constitutiveMatrix] = constitutiveIsotropicMatrix(problemType,E,nu);
%definicion de carga y condiciones de borde

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

        q = 1000/(25^2/2); %*y (distribuida lineal con y) para calcular sigma_o
        
elseif election == 7 %pieza con simetria simple
    
        % Boundary conditions
        boundaryConditionsArray = false(nNodes,nDimensions);    % Boundary conditions array true=fixed
        boundaryConditionsArray(sideWest,1) = true; %apoyo movil restringe x        
        boundaryConditionsArray(sideSouth(1),:) = true; %apoyo fijo
        % Load definition
        distributedLoadsArray = zeros(nNodes,nDimensions);      % Distributed load nodal value for each direction
        pointLoadsArray = zeros(nNodes,nDimensions);            % Point load nodal value for each direction
        fuerzas = readmatrix(forceData);    
        pointLoadsArray(fuerzas(:,1),1) = fuerzas(:,2);
        
else %election == 6 es la pieza del calculo de concentradores pero con otras relaciones
    
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
[elementStressAtNodes]=stressRecovery(elementType,elementNodesArray,nodesPositionArray,constitutiveMatrix,displacementsVector);
% Deformed plot
% meshPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)','b','No');

% Stresses plot
bandPlot(elementNodesArray,nodesPositionArray+reshape(displacementsVector,nDimensions,nNodes)',squeeze(elementStressAtNodes(:,:,1)));


%% calculo de concentrador y SCL(según corresponda)

fprintf('Caso %d\n',election)
smaxNx = [117.41, 93.221, 82.76, 75.229, 70.022]; %Mpa se sacaron de los resultados de Nx
% smaxMb = [118.776 93.756 82.765 76.189 71.006]; %Mpa se sacaron de matlab
load('Wd11PlacaMomentoEdM.mat');
load('Wd11PlacaMomento.mat');

if election < 6 %calculo concentradores
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
    figure
    plot(Wd11PlacaMomento(:,1), Wd11PlacaMomento(:,2),'m-.')
    grid on; hold on
    plot(Wd11PlacaMomentoEdM(:,1), Wd11PlacaMomentoEdM(:,2),'c-.')
    plot(r/d, Kf, '*')
    plot(r/d, KfNx, '*')
    legend('Pilkey', 'Budynas','MATLAB','Nx')
end

if election > 5 %casos para SCL
    if election == 6
        sxNx = csvread('Sxx_pieza1_ej2.csv',28);
%         sxNx = sxNx(:,[1 2])
        syNx = csvread('Syy_pieza1_ej2.csv',28);  
        
    else %election ==7
        sxNx = csvread('Sxx_pieza2_ej2.csv',28);
%       sxNx = sxNx(:,[1 2])
        syNx = csvread('Syy_pieza2_ej2.csv',28);  
    end
    stressNx = [sxNx(:,1) sxNx(:,2) syNx(:,2)]; %[pos sigmaX sigmaY] el corte es nulo y no se exportó
    [Sm,Sb,Sf,SmNx,SbNx,SfNx] = SCL(Sy, nodesPositionArray,elementNodesArray, sideWest,elementStressAtNodes,stressNx);
    results = [Sm(end) SmNx(end)
                Sb(end) SbNx(end)
                Sf(end) SfNx(end)];  
    Tension = {'\sigma_{vm,m} [MPa]';'\sigma_{vm,b} [MPa]';'\sigma_{vm,F} [MPa]'};
    MATLAB = [results(:,1)];
    Nx = [results(:,2)];
    tabla = table(Tension, MATLAB, Nx)
end

%% Grafico de todos los concentradores obtenidos (copiados para no tener que correr todo siempre)

% KFEAMatlab = [0.05 2.454; 0.1 1.937; 0.15 1.710; 0.22 1.574; 0.3 1.467];
% KFEANx = [0.05 2.426; 0.1 1.926; 0.15 1.710; 0.22 1.554; 0.3 1.447];
% 
% figure
% plot(Wd11PlacaMomento(:,1),Wd11PlacaMomento(:,2),'m-.')
% grid on; hold on
% plot(Wd11PlacaMomentoEdM(:,1), Wd11PlacaMomentoEdM(:,2),'b-.')
% plot(KFEAMatlab(:,1),KFEAMatlab(:,2),'r*')
% plot(KFEANx(:,1),KFEANx(:,2),'ko')
% title('Concentradores')
% legend('Pilkey', 'Budynas', 'MATLAB','Nx')

%Interpolacion de elements
clear; close all;

%% Generaci�n Simb�lica de funciones de forma
syms x y real

%Difici�n de funciones de forma
X=[1 x y x*y x^2 y^2 x^2*y x*y^2 x^2*y^2];
X= X(1:8); % Q8
 
nodes=[0 0 
       2 0
       2 2
       0 2
       1 0
       2 1
       1 2
       0 1];
elements = [1 2
            2 3
            3 5
            5 8
            8 7
            7 6
            1 4
            4 6];     
numberOfNodes=size(nodes,1);

if numberOfNodes ~= length(X)
    warning('Faltan terminos en la funcionalidad')
end

dofPerNode = 2;

%%Armado de la matriz A

for i=1:numberOfNodes   
    A(i,:)=subs(X,{x,y},{nodes(i,1),nodes(i,2)});
end

N=X*inv(A); %% Uso inv porque \ no funciona con simbolico
fprintf('Las funciones de forma\n')
disp(N')

fprintf('Las funciones de forma valen 1 en el nodo y 0 en los demas?\n')
subs(N,{x,y},{nodes(:,1),nodes(:,2)})

%% Ploteo Funciones de Forma
figure(1)
ezsurf(N(3),[0 2 0 2]); hold on

%text(nodes(:,1),nodes(:,2),-ones(numberOfNodes,1)*0.25,num2str([1:numberOfNodes]'))

%% Element Formulation    
%Ordeno la matriz de funciones de forma para armar B
for i=1:numberOfNodes
    orderedN(1,i*2-1)=N(i);
    orderedN(2,i*2)=N(i);
end

B=[diff(orderedN(1,:),x)
   diff(orderedN(2,:),y)
   diff(orderedN(1,:),y)+diff(orderedN(2,:),x)];


%% Stifness Matrix
E = 200000;%Mpa
nu = 0.25;
% C=E/(1+nu)*[(1-nu)/(1-2*nu)     nu/(1-2*nu) 0
%                 nu/(1-2*nu) (1-nu)/(1-2*nu) 0
%                           0               0 0.5];
                      
C = E/(1 - nu^2)*[ 1.0     nu         0.0
                    nu    1.0         0.0
                   0.0    0.0     (1 - nu)/2 ];

               
%% Integro para armar la matriz de rigidez. VER LIMITES DE INTEGRACION
K=int(int(B'*C*B,y,0,2),x,0,2);

R = zeros(numberOfNodes,dofPerNode);        % Vector de cargas
R(3,1) = -1000;
R(2,1) = 1000;

Rr = reshape(R',[],1);

bc = false(numberOfNodes,dofPerNode);       % Matriz de condiciones de borde
bc([1 8 4],:) = true;
isFixed = reshape(bc',[],1);
isFree = ~isFixed;
reducedDisplacements = K(isFree,isFree)\Rr(isFree);

% Reconstrucci�n
D = zeros(numberOfNodes*dofPerNode,1);
D(isFree) = D(isFree) + reducedDisplacements;

% Reacciones
Rv = K(isFixed,isFree)*D(isFree);
reacciones=zeros(dofPerNode*numberOfNodes,1);
reacciones(isFixed) = Rv;
reacciones = (reshape(reacciones,dofPerNode,[]))';

% %Recuperaci�n de tensiones
stress = C*B*D; % Sx,Sy,tau_xy

%%Plots 
% Configuraci�n deformada
D = (reshape(D,dofPerNode,[]))';
nodePosition = nodes + D(:,1:2);

% %Graficación
% bandplot(elements,nodePosition,stress(:,:,2),[-1 1],'k');
figure(2)
hold on
meshPlot([1 2 3 4],nodePosition,'r')
meshPlot([1 2 3 4],nodes,'k')



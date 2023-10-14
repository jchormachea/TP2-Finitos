clc
clear all
close all

forma=5; %1 a 6
if forma==1 %Q4
    eleType=4;
    nodos=[2,2;
           0,2
           0,0;
           2,0;];
elseif forma==2 %Q4
    eleType=4;
    nodos=[0,0;
        1,-1.6;
        2,1.2;
        0.6,1];
elseif forma==3 %Q4
    eleType=4;
    nodos=[0,0;
        0.5,-5.5;
        1,1;
        -0.5,2.5];
elseif forma==4 %Q4
    eleType=4;
    nodos=[-0.5,2.5;
        0.25,0.25;
        0.5,-2.5;
        1.0,1.0];
elseif forma==5 %Q8
    eleType=8;
    nodos=[0,0;
        0.05,-0.1;
        0.2,-0.1;
        0.2,0.2;
        0.025,-0.05;
        0.15,-0.2;
        0.2,0.1;
        0.1,0.1];
elseif forma==6 % Q9
    eleType=9;
    nodos=[-1,-1;
        1,-1;
        1,1;
        -1,1;
        0,-1;
        1,0;
        0,1;
        -1,0;
        -0.9,0];
end

%% Funciones de forma
if eleType==4
    N4 = @(ksi,eta) 0.25*(1 - ksi)*(1 + eta);
    N3 = @(ksi,eta) 0.25*(1 + ksi)*(1 + eta);
    N2 = @(ksi,eta) 0.25*(1 + ksi)*(1 - eta);
    N1 = @(ksi,eta) 0.25*(1 - ksi)*(1 - eta);
    % Derivadas de las funciones de forma respecto de ksi, eta
    dN(:,:) = @(ksi,eta) [  % derivadas respecto de ksi
        -0.25*(1 - eta),  0.25*(1 - eta), 0.25*(1 + eta), -0.25*(1 + eta)
        % derivadas respecto de eta
        -0.25*(1 - ksi), -0.25*(1 + ksi), 0.25*(1 + ksi),  0.25*(1 - ksi) ];
elseif eleType==8
    N8 = @(ksi,eta) 0.50*(1 - ksi  )*(1 - eta^2);
    N7 = @(ksi,eta) 0.50*(1 - ksi^2)*(1 + eta  );
    N6 = @(ksi,eta) 0.50*(1 + ksi  )*(1 - eta^2);
    N5 = @(ksi,eta) 0.50*(1 - ksi^2)*(1 - eta  );
    N4 = @(ksi,eta) 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7(ksi,eta) + N8(ksi,eta));
    N3 = @(ksi,eta) 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6(ksi,eta) + N7(ksi,eta));
    N2 = @(ksi,eta) 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5(ksi,eta) + N6(ksi,eta));
    N1 = @(ksi,eta) 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5(ksi,eta) + N8(ksi,eta));
    % Derivadas de las funciones de forma respecto de ksi, eta
    dN(:,:) = @(ksi,eta) [  % derivadas respecto de ksi
        -0.25*(-1+eta)*(eta+2*ksi),  -0.25*(-1+eta)*(-eta+2*ksi),    0.25*(1+eta)*(eta+2*ksi),   0.25*(1+eta)*(-eta+2*ksi),...
        ksi*(-1+eta),        -0.5*(-1+eta)*(1+eta),                -ksi*(1+eta),        0.5*(-1+eta)*(1+eta)
        % derivadas respecto de eta
        -0.25*(-1+ksi)*(ksi+2*eta),   -0.25*(1+ksi)*(ksi-2*eta),    0.25*(1+ksi)*(ksi+2*eta),   0.25*(-1+ksi)*(ksi-2*eta),...
        0.5*(-1+ksi)*(1+ksi),                -(1+ksi)*eta,       -0.5*(-1+ksi)*(1+ksi),               (-1+ksi)*eta ];
elseif eleType==9
    N9 = @(ksi,eta)      (1 - ksi^2)*(1 - eta^2);
    N8 = @(ksi,eta) 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9(ksi,eta);
    N7 = @(ksi,eta) 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9(ksi,eta);
    N6 = @(ksi,eta) 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9(ksi,eta);
    N5 = @(ksi,eta) 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9(ksi,eta);
    N4 = @(ksi,eta) 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7(ksi,eta) + N8(ksi,eta) + 0.5*N9(ksi,eta));
    N3 = @(ksi,eta) 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6(ksi,eta) + N7(ksi,eta) + 0.5*N9(ksi,eta));
    N2 = @(ksi,eta) 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5(ksi,eta) + N6(ksi,eta) + 0.5*N9(ksi,eta));
    N1 = @(ksi,eta) 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5(ksi,eta) + N8(ksi,eta) + 0.5*N9(ksi,eta));
    dN = @(ksi,eta) [ % derivadas respecto de ksi
        0.25*eta*(-1+eta)*(2*ksi-1),      0.25*eta*(-1+eta)*(2*ksi+1),       0.25*eta*(1+eta)*(2*ksi+1),...
        0.25*eta*( 1+eta)*(2*ksi-1),                -ksi*eta*(-1+eta),  -1/2*(-1+eta)*(1+eta)*(2*ksi+1),...
        -ksi*eta*(1+eta),  -1/2*(-1+eta)*(1+eta)*(2*ksi-1),           2*ksi*(-1+eta)*(1+eta)
        % derivadas respecto de eta
        0.25*ksi*(-1+2*eta)*(ksi-1),      0.25*ksi*(-1+2*eta)*(1+ksi),       0.25*ksi*(2*eta+1)*(1+ksi),...
        0.25*ksi*(2*eta+1)*(ksi-1),  -0.5*(ksi-1)*(1+ksi)*(-1+2*eta),                 -ksi*eta*(1+ksi),...
        -0.5*(ksi-1)*(1+ksi)*(2*eta+1),                 -ksi*eta*(ksi-1),           2*(ksi-1)*(1+ksi)*eta ];
end

%% Cálculo del jacobiano
% Elegir puntos de evaluación
Ksi=[0.3,-0.9,0.8,0.25,-0.75];
Eta=[0,-0.9,0,-0.5,0.75];
for p=1:length(Ksi)
    % Derivadas de x,y, respecto de ksi, eta
    jac = dN(Ksi(p),Eta(p))*nodos;
    % Derivadas de las funciones de forma respecto de x,y.
    dNxy = jac\dN(Ksi(p),Eta(p));          % dNxy = inv(jac)*dN
    djac = det(jac);
    str = {[' |',num2str(jac(1,1)),' ',num2str(jac(1,2)),'| '],[' |',num2str(jac(2,1)),' ',num2str(jac(2,2)),'| '],[' det(jac) = ',num2str(djac)]};
    if eleType==4
        X = N1(Ksi(p),Eta(p))*nodos(1,1)+N2(Ksi(p),Eta(p))*nodos(2,1)+N3(Ksi(p),Eta(p))*nodos(3,1)+N4(Ksi(p),Eta(p))*nodos(4,1);
        Y = N1(Ksi(p),Eta(p))*nodos(1,2)+N2(Ksi(p),Eta(p))*nodos(2,2)+N3(Ksi(p),Eta(p))*nodos(3,2)+N4(Ksi(p),Eta(p))*nodos(4,2);
    elseif eleType==8
        X = N1(Ksi(p),Eta(p))*nodos(1,1)+N2(Ksi(p),Eta(p))*nodos(2,1)+N3(Ksi(p),Eta(p))*nodos(3,1)+N4(Ksi(p),Eta(p))*nodos(4,1)+...
            N5(Ksi(p),Eta(p))*nodos(5,1)+N6(Ksi(p),Eta(p))*nodos(6,1)+N7(Ksi(p),Eta(p))*nodos(7,1)+N8(Ksi(p),Eta(p))*nodos(8,1);
        Y = N1(Ksi(p),Eta(p))*nodos(1,2)+N2(Ksi(p),Eta(p))*nodos(2,2)+N3(Ksi(p),Eta(p))*nodos(3,2)+N4(Ksi(p),Eta(p))*nodos(4,2)+...
            N5(Ksi(p),Eta(p))*nodos(5,2)+N6(Ksi(p),Eta(p))*nodos(6,2)+N7(Ksi(p),Eta(p))*nodos(7,2)+N8(Ksi(p),Eta(p))*nodos(8,2);
    elseif eleType==9
        X = N1(Ksi(p),Eta(p))*nodos(1,1)+N2(Ksi(p),Eta(p))*nodos(2,1)+N3(Ksi(p),Eta(p))*nodos(3,1)+N4(Ksi(p),Eta(p))*nodos(4,1)+...
            N5(Ksi(p),Eta(p))*nodos(5,1)+N6(Ksi(p),Eta(p))*nodos(6,1)+N7(Ksi(p),Eta(p))*nodos(7,1)+N8(Ksi(p),Eta(p))*nodos(8,1)+N9(Ksi(p),Eta(p))*nodos(9,1);
        Y = N1(Ksi(p),Eta(p))*nodos(1,2)+N2(Ksi(p),Eta(p))*nodos(2,2)+N3(Ksi(p),Eta(p))*nodos(3,2)+N4(Ksi(p),Eta(p))*nodos(4,2)+...
            N5(Ksi(p),Eta(p))*nodos(5,2)+N6(Ksi(p),Eta(p))*nodos(6,2)+N7(Ksi(p),Eta(p))*nodos(7,2)+N8(Ksi(p),Eta(p))*nodos(8,2)+N9(Ksi(p),Eta(p))*nodos(9,2);
    end
    hold on
    line(X,Y,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k')
    text(X,Y,str,'FontSize',12)
end

% Líneas eta fijo
figure(1)
hold on
for eta=[-1:0.1:1]
    x=[];
    y=[];
    for ksi=[-1:0.1:1]
        if eleType==4
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)];
        elseif eleType==8
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)+...
                N5(ksi,eta)*nodos(5,1)+N6(ksi,eta)*nodos(6,1)+N7(ksi,eta)*nodos(7,1)+N8(ksi,eta)*nodos(8,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)+...
                N5(ksi,eta)*nodos(5,2)+N6(ksi,eta)*nodos(6,2)+N7(ksi,eta)*nodos(7,2)+N8(ksi,eta)*nodos(8,2)];
        elseif eleType==9
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)+...
                N5(ksi,eta)*nodos(5,1)+N6(ksi,eta)*nodos(6,1)+N7(ksi,eta)*nodos(7,1)+N8(ksi,eta)*nodos(8,1)+N9(ksi,eta)*nodos(9,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)+...
                N5(ksi,eta)*nodos(5,2)+N6(ksi,eta)*nodos(6,2)+N7(ksi,eta)*nodos(7,2)+N8(ksi,eta)*nodos(8,2)+N9(ksi,eta)*nodos(9,2)];
        end
    end
    if eta==-1 || eta==1
        line(x,y,'Color','k','LineWidth',2)
    else
        line(x,y,'Color','r')
    end
end

% Líneas ksi fijo
figure(1)
hold on
for ksi=[-1:0.1:1]
    x=[];
    y=[];
    for eta=[-1:0.1:1]
        if eleType==4
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)];
        elseif eleType==8
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)+...
                N5(ksi,eta)*nodos(5,1)+N6(ksi,eta)*nodos(6,1)+N7(ksi,eta)*nodos(7,1)+N8(ksi,eta)*nodos(8,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)+...
                N5(ksi,eta)*nodos(5,2)+N6(ksi,eta)*nodos(6,2)+N7(ksi,eta)*nodos(7,2)+N8(ksi,eta)*nodos(8,2)];
        elseif eleType==9
            x = [x,N1(ksi,eta)*nodos(1,1)+N2(ksi,eta)*nodos(2,1)+N3(ksi,eta)*nodos(3,1)+N4(ksi,eta)*nodos(4,1)+...
                N5(ksi,eta)*nodos(5,1)+N6(ksi,eta)*nodos(6,1)+N7(ksi,eta)*nodos(7,1)+N8(ksi,eta)*nodos(8,1)+N9(ksi,eta)*nodos(9,1)];
            y = [y,N1(ksi,eta)*nodos(1,2)+N2(ksi,eta)*nodos(2,2)+N3(ksi,eta)*nodos(3,2)+N4(ksi,eta)*nodos(4,2)+...
                N5(ksi,eta)*nodos(5,2)+N6(ksi,eta)*nodos(6,2)+N7(ksi,eta)*nodos(7,2)+N8(ksi,eta)*nodos(8,2)+N9(ksi,eta)*nodos(9,2)];
        end
    end
    if ksi==-1 || ksi==1
        line(x,y,'Color','k','LineWidth',2)
    else
        line(x,y,'Color','b')
    end
end
plot(nodos(:,1),nodos(:,2),'o','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','b')
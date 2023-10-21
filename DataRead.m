function [nodes,conectivities] = DataRead(nodesData,conecData)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
nodos = readmatrix(nodesData);
nodos(:,[1,4]) = '';%me quedo solo con x e y
iter = 1;
fid = fopen(conecData,'r');
element = [];%vector
while ~feof(fid)    
    linea = fgetl(fid);
    if iter >= 16
        values = regexp(linea,'\d*','Match');
        if isempty(values)
            break
        end            
        element(str2double(values{1}),:) = [str2double(values{2}) str2double(values{3}) str2double(values{4}) str2double(values{5}) str2double(values{6}) str2double(values{7}) str2double(values{8}) str2double(values{9})];
    end
    iter = iter+1;
    if iter >20000
        error('La lectura de archivos falló');
    end
end
fclose(fid);

for i = 1:size(nodos,1)    
   if abs(nodos(i,1)) < 1e-3
       nodos(i,1) = 0;
   end
   if abs(nodos(i,2)) < 1e-3
       nodos(i,2) = 0;
   end    
end

nodes = nodos;
conectivities = element;

end


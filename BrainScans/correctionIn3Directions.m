function correctionIn3Directions()

resultFile = 'materials_104by123by104.bin';

%% Read the coarse file:
disp('Opening raw data file...');
f = fopen('materials_100by119by100.bin','r');
BIG = fread(f,'float');
fclose(f);

%% Resize the data to a matrix:
disp('Resize raw data inside a matrix...');
Nx = 100;
Ny = 119;
Nz = 100;
for k = 1 : Nz
    for j = 1 : Ny
        for i = 1 : Nx
            rawData(i,j,k) = BIG(...
                i-1 + 100*((j-1)+119*(k-1))+1);
        end
    end
end

%% Ajouter de l'air (4 couches selon chaque direction):
disp('Add air layers...');
rawDataNew = zeros(size(rawData)+4);
for k = 1 : size(rawDataNew,3)
    for j = 1 : size(rawDataNew,2)
        for i = 1 : size(rawDataNew,1)
            
            if i > 2 && i < size(rawDataNew,1)-1
                if j > 2 && j < size(rawDataNew,2)-1
                    if k > 2 && k < size(rawDataNew,3)-1
                        rawDataNew(i,j,k) = rawData(i-2,j-2,k-2);
                    end
                end
            end
            
        end
    end
end
%% Save the raw data with air layers:
rawDataWithAirLayers = rawDataNew;
save('rawDataWithAirLayers','rawDataWithAirLayers');
load('rawDataWithAirLayers');
rawData = rawDataWithAirLayers;

%% Correction in the 3 directions (automatic filling of holes):
for dir = 1 : 3
    disp(['Apply automatic filling in direction ' num2str(dir)]);
    % We must do the correction for i or j or k kept constant !
    % -> We thus have three directions.
    [rawData] = correctionDir(rawData,dir);
end
%% Save the raw data corrected with auto-correction:
rawDataAutomaticCorrection = rawData;
save('rawDataAutomaticCorrection','rawDataAutomaticCorrection');
load('rawDataAutomaticCorrection');
rawData = rawDataAutomaticCorrection;
%% Apply a correction consisting in removing lines of material, or
%% removing alone material nodes or
%% removing alone air nodes.
%% Apply this correction until no node is corrected by the algorithm.
counter = 1;
phase = 1;
while counter ~= 0
    disp(['Correction of lines/alone number ' num2str(phase)]);
    [rawData,counter] = applyCorrectionLinesAlone(rawData);
    disp(['   > ' num2str(counter) ' nodes were modified...']);
    phase = phase + 1;
end
%% Do the same as in the C++ code to remove annoying nodes:
counter = 1;
phase = 1;
while counter ~= 0
    [rawData,counter] = applyCommeCPP(rawData);
    disp(['Applying correction as in C++ code number ' num2str(phase)]);
    disp(['    > ' num2str(counter) ' nodes were modified...']);
    phase = phase + 1;
end

%% Save the data:
rawDataNew_Modifier = rawData;
save('rawDataNew_Modifier','rawDataNew_Modifier');

%% Display in the 3 directions:
% for dir = 1 : 3
%     close all;
%     figure;
%     clear I;
%     for slice = 1 : size(rawData,dir)
%         if dir == 1
%             I(:,:) = rawData(slice,:,:);
%         elseif dir == 2
%             I(:,:) = rawData(:,slice,:);
%         else
%             I(:,:) = rawData(:,:,slice);
%         end
%         imagesc(I);
%         colormap;
%         drawnow;
%     end
% end

%% Write results to a file:
f = fopen(resultFile,'wb');
disp(['Writing results to ' resultFile]);
for k = 1 : size(rawData,3)
    for j = 1 : size(rawData,2)
        for i = 1 : size(rawData,1)
            BIG(i-1 + size(rawData,1)*((j-1)++size(rawData,2)*(k-1))+1) ...
                = rawData(i,j,k);
        end
    end
end

fwrite(f,BIG,'float');

fclose(f);

end

%% Apply automatic filling of holes:
function [A] = correctionDir(A,dir)
% Loop over slices in the given direction:
for k = 1 : size(A,dir)
    if dir == 1
        I(:,:) = A(k,:,:);
    elseif dir == 2
        I(:,:) = A(:,k,:);
    else
        I(:,:) = A(:,:,k);
    end
    % Open the image for groups of less than 50 pixels:
    BW = bwareaopen(I, 50);
    
    %% Inversion de l'image noir et blanc:
    BW2 = BW;%im2bw(BW);%#ok
    indNoir  = find(BW2 == 1);
    indBlanc = find(BW2 == 0);
    BW2(indNoir)  = 0;
    BW2(indBlanc) = 1;
    BW2 = bwareaopen(BW2,10);
    indNoir  = find(BW2 == 1);
    indBlanc = find(BW2 == 0);
    BW2(indNoir)  = 0;
    BW2(indBlanc) = 1;
    
    %% Creation du mask pour regionfill:
    Mask = abs(BW2-BW);
    
    %% Keep only bwareaopen de I:
    for i = 1 : size(I,1)
        for j = 1 : size(I,2)
            if BW(i,j) == 1
                I(i,j) = I(i,j);
            else
                I(i,j) = 0;
            end
        end
    end
    %% Interpolation
    I = regionfill(I,Mask);
    
    %% Faire round là où on a interpolé:
    for j = 1 : size(I,2)
        for i = 1 : size(I,1)
            if Mask(i,j) == 1
                I(i,j) = round(I(i,j));
            end
        end
    end
    if dir == 1
        A(k,:,:) = I;
    elseif dir == 2
        A(:,k,:) = I;
    else
        A(:,:,k) = I;
    end
end

end

function [rawData,counter_mod] = applyCorrectionLinesAlone(rawData)
%% Si un materiau a un seul materiau voisin, on supprime.
%% Si noeud air a que des voisins materiau, on remplit.
counter_mod = 0;
for k = 2 : size(rawData,3)-1
    for j = 2 : size(rawData,2)-1
        for i = 2 : size(rawData,1)-1
            
            if rawData(i,j,k) ~= 0 % materiau
                voisins = compterVoisin(rawData,i,j,k);
                if voisins == 1 || voisins == 0
                    rawData(i,j,k) = 0;
                    counter_mod = counter_mod + 1;
                end
            end
            
            if rawData(i,j,k) == 0 % air
                voisins = compterVoisin(rawData,i,j,k);
                if voisins == 6
                    rawData(i,j,k) = rawData(i,j,k-1);
                    counter_mod = counter_mod + 1;
                end
            end
            
        end
    end
end
disp(['Number of nodes modified ' num2str(counter_mod)]);
end

function voisins = compterVoisin(rawData,i,j,k)
voisins = 0;
if rawData(i+1,j,k) ~= 0
    voisins = voisins + 1;
end
if rawData(i-1,j,k) ~= 0
    voisins = voisins + 1;
end
if rawData(i,j+1,k) ~= 0
    voisins = voisins + 1;
end
if rawData(i,j-1,k) ~= 0
    voisins = voisins + 1;
end
if rawData(i,j,k+1) ~= 0
    voisins = voisins + 1;
end
if rawData(i,j,k-1) ~= 0
    voisins = voisins + 1;
end
end


function [rawData,counterError] = applyCommeCPP(rawData)

nodes_near_brain = zeros(size(rawData));
%% Detecter les noeuds air à la frontière:
for k = 2 : size(rawData,3)-1
    for j = 2 : size(rawData,2)-1
        for i = 2 : size(rawData,1)-1
            
            if rawData(i+1,j,k) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
            if rawData(i-1,j,k) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
            if rawData(i,j+1,k) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
            if rawData(i,j-1,k) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
            if rawData(i,j,k-1) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
            if rawData(i,j,k+1) ~= 0
                nodes_near_brain(i,j,k) = 1;
            end
        end
    end
end
%% Imposer convection:
counterError = 0;
for k = 2 : size(rawData,3)-1
    for j = 2 : size(rawData,2)-1
        for i = 2 : size(rawData,1)-1
            
            if nodes_near_brain(i,j,k) == 1
                
                a = 0;
                
                if rawData(i+1,j,k) ~= 0 && rawData(i+2,j,k) ~= 0
                    a = 1;
                end
                if rawData(i-1,j,k) ~= 0 && rawData(i-2,j,k) ~= 0 && a == 0
                    a = 1;
                end
                if rawData(i,j+1,k) ~= 0 && rawData(i,j+2,k) ~= 0
                    a = 1;
                end
                if rawData(i,j-1,k) ~= 0 && rawData(i,j-2,k) ~= 0 && a == 0
                    a = 1;
                end
                if rawData(i,j,k+1) ~= 0 && rawData(i,j,k+2) ~= 0
                    a = 1;
                end
                if rawData(i,j,k-1) ~= 0 && rawData(i,j,k-2) ~= 0 && a == 0
                    a = 1;
                end
                
                if a == 0
                    %error('Pas possible convection');
                    counterError = counterError + 1;
                    if rawData(i+1,j,k) ~= 0
                        rawData(i+1,j,k) = 0;
                        
                    elseif rawData(i-1,j,k) ~= 0
                        rawData(i-1,j,k) = 0;
                        
                    elseif rawData(i,j+1,k) ~= 0
                        rawData(i,j+1,k) = 0;
                        
                    elseif rawData(i,j-1,k) ~= 0
                        rawData(i,j-1,k) = 0;
                        
                    elseif rawData(i,j,k+1) ~= 0 
                        rawData(i,j,k+1) = 0;
                        
                    elseif rawData(i,j,k-1) ~= 0
                        rawData(i,j,k-1) = 0;
                    end
                end
                
                
            end
            
        end
    end
end
end

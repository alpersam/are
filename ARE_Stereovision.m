% Charge les images stéréo.
left = imread('left_450.jpg');
right = imread('right_450.jpg');

% ====================================
%        "Basic Block Matching"
% ====================================
% Calculer la disparité à l'aide d'une correspondance de bloc 
% de base avec une estimation de sous-pixel.

fprintf('Performing basic block matching...\n');

% On lance le minuteur
tic();

% Convertir les images de RVB grayscale 
% en faisant la moyenne des trois canaux de couleur.
leftI = mean(left, 3);
rightI = mean(right, 3);

% DbasicSubpixel contiendra le résultat de la correspondance des blocs. 
% Les valeurs seront en virgule flottante "simple précision" (32 bits).
DbasicSubpixel = zeros(size(leftI), 'single');

% Le "disparityRange" définit le nombre de pixels éloignés de l'emplacement du bloc
% dans la première image pour rechercher un bloc "correspondant" dans l'autre image.
disparityRange = 150;

% Définir la taille des blocs pour la correspondance des blocs.
halfBlockSize = 8;
blockSize = 2 * halfBlockSize + 1;

% Obtenir les dimensions de l'image
[imgHeight, imgWidth] = size(leftI);

% For each row 'm' of pixels in the image...
% On touche pas ici
for (m = 1 : imgHeight)
    	
	% Set min/max row bounds for the template and blocks.
	% e.g., for the first row, minr = 1 and maxr = 4
    minr = max(1, m - halfBlockSize);
    maxr = min(imgHeight, m + halfBlockSize);
	
    % For each column 'n' of pixels in the image...
    for (n = 1 : imgWidth)
        
		% Set the min/max column bounds for the template.
		% e.g., for the first column, minc = 1 and maxc = 4
		minc = max(1, n - halfBlockSize);
        maxc = min(imgWidth, n + halfBlockSize);
        
		% Define the search boundaries as offsets from the template location.
		% Limit the search so that we don't go outside of the image. 
		% 'mind' is the the maximum number of pixels we can search to the left.
		% 'maxd' is the maximum number of pixels we can search to the right.
		%
		% In the "Cones" dataset, we only need to search to the right, so mind
		% is 0.
		%
		% For other images which require searching in both directions, set mind
		% as follows:
        %   mind = max(-disparityRange, 1 - minc);
		mind = 0;
        maxd = min(disparityRange, imgWidth - maxc);

		% Select the block from the right image to use as the template.
        template = rightI(minr:maxr, minc:maxc);
		
		% Get the number of blocks in this search.
		numBlocks = maxd - mind + 1;
		
		% Create a vector to hold the block differences.
		blockDiffs = zeros(numBlocks, 1);
		
		% Calculate the difference between the template and each of the blocks.
		for (i = mind : maxd)
		
			% Select the block from the left image at the distance 'i'.
			block = leftI(minr:maxr, (minc + i):(maxc + i));
		
			% Compute the 1-based index of this block into the 'blockDiffs' vector.
			blockIndex = i - mind + 1;
		
			% Take the sum of absolute differences (SAD) between the template
			% and the block and store the resulting value.
			blockDiffs(blockIndex, 1) = sum(sum(abs(template - block)));
		end
		
		% Sort the SAD values to find the closest match (smallest difference).
		% Discard the sorted vector (the "~" notation), we just want the list
		% of indices.
		[temp, sortedIndeces] = sort(blockDiffs);
		
		% Get the 1-based index of the closest-matching block.
		bestMatchIndex = sortedIndeces(1, 1);
		
		% Convert the 1-based index of this block back into an offset.
		% This is the final disparity value produced by basic block matching.
		d = bestMatchIndex + mind - 1;
			
		% Calculate a sub-pixel estimate of the disparity by interpolating.
		% Sub-pixel estimation requires a block to the left and right, so we 
		% skip it if the best matching block is at either edge of the search
		% window.
		if ((bestMatchIndex == 1) || (bestMatchIndex == numBlocks))
			% Skip sub-pixel estimation and store the initial disparity value.
			DbasicSubpixel(m, n) = d;
		else
			% Grab the SAD values at the closest matching block (C2) and it's 
			% immediate neighbors (C1 and C3).
			C1 = blockDiffs(bestMatchIndex - 1);
			C2 = blockDiffs(bestMatchIndex);
			C3 = blockDiffs(bestMatchIndex + 1);
			
			% Adjust the disparity by some fraction.
			% We're estimating the subpixel location of the true best match.
			DbasicSubpixel(m, n) = d - (0.5 * (C3 - C1) / (C1 - (2*C2) + C3));
		end
    end

	% Update progress every 10th row.
	if (mod(m, 10) == 0)
		fprintf('  Image row %d / %d (%.0f%%)\n', m, imgHeight, (m / imgHeight) * 100);
	end
		
end

% Display compute time.
elapsed = toc();
fprintf('Calculating disparity map took %.2f min.\n', elapsed / 60.0);

% =========================================
%        Visualiser le Disparity Map
% =========================================

fprintf('Displaying disparity map...\n');


% Passer à la figure 1.
figure(1);

% Effacer la fenêtre de la figure actuelle.
clf;

% Afficher la carte des disparités. 
image(DbasicSubpixel);
axis image;

% Pour le color map.
colormap jet
clim([0 50]);

% Afficher la légende de la carte des couleurs.
colorbar;

% Le titre de la fenêtre.
title(strcat('Basic block matching, Sub-px acc., Search right, Block size = ', num2str(blockSize)));

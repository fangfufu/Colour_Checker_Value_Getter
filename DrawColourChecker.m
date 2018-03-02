function [ ] = DrawColourChecker(inMat, ncol, nrow)
%% DrawColourChart Given XYZ values, draw a colour chart
%   Mandatory Parameters : 
%       triplets : The triplets that describe a colour checker - could be 
%                   RGB or XYZ
%       ncols : The number of columns in this colour checker
%       nrows : The nubmer of rows in this colour checker
%


%% Sanity checks
if size(inMat, 1) ~= ncol * nrow
    error('Invalid matrix dimension.');
end

%% Convert XYZ to sRGB
% Note that we needed the whitepoint reading off the tiles
wp = GetWpFromColourChecker(inMat);
rgb = xyz2rgb(inMat, 'WhitePoint', wp);
rgb(rgb < 0) = 0;
rgb(rgb > 1) = 1;


%% Draw the rectangles
axis equal;
axis off;
k = 1;
for j = 1:nrow
    for i = 1:ncol
        rectangle('Position', [i * 10, -j * 10, 10, 10], 'FaceColor', ...
            rgb(k,:));
        k = k + 1;
    end
end

end

function [ wp, ind ] = GetWpFromColourChecker( mat )
%% GET_WP_FROM_COLOUR_CHART_MAT 
% Get the whitepoint from a matrix containing colour checker data.

matsum = sum(mat, 2);
ind = find(matsum == max(matsum));
% Make sure that we only output a unique whitepoint
ind = ind(1);
wp = mat(ind,:);
end



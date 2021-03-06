function [ ] = DrawColourChecker(inMat, ncol, nrow)
%% DrawColourChart Given XYZ values, draw a colour chart
%   Mandatory Parameters : 
%       triplets : The triplets that describe a colour checker - could be 
%                   RGB or XYZ
%       ncols : The number of columns in this colour checker
%       nrows : The nubmer of rows in this colour checker

%% Sanity checks
if size(inMat, 1) ~= ncol * nrow
    error('Invalid matrix dimension.');
end

%% Transpose colour checker 
% These functions were designed to be row-then-column, but everyone else's
% data seem to be column-then-row.
inMat = TransposeColourChecker(inMat, ncol, nrow);

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

function [ outImg, inPts ] = RegulariseImg(inImg, varargin)
%% RegulariseImg Transform an image so the camera is above the image.
%   Transform a colour checker image in such a way so that the camera is 90
%   degrees above the plane being imaged. This function is designed to work
%   with colour checker boards.
%
%   Example:
%       outImg = UntransformImg(inImg);
%       The in_img will be displayed. The user then has to select 4 corners
%       specified by the example image.
%
%   Mandatory parameters:
%       inImg : The image of the colour checker
%
%   Optional name-value pair parameters:
%       inPts : The coordinate of the corners of the colour checker.
%       refImg : Optional reference colour checker image.
%       refPts : The coordinate of the corners of the reference colour
%       checker
%       imgGamma : The gamma of the image. 
%
%   Output:
%       outImg : The image of the colour checker board after
%       transformation.
%       inPts : The coordinate of the corners of the colour checker.
%

%% Parse the input parser
p = inputParser;
addOptional(p, 'inPts', [], @(x) isnumeric(x) && numel(x) == 8);
addOptional(p, 'refImg', [], @(x) isnumeric(x) && ndims(x) == 3);
addOptional(p, 'refPts', [], @(x) isnumeric(x) && numel(x) == 8);
addParameter(p, 'imgGamma', 1/2.2, @(x) isscalar(x));

parse(p, varargin{:});

inPts = p.Results.inPts;
refImg = p.Results.refImg;
refPts = p.Results.refPts;
imgGamma = p.Results.imgGamma;

%% Constants

% These are the hardcoded reference coordinates taken from a random colour
% checker image off the Internet.
REFPTS = [ ...
       94.521        70.42; ...
       1021.3       71.334; ...
       1021.3       731.24; ...
       93.607       731.24; ...
       146.62       124.35; ...
       969.21       126.17; ...
       969.21       677.31; ...
       149.36        676.4; ...
        ];

% We want to re-write the coordinate system so it is in terms of the length
% of the colour checker
REFLEN = mean([sqrt(sum((REFPTS(1,:) - REFPTS(2,:)).^2)), ...
                sqrt(sum((REFPTS(3,:) - REFPTS(4,:)).^2))]);
REFPTS = REFPTS./REFLEN;

%% Checking what the user has supplied
if ~isempty(refImg)
    if ~isempty(refPts)
        % The user has supplied a non-empty refImg, a non-empty
        % ref_match_pts, we are drawing the matched points on the
        % refImg
        figure;
        imshow(refImg);
        hold on;
        plot(refPts(:,1), refPts(:,2), ...
            'rx', 'MarkerSize', 15, 'LineWidth', 2);
        hold off;
    else
        % refPts exist, but empty
        refPts = getImgCoordWrapper(refImg);
    end
end

% If refPts are still unset, these are the default ratio for the
% xrite colour chart
if isempty(refPts)
    refPts = REFPTS;
end

%% Process the input image
if isempty(inPts)
    % The user hasn't supplied matching points (this should normally be the
    % case)
    imshow(inImg.^imgGamma);
    waitfor(msgbox([...
        'Please click on the four corners of the colour checker in ', ...
        ' following order: top-left, top-right, bottom-right, ', ...
        'bottom-left. ']));
    [inPts, inLen] = GetCoordFromImg(4);
else
    inLen = mean([sqrt(sum((inPts(1,:) - inPts(2,:)).^2)), ...
                sqrt(sum((inPts(3,:) - inPts(4,:)).^2))]);
end

%% Calculate the geometric transform between two sets of points
% Scale up the refPts
refPts = refPts * inLen;
tform_in_to_ref = estimateGeometricTransform(inPts, refPts(1:4,:), ...
    'projective');
outImg = imwarp(inImg, tform_in_to_ref);
end

function [refPts, refLen] = getImgCoordWrapper(refImg)
% This function is called when the user has supplied a reference colour 
% checker, but without the coordinates for the corners. 
imshow(refImg);
waitfor(msgbox([...
        'Please click on the four corners of the colour checker in ', ...
        'following order: top-left, top-right, bottom-right, ', ...
        'bottom-left. ']));
[refPts, refLen] = GetCoordFromImg(4);
refPts = refPts ./ refLen;
end

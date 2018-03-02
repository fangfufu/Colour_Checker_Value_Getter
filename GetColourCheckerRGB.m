function [RGB, CCP] = GetColourCheckerRGB(varargin)
%GetColourCheckerRGB Get the RGB values of a colour checker
%   Note that this function outputs the colour checker in row-then-column 
%   format.
%   
%   Parameters:
%       Mandatory:
%       CI : Colour checker Image
%       Optional:
%       DI : Dark current Image
%       SI : Shading Image
%       Pairwise:
%       CCP : Colour checker parameter - to run the whole process again 
%       without manual selection again
%
%   Returns:
%       RGB : Colour checker RGB values
%       CCP : Colour checker parameter
%

%% Input parsers
p = inputParser;
addRequired(p, 'CI', @(x) size(x, 3) == 3);
addOptional(p, 'DI', [], @(x) size(x,3) == 3);
addOptional(p, 'SI', [], @(x) size(x,3) == 3);
addParameter(p, 'CCP', [], @(x) isstruct(x));

parse(p, varargin{:});
CI = p.Results.CI;
DI = p.Results.DI;
SI = p.Results.SI;
CCP = p.Results.CCP;

%% Pre-process the input image
% Remove dark current, if dark current image exists
if ~isempty(DI)
    CI = CI - DI;
    CI = im2double(CI);
else
    warning('GetColourCheckerValues:NoDarkCurrentImage', ...
        'Dark current image was not supplied');
    CI = im2double(CI);
end

if ~isempty(SI)
    disp('Shading image provided, dividing out the shading image');
    SI = im2double(SI);
    CI = CI ./ SI;
end

%% Extract the RGB values from the colour checker.
% If Colour Checker Parameter is not supplied
if isempty(CCP)
    % Set the gamma
    imgGamma = adjustGamma(CI);

    % We have to do the crop in 2 steps, because we can't really see the 
    % image without gamma correction, and we don't want gamma corrected
    % image as our output. 

    waitfor(msgbox(...
        'Let''s crop the region that contains the colour checker'));
    % CCBB : Colour Checker Bounding Box
    [~, CCBB] = imcrop(CI.^imgGamma);
    CI = imcrop(CI, CCBB);

    waitfor(msgbox(...
        'Let''s regularise colour checker image'));
    % CCCC : Colour Checker Corner Coordinates
    [CI, CCC] = RegulariseImg(CI, 'imgGamma', imgGamma);

    waitfor(msgbox(...
        'Let''s trim the regularised colour checker image'));
    % CCTBB : Colour Checker Trimming Bounding Box
    [~, CCTBB] = imcrop(CI.^imgGamma);
    CI = imcrop(CI, CCTBB);

    waitfor(msgbox(...
        'Finally, let''s read the values off the colour checker.'));
    % CCSS
    [RGB, CCSS] = ReadRGB(CI, imgGamma);

    % Save the Colour Checker Paramters
    CCP.CCBB = CCBB;
    CCP.CCC = CCC;
    CCP.CCTBB = CCTBB;
    CCP.CCSS = CCSS;
    CCP.imgGamma = imgGamma;
    
% If the Colour Checker Parameter is supplied
else
    % Unpack the parameters
    CCBB = CCP.CCBB;
    CCC = CCP.CCC;
    CCTBB = CCP.CCTBB;
    CCSS = CCP.CCSS;
    imgGamma = CCP.imgGamma;
    
    CI = imcrop(CI, CCBB);
    CI = RegulariseImg(CI, CCC, 'imgGamma', imgGamma);
    CI = imcrop(CI, CCTBB);
    RGB = ReadRGB(CI, imgGamma, CCSS);
end

% These functions were designed to be row-then-column, but everyone else's
% data seem to be column-then-row.
RGB = TransposeColourChecker(RGB, CCSS.yCount, CCSS.xCount);
end

function [ outImg, inPts ] = RegulariseImg(inImg, varargin)
%% UntransformImg Transform an image so the camera is above the image.
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

function [ pts, len ] = GetCoordFromImg(npts)
%GetCoordFromImg Get coordinates from an image.
%   Parameters
%       img : The input image
if ~exist('npts', 'var')
    npts = 4;
end

hold on;
pts = zeros(npts, 2);
for i = 1:npts
    this_pt = ginput(1);
    plot(this_pt(1), this_pt(2), 'rx', 'MarkerSize', 15, 'LineWidth', 2);
    pts(i, :) = this_pt;
end
hold off;

len = sqrt(sum((pts(1,:) - pts(2,:)).^2));

end

function [ imgGamma, figH ] = adjustGamma(img)
%% IMSHOWGAMMAADJUST Set a desired gamma value for an image

%% Constant definitions
% Gamma correction value for images - unlikely to be changed.
IMGGAMMA = 1/2.2;

img = im2double(img);
while(true)
    str = inputdlg('Please input a gamma value:',...
                 'Gamma', 1, {num2str(round(IMGGAMMA, 2))});
    if isempty(str)
        error('imshowGammaAdjust:Cancelled', ...
            'Cancelled by user.');
    end
    imgGamma = str2double(str);
    figH = figure;
    imshow(img.^imgGamma);
    choice = questdlg('Are you happy with the current gamma value?', ...
	'Gamma Correction', ...
	'Yes', 'No', 'Yes');
    if strcmp(choice, 'Yes')
        break;
    else
        close(figH);
    end
end

end

function [rgb, CCSStruct] = ReadRGB(img, imgGamma, CCSStruct)
%% GetRGB Get the RGB values off a colour checker
%   Parameters:
%       img - the raw image of a colour checker
%
%   Outputs:
%       CCSStruct : Colour Checker Structure Struct - a struct containing
%       the information about the colour checker
%       rgb : A matrix with the rgb values of the colour checker

%% Sanity check
if ndims(img) ~= 3
    error('The parameter img is not an image.');
end

if exist('CCSStruct', 'var')
    % The user supplied a colour checker settings struct
    if ~isa(CCSStruct, 'struct')
        % The supplied struct has a strong format.
        error('CCSStruct has a wrong format');
    end

    % We haven't errored out, so let's copy out CCSStruct's settings.
    xStart = CCSStruct.xStart;
    xSpacing = CCSStruct.xSpacing;
    xEnd = CCSStruct.xEnd;

    yStart = CCSStruct.yStart;
    ySpacing = CCSStruct.ySpacing;
    yEnd = CCSStruct.yEnd;

    imgGamma = CCSStruct.imgGamma;
    ps = CCSStruct.ps;

    xCount = CCSStruct.xCount;
    yCount = CCSStruct.yCount;

    cpc = CCSStruct.cpc;

    imshow(img.^imgGamma);

    hold on;
    for i = yStart:ySpacing:yEnd
        for j = xStart:xSpacing:xEnd
            rectangle('Position', [j-ps/2, i-ps/2, ps, ps], ...
                'EdgeColor', 'r', 'LineWidth', 2);
        end
    end
    plot(cpc(:,1), cpc(:,2), 'rx', 'MarkerSize', 15, 'LineWidth', 2);
    hold off;

else
    %% Gather some user input on how big the colour checker is.
    prompt = {['Please enter the number of rows and columns for your', ...
        'colour checker:'],...
        'Column:','Row:'};
    dlgTitle = 'Input';
    nLine = [0, 1, 1];
    defaultans = {'0','12','8'};
    answer = inputdlg(prompt,dlgTitle,nLine,defaultans);

    xCount = str2double(answer{2});
    yCount = str2double(answer{3});

    %% Start of colour patch selection loop
    while true
        %% Collecting coordinates of the corners of the colour checker.
        figH = figure;
        imshow(img.^imgGamma);

        waitfor(msgbox([...
        'Please click on the centre of the four corner patches of the ', ...
        'colour checker in the following order: top-left, top-right, ', ...
        'bottom-right, bottom-left. ']));

        % Note that the coordinates are in x,y format, even though the
        % matrices are addressed in y,x format!!!
        cpc = GetCoordFromImg(4); % Corner Patches Centre

        % Calculate the size of the colour chart
        xSize = mean([cpc(2,1) - cpc(1,1), cpc(3,1) - cpc(4,1)]);
        ySize = mean([cpc(3,2) - cpc(2,2), cpc(4,2) - cpc(1,2)]);

        % Calculate the size of spacing
        xSpacing = xSize / (xCount-1);
        ySpacing = ySize / (yCount-1);

        % Define the starting and ending positions
        xStart = cpc(1,1);
        yStart = cpc(1,2);
        xEnd = cpc(3,1);
        yEnd = cpc(4,2);

        %% Gather information on how big each colour patch is.
        waitfor(msgbox(['Please choose the sampling region by drawing', ...
            'a rectangle around the top left colour patch']));
        pRect = getrect();
        ps = mean([pRect(3), pRect(4)]);

        %% Looping through the colour checker to draw the rectangles
        hold on;

        % Make sure that we draw enough rectangles
        if size(xStart:xSpacing:xEnd) < xCount
            xEnd = xEnd + xSpacing;
        end

        if size(yStart:ySpacing:yEnd) < yCount
            yEnd = yEnd + ySpacing;
        end

        % The actual drawing loop
        for i = yStart:ySpacing:yEnd
            for j = xStart:xSpacing:xEnd
                rectangle('Position', [j-ps/2, i-ps/2, ps, ps], ...
                    'EdgeColor', 'r','LineWidth', 2);
            end
        end
        hold off;
        choice = questdlg('Are you happy with the selection?', ...
            'Selection confirmation', ...
            'Yes', 'No', 'Yes');
        if strcmp(choice, 'Yes')
            break;
        else
            close(figH);
        end
    end

    % Write the CCSStruct
    CCSStruct.xStart = xStart;
    CCSStruct.xSpacing = xSpacing;
    CCSStruct.xEnd = xEnd;

    CCSStruct.yStart = yStart;
    CCSStruct.ySpacing = ySpacing;
    CCSStruct.yEnd = yEnd;

    CCSStruct.imgGamma = imgGamma;
    CCSStruct.ps = ps;

    CCSStruct.xCount = xCount;
    CCSStruct.yCount = yCount;

    CCSStruct.cpc = cpc;

    %% End of Selection loop
end

%% Collecting colour patch samples
k = 1;
rgb = zeros(xCount * yCount, 3);
for i = yStart:ySpacing:yEnd
    for j = xStart:xSpacing:xEnd
        this_patch = img(floor(i-ps/2):ceil(i+ps/2), ...,
            floor(j-ps/2):ceil(j+ps/2), :);
        patch_mean = mean(mean(this_patch,1),2);
        rgb(k,:) = patch_mean;
        k = k + 1;
    end
end

end





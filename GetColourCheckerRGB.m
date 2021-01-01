function [RGB, CCP, CI] = GetColourCheckerRGB(varargin)
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
%       CI : Cropped colour checker image

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
    % CCCC : Colour Checker Coordinates
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

end

function [ imgGamma, figH ] = adjustGamma(img)
%% IMSHOWGAMMAADJUST Set a desired gamma value for an image

%% Constant definitions
% Gamma correction value for images - unlikely to be changed.
IMGGAMMA = 1/2.2;

img = im2double(img);
while(true)
    str = inputdlg('Please input a gamma value:',...
                 'Gamma', 1, {num2str(round(IMGGAMMA, 3))});
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






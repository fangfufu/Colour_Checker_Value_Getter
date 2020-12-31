function [rgb, CCSStruct] = ReadRGB(img, imgGamma, CCSStruct)
%% GetRGB Get the RGB values off a colour checker
%   Parameters:
%       img - the raw image of a colour checker
%       imgGamma - the gamma correction value
%       CCSStruct - the configuration of the colour checker
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

% These functions were designed to be row-then-column, but everyone else's
% data seem to be column-then-row.
rgb = TransposeColourChecker(rgb, CCSStruct.yCount, CCSStruct.xCount);

end 

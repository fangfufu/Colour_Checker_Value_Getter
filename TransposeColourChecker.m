function [ out_mat ] = TransposeColourChecker( in_mat, ncol, nrow)
%TransposeColourChecker In-place transposition of a nx3 matrix containing colour
%checker values.
%   You can use this function to interchange the format of the nx3 matrix which
%   contains the readings of the colour checker. You can inter-convert between 
%   row-then-column format and column-then-row format. 

% Note the col / row convention here.
tmp = zeros(ncol, nrow, 3);

patch_count = size(in_mat, 1);

% Current row and column
row = 1;
col = 1;
for i = 1:patch_count
    % Actual transposition happens here, note the row/col convention here
    tmp(row, col, :) = in_mat(i, :);
    if ~mod(col, nrow)
        row = row + 1;
        col = 1;
    else
        col = col + 1;
    end
end

out_mat = zeros(size(in_mat));
k = 1;
for i = 1:nrow
    for j = 1:ncol
        out_mat(k, :) = tmp(j,i,:);
        k = k + 1;
    end
end

end

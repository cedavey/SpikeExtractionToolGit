% cols = getColourMatrix(n)
% Gets a colour matrix. Differs to getUniqueCols because this returns the
% same colour matrix each time you run, whereas getUniqueCols randomises
% its colours. Also, getUniqueColour gets a single colour each time. Get
% colour matrix returns a matrix of unique colours, the size of which is
% n x 3, where n is the number of colours required.
% Inputs:
%	n    - specifies the number of colours required in the matrix
% Outputs:
%	cols - n x 3 matrix of RGB colours
function cols = getColourMatrix(n)
    if nargin<1
        disp(['usage: cols = getColourMatrix(n)']);
        return;
    end
    cols = [0         0     1.0000; ...
             0    0.5000         0; ...
             0    0.7500    0.7500; ...
        0.2500    0.2500    0.2500; ...
        0.7500         0    0.7500; ...
        0.7500    0.7500         0; ...
        1.0000         0         0];

    while n > size(cols,1)
        newcols  = ( cols(2:end,:) + cols(1:end-1,:) ) / 2;
        ncols = length( cols ) + length( newcols );
        tmp   = zeros( ncols, 3 );
        tmp(1:2:end,:) = cols; 
        tmp(2:2:end,:) = newcols; 
        cols = tmp;
    end
    
    cols = cols(1:n,:);
end
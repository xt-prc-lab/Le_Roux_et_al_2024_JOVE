%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%         Function used to apply a 2D shift to a 2D stack. Subpixel accuracy is achieved through linear interpolation.         %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im [2D matrix]: 2D matrix to be shifted.                                                                               %
%       Xj [double]: sub-pixel shift in the x direction, or the columns of the 2D matrix.                                      %
%       Xi [double]: sub-pixel shift in the y direction, or the rows of the 2D matrix.                                         %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       im1 [2D matrix]: shifted 2D matrix.                                                                                    %
%                                                                                                                              %
%   Last Revison Date: 29/07/2024                                                                                              %
%   Based on the code created by Prof. Xavier Trepat.                                                                          %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function im1 = xsubpix_shift_2D(im, Xj, Xi)

    I = floor(-Xi);
    J = floor(-Xj);

    xi = -Xi-I;
    xj = -Xj-J;

    x0 = max(1, 2-I);
    y0 = max(1, 2-J);
    xf = min(size(im, 1), size(im, 1)-I-2);
    yf = min(size(im, 2), size(im, 2)-J-2);

    im1 = zeros(size(im));

    im1(x0:xf, y0:yf) = (1-xi)*(1-xj)*im((x0:xf)+I  , (y0:yf)+J  ) + ...
                                xi*xj*im((x0:xf)+I+1, (y0:yf)+J+1) + ...
                            xi*(1-xj)*im((x0:xf)+I+1, (y0:yf)+J  ) + ...
                            xj*(1-xi)*im((x0:xf)+I  , (y0:yf)+J+1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

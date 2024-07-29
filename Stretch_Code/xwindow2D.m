%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                               Function used to apply a multiplicative 2D window to a 2D matrix.                              %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im [2D matrix]: 2D matrix to be windowed.                                                                              %
%       method [int, optional]: type of window applied. 1 means Hanning window, otherwise means no windowing.                  %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       im1 [2D matrix]: windowed 2D matrix.                                                                                   %
%                                                                                                                              %
%   Last Revison Date: 29/07/2024                                                                                              %
%   Based on the code created by Prof. Xavier Trepat.                                                                          %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function im = xwindow2D(im,method)

    [n, m] = size(im);
    w = zeros(n, m);

    if method==1

        [xi,xj] = meshgrid(1:m,1:n);
        w = 0.25*(1-cos(2*pi*xi/m)).*(1-cos(2*pi*xj/n));

    elseif method==2

        dev=8;

        [wn,wm] = meshgrid(0:m-1,0:n-1);
        wn = (1-cos(pi*wn/(m-1)).^dev);
        wm = (1-cos(pi*wm/(n-1)).^dev);
        w = wn.*wm;

    end

    im = im.*w;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%               Function that applies a fast PIV to calculate the displacement field of im1 with respect to im2.               %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im1 [2D matrix]: reference 2D image.                                                                                   %
%       im2 [2D matrix]: deformed 2D image on which the displacement field will be computed.                                   %
%       iblocksize [int]: number of rows of the PIV window.                                                                    %
%       jblocksize [int]: number of columns of the PIV window.                                                                 %
%       overlap [scalar]: overlap, in rows and columns, of the PIV  windows.                                                   %
%       padding [int]: zero-padding, in rows and columns, around the PIV  windows.                                             %
%       window [int]: selects the multiplicative window applied on the PIV  windows.                                           %
%                     1 means Hanning window, otherwise means no windowing.                                                    %
%       method [int]: Method used for the location of the peak of the 2D correlation.                                          %
%                     1 means the 2D centroid, 2 means Gaussian fit.                                                           %
%       N [int]: number of points around the pixel peak location used to locate it with subpixel accuracy.                     %
%       Ni [int]: minimum number of iterations. The algorithm will iterate until the number of iterations is larger than Ni    %
%                 and the desired accuracy is reached.                                                                         %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       x [2D matrix]: x-coordinate of the displacement field.                                                                 %
%       y [2D matrix]: y-coordinate of the displacement field.                                                                 %
%       dx [2D matrix]: x-component of the displacement field.                                                                 %
%       dy [2D matrix]: y-component of the displacement field.                                                                 %
%       pkh [2D matrix]: value of the cross-correlation at the peak.                                                           %
%                                                                                                                              %
%   Last Revison Date: 29/07/2024                                                                                              %
%   Based on the 2D PIV created by Xavier Trepat (03/2008).                                                                    %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       "A correlation-based continuous window-shift technique to reduce the peak-locking effect in digital PIV image          %
%       evaluation"; Lichuan Gui and Steven T. Wereley; Experiments in Fluids, April 2002, Volume 32, Issue 4, pp 506-517.     %
%       https://doi.org/10.1007/s00348-001-0396-1                                                                              %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, dx, dy, pkh] = disp_on_blocks_fast_2D(im1, im2, iblocksize, jblocksize, overlap, padding, window, method, N, Ni)

    sz = size(im1);     % Size of the image.

    xaux = (sz(1)/iblocksize-1)*(1/(1-overlap)) + 1;
    yaux = (sz(2)/jblocksize-1)*(1/(1-overlap)) + 1;

    if xaux>2048 || xaux<0 || yaux>2048 || yaux<0

        x = NaN;
        y = NaN;
        dx = NaN;
        dy = NaN;
        pkh = NaN;

    else

        im1 = xExpandMatrix_3D(im1, 1, 1, 1, padding, padding, padding, padding, 0, 0, 0); 
        im2 = xExpandMatrix_3D(im2, 1, 1, 1, padding, padding, padding, padding, 0, 0, 0);

        inci = round(iblocksize*(1-overlap));           % Define increment in X.
        incj = round(jblocksize*(1-overlap));           % Define increment in Y.

        if (inci<1 || inci>iblocksize)
            error('Wrong Overlap in Correlation Algorithm')
        end

        if (incj<1 || incj>jblocksize)
            error('Wrong Overlap in Correlation Algorithm')
        end

        [x, y] = meshgrid( ...              % x and y coordinates of the  windows.
            (1:inci:sz(1)-iblocksize+1) + iblocksize/2, ...
            (1:incj:sz(2)-jblocksize+1) + jblocksize/2);

        dx  = zeros(size(x));               % Displacements in x.
        dy  = zeros(size(x));               % Displacements in y.
        pkh = zeros(size(x));               % Height of the xcorr peak.
        c   = zeros(size(x));               % Cross correlation of the images.

        % Loop along the rows.
        for ki = 1:inci:sz(1)-iblocksize+1

            fprintf('.');

            % Loop along the columns.
            for kj = 1:incj:sz(2)-jblocksize+1

                % Initialize iterative process.
                niterations = 0;
                DX = inf;
                DY = inf;

                while niterations<=Ni && (abs(dx((ki+inci-1)/inci, (kj+incj-1)/incj) - DX)>0.02 || ...
                                          abs(dy((ki+inci-1)/inci, (kj+incj-1)/incj) - DY)>0.02)    % Set iteration conditions.  

                    niterations = niterations+1;
                    DX = dx((ki+inci-1)/inci, (kj+incj-1)/incj);
                    DY = dy((ki+inci-1)/inci, (kj+incj-1)/incj);

                    im11 = im1(ki:ki+iblocksize+2*padding-1, kj:kj+jblocksize+2*padding-1);
                    im22 = im2(ki:ki+iblocksize+2*padding-1, kj:kj+jblocksize+2*padding-1);

                    if DX || DY         % Skip first iteration.

                        im11 = xsubpix_shift_2D(im11, DX/2, DY/2);
                        im22 = xsubpix_shift_2D(im22,-DX/2,-DY/2);

                    end

                    x_sect = (1:iblocksize)+padding;
                    y_sect = (1:jblocksize)+padding;

                    % Subtract the mean.
                    im11(x_sect, y_sect) = im11(x_sect, y_sect) - mean(mean(im11(x_sect, y_sect)));
                    im22(x_sect, y_sect) = im22(x_sect, y_sect) - mean(mean(im22(x_sect, y_sect)));            

                    if window==1

                        % Multiply each block by a window.
                        im11(x_sect, y_sect) = xwindow2D(im11(x_sect, y_sect), window);
                        im22(x_sect, y_sect) = xwindow2D(im22(x_sect, y_sect), window);

                    end

                    c = real(xcorrf3(im11(x_sect, y_sect), im22(x_sect, y_sect), 'no'));

                    [xc, yc, ~] = x_cntr_3D(c, N, method);
                    dx((ki+inci-1)/inci, (kj+incj-1)/incj) = real(DX+xc);
                    dy((ki+inci-1)/inci, (kj+incj-1)/incj) = real(DY+yc);

                end

                pkh((ki+inci-1)/inci, (kj+incj-1)/incj) = max(abs(c(:)));       % Store peak height.     
 
            end

        end

    end

    fprintf('\n');

end

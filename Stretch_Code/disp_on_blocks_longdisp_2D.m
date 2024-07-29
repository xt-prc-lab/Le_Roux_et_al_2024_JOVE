%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%        Function that applies an accurate (slower) PIV to calculate the displacement field of im1 with respect to im2.        %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       im1 [2D matrix]: reference 2D image.                                                                                   %
%       im2 [2D matrix]: deformed 2D image on which the displacement field will be computed.                                   %
%       iblocksize [int scalar]: number of rows of the PIV window.                                                             %
%       jblocksize [int scalar]: number of columns of the PIV window.                                                          %
%       overlap [double scalar]: overlap in rows and columns of the PIV windows.                                               %
%       padding [int scalar]: zero-padding, in rows and columns, around the PIV window.                                        %
%       blocksizedec [int scalar]: amount by which each side of the interrogation window, in xy, shrinks at each iteration.    %
%       window [int scalar]: selects the multiplicative window applied on the PIV windows.                                     %
%                            1 means Hanning window, otherwise means no windowing.                                             %
%       method [int scalar]: Method used for the location of the peak of the 2D correlation.                                   %
%                            1 means the 2D centroid, 2 means Gaussian fit.                                                    %
%       N [int scalar]: number of points around the pixel peak location used to locate it with subpixel accuracy.              %
%       N_lim [int scalar]: maximum number of iterations after decrease of window size. The algorithm will iterate until the   %
%                           number of iterations, after the initial window size decrease, is larger than N_lim or all windows  %
%                           have converged, whatever happens first.                                                            %
%       epsilon [double scalar]: maximum acceptable relative error of the iteration for each window to be considered as        %
%                                converged.                                                                                    %
%       N_conver [int scalar]: stop iterating after N_conver iterations without new windows converging.                        %
%       del [int scalar]: For each window where longdisp is required, apply longdisp to the neighbouring "del" windows around  %
%                         it.                                                                                                  %
%       dynamic_window [int scalar]: when determining the initial window size, we could use different criteria for the         %
%                                    different windows. 1 if the sizes of the PIV windows are determined dynamically; 2 if all %
%                                    the windows larger than blocksize have the maximum initial window size; 0 if all the      %
%                                    windows have the maximum window size.                                                     %
%       filt_size [int scalar]: XY-size of the Predictor filter, in window lengths, i.e. filter_size = filt_size*window_size.  %
%       filt_pow [double scalar]: power of the Predictor filter.                                                               %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       x [2D matrix]: x-coordinate of the displacement field.                                                                 %
%       y [2D matrix]: y-coordinate of the displacement field.                                                                 %
%       dx [2D matrix]: x-component of the displacement field.                                                                 %
%       dy [2D matrix]: y-component of the displacement field.                                                                 %
%       pkh [2D matrix]: value of the cross-correlation at the peak.                                                           %
%       conver [vector]: percentage of PIV windows that have reached convergence at each iteration.                            %
%       iter_time [vector]: time taken to run each iteration.                                                                  %
%       ws_0 [2D matrix]: initial window size, in xy, of each PIV window, before any shrinking.                                %
%       Convergence_Mask [2D matrix]: for each PIV window, 1 if it has reached convergence, 0 otherwise.                       %
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
%       "Particle Image Velocimetry, A Practical Guide"; Markus Raffel, Christian E. Willert, Steve T. Wereley,                %
%       Jürgen Kompenhans; Springer, 2007. https://doi.org/10.1007/978-3-540-72308-0                                           %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, y, dx, dy, pkh, conver, iter_time, ws0, Convergence_Mask] = ...
    disp_on_blocks_longdisp_2D(im1, im2, iblocksize, jblocksize, overlap, padding, blocksizedec, window, method, ...
                               N, N_lim, epsilon, N_conver, del, dynamic_window, filt_size, filt_pow)

    pad_filter = round(filt_size/(2*(1-overlap)));          % The filter used on the displacement will create border effects. 
                                                            % In order to avoid them, we create a rim around the matrix, with 
                                                            % size pad_filter, that is the specular image of the displacement
                                                            % field.

    % Size of the image. Keep in mind that the first index (rows) corresponds to y, the second index (columns) to x.
    sz = size(im1);
    
    method = 1;

    im1 = xExpandMatrix_3D(im1, 1, 1, 1, padding, padding, padding, padding, 0, 0, 0);
    im2 = xExpandMatrix_3D(im2, 1, 1, 1, padding, padding, padding, padding, 0, 0, 0);

    inci = round(iblocksize*(1-overlap));           % Define increment in Y.
    incj = round(jblocksize*(1-overlap));           % Define increment in X.
    
    if (inci<1 || inci>iblocksize)
        error('Wrong Overlap in Correlation Algorithm')
    end

    if (incj<1 || incj>jblocksize)
        error('Wrong Overlap in Correlation Algorithm')
    end

    % x and y coordinates of the windows.
    [x, y] = meshgrid((1:inci:sz(2)-iblocksize+1) + iblocksize/2, ...
                      (1:incj:sz(1)-jblocksize+1) + jblocksize/2);

    dx  = zeros(size(x));               % Displacements in x.
    dy  = zeros(size(x));               % Displacements in y.
    pkh = zeros(size(x));               % Height of the xcorr peak.
    % c   = zeros(size(x));               % Cross correlation of the images.

    ws  = ones(size(x))*iblocksize;     % Size of each window for the longdigsp calculations.

    dx_filt = zeros(size(x) + 2*pad_filter);            % Displacements in x, filtered.
    dy_filt = zeros(size(x) + 2*pad_filter);            % Displacements in y, filtered.

    dx_filt_2 = zeros(size(x) + 2*pad_filter);          % Displacements in x, filtered.
    dy_filt_2 = zeros(size(x) + 2*pad_filter);          % Displacements in y, filtered.

    Convergence_Mask = zeros(size(x));          % 1 for PIV windows that reach convergence, 0 otherwise.
    e = zeros(size(x));
    inc_err0 = zeros(size(x));

    niterations = 1;
    exit = 0;
    conver = [0];

    % Index of the non-converged windows.
    conv_vect = cell(1, size(Convergence_Mask, 1));

    for ki=1:inci:sz(1)-iblocksize+1
        conv_vect{(ki+inci-1)/inci} = 1:incj:sz(2)-jblocksize+1;
    end
    
    % Make a first PIV run with the smallest window size. It should be very fast. Use the displacement of each window to decide 
    % the window size needed to resolve the displacement of this window, i.e. the one-quarter rule: disp < window-size / 4.
    % See "Particle Image Velocimetry, a Practical Guide", by M. Raffel et al.

    % Coordinates of each point in the windows.
    x_sect = 1:iblocksize + 2*padding;
    y_sect = 1:jblocksize + 2*padding;

    time_iter0 = tic;

    for ki = 1:inci:sz(1)-iblocksize+1

        ii = (ki+inci-1)/inci;

        fprintf('*');

        for kj = 1:incj:sz(2)-jblocksize+1

            jj = (kj+incj-1)/incj;

            % Interrogation windows, cropped to the correct size for each itaration.
            im11 = im1(y_sect + ki-1, x_sect + kj-1);
            im22 = im2(y_sect + ki-1, x_sect + kj-1);

            % Normalize the images.
            if std_flat(im11) && std_flat(im22)

                im11 = (im11 - mean(im11, 'all'))/std_flat(im11);
                im22 = (im22 - mean(im22, 'all'))/std_flat(im22);

            else

                im11 = im11 - mean(im11, 'all');
                im22 = im22 - mean(im22, 'all');

            end

            % Multiply each block by a window.
            im11 = xwindow2D(im11, window);
            im22 = xwindow2D(im22, window);

            c = real(xcorrf3(im11, im22, 'no'));

            [xc, yc, ~] = x_cntr_3D(c, N, method);

            dx(ii, jj) = real(xc);          % Accumulated displacement across iterations.
            dy(ii, jj) = real(yc);

        end

    end

    iter_time(1) = toc(time_iter0);
    time_iter0 = iter_time(1)/60

    [n, m] = size(ws);

    for ii=1:n
        for jj=1:m

            % If the displacement of a box (or neighboutring boxes) is larger than one fourth of the user provided 
            % window size, it cannot be resolved with the user provided box size. In that case, enlarge the initial box size 
            % (longdisp algorithm).
            tmp = max(max(sqrt(dx(max(1, ii-del):min(n, ii+del), max(1, jj-del):min(m, jj+del)).^2 + ...
                               dy(max(1, ii-del):min(n, ii+del), max(1, jj-del):min(m, jj+del)).^2)));
            ws(ii, jj) = ws(ii, jj) + ceil(max(0, 4*tmp - ws(ii, jj))/(2*blocksizedec))*2*blocksizedec;

        end
    end

    if ~dynamic_window 
        ws(:) = max(ws(:));
    end

    ws0 = ws;
    
    pad_dec = max(padding, (max(ws(:)) - iblocksize)/2);            % Padding for the increasing of the edges in the longdisp.

    im1 = xExpandMatrix_3D(im1(1+padding:end-padding, 1+padding:end-padding), 1, 1, 1, pad_dec, pad_dec, ...
                           pad_dec, pad_dec, 0, 0, 0);
    im2 = xExpandMatrix_3D(im2(1+padding:end-padding, 1+padding:end-padding), 1, 1, 1, pad_dec, pad_dec, ...
                           pad_dec, pad_dec, 0, 0, 0);

    % Coordinates of each point in the largest windows.
    x_sect = 1:(iblocksize + 2*pad_dec);
    y_sect = 1:(jblocksize + 2*pad_dec);

    % Perform the iterations with box size decrease.
    while (exit~= 1)
 
        niterations
        time_iter = tic;

        for ki = 1:inci:sz(1)-iblocksize+1        

            ii = (ki+inci-1)/inci;

            fprintf('.');

            for kk = numel(conv_vect{ii}):-1:1

                kj = conv_vect{ii}(kk);
                jj = (kj+incj-1)/incj;

                % This is the subpixel amount of shifting until getting peak locking conditions.
            	DX = dx_filt_2(ii + pad_filter, jj + pad_filter);
                DY = dy_filt_2(ii + pad_filter, jj + pad_filter);

                x_wind_dec = (x_sect(end) - ws(ii, jj))/2.;             % Reduction of the cropped window.
                y_wind_dec = (y_sect(end) - ws(ii, jj))/2.;

                % Coordinates of the window where we will perform the cross-correlation. It is at the center of x_sect, y_sect.
                x_wind = x_sect(1+x_wind_dec:x_sect(end)-x_wind_dec);
                y_wind = y_sect(1+y_wind_dec:y_sect(end)-y_wind_dec);

                x_wind_shift = x_sect(1+max(0, x_wind_dec-1):x_sect(end)-max(0, x_wind_dec-1));
                y_wind_shift = y_sect(1+max(0, y_wind_dec-1):y_sect(end)-max(0, y_wind_dec-1));

                x_wind_shift = x_wind_shift((x_wind_shift-round( DX/2)>=1)&(x_wind_shift-round( DX/2)<=x_sect(end))& ...
                                            (x_wind_shift-round(-DX/2)>=1)&(x_wind_shift-round(-DX/2)<=x_sect(end)));
                y_wind_shift = y_wind_shift((y_wind_shift-round( DY/2)>=1)&(y_wind_shift-round( DY/2)<=y_sect(end))& ...
                                            (y_wind_shift-round(-DY/2)>=1)&(y_wind_shift-round(-DY/2)<=y_sect(end))) ;

                % Interrogation windows, cropped to the correct size for each iteration.
                im11 = im1(y_sect+ki-1, x_sect+kj-1);
                im22 = im2(y_sect+ki-1, x_sect+kj-1);

                % Error parameters of the first iteration (it is used after the second iteration).
                if niterations == 2

                	I0 = sum((im11 - im22).^2, 'all');

                    err0 = sqrt(I0/((x_sect(end)-x_sect(1)+1)*(y_sect(end)-y_sect(1)+1)))/...
                    	   ((std_flat(im11) + std_flat(im22))/2);

                end

                % Subpixel shift of the image.
                im11(y_wind_shift, x_wind_shift) = xsubpix_shift_2D(...
                            im11(y_wind_shift-round( DY/2), x_wind_shift-round( DX/2)),  DX/2-round( DX/2),  DY/2-round( DY/2));
                im22(y_wind_shift, x_wind_shift) = xsubpix_shift_2D(...
                            im22(y_wind_shift-round(-DY/2), x_wind_shift-round(-DX/2)), -DX/2-round(-DX/2), -DY/2-round(-DY/2));

                % Test the convergence of the previous step. If we test the convergence of the current step, we have to
                % perform again xsubpix_shift_3D_mex(  ), which is a very expensive function. This way, reduse the calls
                % to the function by half.
                I = sum((im11(y_wind, x_wind) - im22(y_wind, x_wind)).^2, 'all');

                err = sqrt(I/((x_wind(end)-x_wind(1)+1)*(y_wind(end)-y_wind(1)+1)))/...
                      ((std_flat(im11(y_wind, x_wind)) + std_flat(im22(y_wind, x_wind)))/2);

                if niterations == 2
                	inc_err0(ii, jj) = err0 - err;
                end

                if (abs((err - e(ii, jj))/inc_err0(ii, jj)) <= epsilon) && (ws(ii, jj) == iblocksize)

                	Convergence_Mask(ii, jj) = 1;
                    conv_vect{ii}(kk) = [];
                    dx(ii, jj) = dx_filt(ii + pad_filter, jj + pad_filter); 
                    dy(ii, jj) = dy_filt(ii + pad_filter, jj + pad_filter);

                else

                	% Normalize the images.
                    if std_flat(im11(y_wind, x_wind)) && std_flat(im22(y_wind, x_wind))

                        im11(y_wind, x_wind) = (im11(y_wind, x_wind) - mean(reshape(im11(y_wind, x_wind), [], 1)))/...
                                               std_flat(im11(y_wind, x_wind));
                        im22(y_wind, x_wind) = (im22(y_wind, x_wind) - mean(reshape(im22(y_wind, x_wind), [], 1)))/...
                                               std_flat(im11(y_wind, x_wind));

                    else

                        im11(y_wind, x_wind) = im11(y_wind, x_wind) - mean(reshape(im11(y_wind, x_wind), [], 1));
                        im22(y_wind, x_wind) = im22(y_wind, x_wind) - mean(reshape(im22(y_wind, x_wind), [], 1));

                    end

                    % Multiply each block by a window.
                    im11(y_wind, x_wind) = xwindow2D(im11(y_wind, x_wind), window);
                    im22(y_wind, x_wind) = xwindow2D(im22(y_wind, x_wind), window);

                    c = real(xcorrf3(im11(y_wind, x_wind), im22(y_wind, x_wind), 'no'));

                    [xc, yc, ~] = x_cntr_3D(c, N, method);

                    dx(ii, jj) = real(DX + xc);         % Accumulated displacement across iterations.
                    dy(ii, jj) = real(DY + yc);
                    pkh(ii, jj) = max(abs(c(:)));       % Store peak height.

                end

                e(ii, jj) = err;

                % Determine whether the window need to be further shrunk.
                ws(ii, jj) = max(iblocksize, ws(ii, jj) - 2*blocksizedec);

            end

        end 

        time_iter1 = toc(time_iter )/60

        % Add the specular image of the borders of the displacement field before filtering.
        dx_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter) = dx;
        dy_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter) = dy;

        dx_filt(1:pad_filter, :) = dx_filt(1+2*pad_filter:-1:2+pad_filter, :);
        dy_filt(1:pad_filter, :) = dy_filt(1+2*pad_filter:-1:2+pad_filter, :);

        dx_filt(end-pad_filter+1:end, :) = dx_filt(end-pad_filter-1:-1:end-2*pad_filter, :);
        dy_filt(end-pad_filter+1:end, :) = dy_filt(end-pad_filter-1:-1:end-2*pad_filter, :);

        dx_filt(:, 1:pad_filter) = dx_filt(:, 1+2*pad_filter:-1:2+pad_filter);
        dy_filt(:, 1:pad_filter) = dy_filt(:, 1+2*pad_filter:-1:2+pad_filter);

        dx_filt(:, end-pad_filter+1:end) = dx_filt(:, end-pad_filter-1:-1:end-2*pad_filter);
        dy_filt(:, end-pad_filter+1:end) = dy_filt(:, end-pad_filter-1:-1:end-2*pad_filter);

        % Filtering.
        dx_filt_2 = inpaint_nans(dx_filt, 3);
        dy_filt_2 = inpaint_nans(dy_filt, 3);

        dx_filt_2(dx_filt_2> max(ws0(:))/2) =  max(ws0(:))/2;
        dx_filt_2(dx_filt_2<-max(ws0(:))/2) = -max(ws0(:))/2;

        dy_filt_2(dy_filt_2> max(ws0(:))/2) =  max(ws0(:))/2;
        dy_filt_2(dy_filt_2<-max(ws0(:))/2) = -max(ws0(:))/2;

        % Outlier removal.
        dx_filt_2 = removeOutliers_2D(dx_filt_2);
        dy_filt_2 = removeOutliers_2D(dy_filt_2);

        dx_filt_2 = dx_filt_2{1};
        dy_filt_2 = dy_filt_2{1};

        dx_filt_2 = filter_2D_Disp(dx_filt_2, filt_size*[1/(1 - overlap), 1/(1 - overlap)], filt_pow);
        dy_filt_2 = filter_2D_Disp(dy_filt_2, filt_size*[1/(1 - overlap), 1/(1 - overlap)], filt_pow);

        % Deform images and check convergence and decide if the window needs to be shrunk or not and for which magnitude.
        fprintf('\n');

        if niterations >= 2

            conver(niterations) = sum(Convergence_Mask(:))/numel(Convergence_Mask)*100;
            disp(['At iteration number ', num2str(niterations-1), ',   ', num2str(conver(end)), ...
                  ' % of points have reached convergence'])

            exit = (niterations >= pad_dec/blocksizedec + N_lim) || (conver(end)==100);

            if niterations >= N_conver+1
                exit = exit || (sum(abs(diff(conver(min(numel(conver), end-N_conver+1):end))))==0);
            end

        end

        iter_time(niterations + 1) = toc(time_iter);       
               
        niterations = niterations + 1;

    end

    fprintf('\n');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

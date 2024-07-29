%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                         Function used to calculate the 2D strain field from a 2D displacement field.                         %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [structure]: structure containing the different parameters for the analysis.                                  %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Last Revison Date: 29/07/2024                                                                                              %
%   Based on codes created by Xavier Trepat's group.                                                                           %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function File = strain_finder_2D(Settings, File)

    disp(' ');
    disp('         *** Start strain finder ***       ');
    disp(' ');

    pad_filter = round(Settings.filt_size_strain/(2*(1-Settings.Overlap)));         % The filter used on the displacement will 
                                                                                    % create border effects. In order to avoid 
                                                                                    % them, we create a rim around the matrix, 
                                                                                    % with size pad_filter, that is the specular 
                                                                                    % image of the displacement field.

    frames = File.NFiles.Beads;

    % Loop over all the timepoints.
    for k=1:frames

        [~, temp, ext] = fileparts(File.Name(k).Beads);

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        temp = [temp, ext];
        temp = temp(1:end-length(Settings.Imgfmt)-1);

        tempFile = [File.resultsname filesep 'Displacements' filesep 'GelDisp_', temp, '.csv'];

        if exist(tempFile, 'file') == 2

            TempData = readmatrix(tempFile);

            x = reshape(TempData(:,1), File.TractionSize(k).i, File.TractionSize(k).j);
%             y = reshape(TempData(:,2), File.TractionSize(k).i, File.TractionSize(k).j);
            DispX = reshape(TempData(:,3), File.TractionSize(k).i, File.TractionSize(k).j);
            DispY = reshape(TempData(:,4), File.TractionSize(k).i, File.TractionSize(k).j);
%             c = reshape(TempData(:,5), File.TractionSize(k).i, File.TractionSize(k).j);

%             conv_th = 0.80;
            conv_th = Settings.MinCorrelationTrheshold;

            pkh = reshape(TempData(:,5), File.TractionSize(k).i, File.TractionSize(k).j);
            conv_mask = reshape(TempData(:,end), File.TractionSize(k).i, File.TractionSize(k).j);
            mask = (conv_mask) & (pkh >= conv_th);

            DispX(~mask) = NaN;
            DispY(~mask) = NaN;

            DispX_filt = zeros(size(x) + 2*pad_filter);             % Displacements in x, filtered.
            DispY_filt = zeros(size(x) + 2*pad_filter);             % Displacements in y, filtered.

            % Add the specular image of the borders of the displacement field before filtering.
            DispX_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter) = DispX;
            DispY_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter) = DispY;

            DispX_filt(1:pad_filter, :) = DispX_filt(1+2*pad_filter:-1:2+pad_filter, :);
            DispY_filt(1:pad_filter, :) = DispY_filt(1+2*pad_filter:-1:2+pad_filter, :);

            DispX_filt(end-pad_filter+1:end, :) = DispX_filt(end-pad_filter-1:-1:end-2*pad_filter, :);
            DispY_filt(end-pad_filter+1:end, :) = DispY_filt(end-pad_filter-1:-1:end-2*pad_filter, :);

            DispX_filt(:, 1:pad_filter) = DispX_filt(:, 1+2*pad_filter:-1:2+pad_filter);
            DispY_filt(:, 1:pad_filter) = DispY_filt(:, 1+2*pad_filter:-1:2+pad_filter);

            DispX_filt(:, end-pad_filter+1:end) = DispX_filt(:, end-pad_filter-1:-1:end-2*pad_filter);
            DispY_filt(:, end-pad_filter+1:end) = DispY_filt(:, end-pad_filter-1:-1:end-2*pad_filter);

            % Applying the filtering and outlier removal to the final solution.

            % Filtering.
            DispX_filt = inpaint_nans(DispX_filt, 3);
            DispY_filt = inpaint_nans(DispY_filt, 3);

            % Outlier removal.
            DispX = removeOutliers_2D(DispX);
            DispY = removeOutliers_2D(DispY);

            DispX = DispX{1};
            DispY = DispY{1};

            DispX_filt = filter_2D_Disp(DispX_filt, ...
                [Settings.filt_size_strain/(1 - Settings.Overlap), Settings.filt_size_strain/(1 - Settings.Overlap)], ...
                Settings.filt_pow_strain);
            DispY_filt = filter_2D_Disp(DispY_filt, ...
                [Settings.filt_size_strain/(1 - Settings.Overlap), Settings.filt_size_strain/(1 - Settings.Overlap)], ...
                Settings.filt_pow_strain);

            DispX = DispX_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter);
            DispY = DispY_filt(1+pad_filter:end-pad_filter, 1+pad_filter:end-pad_filter);

            % Radial displacement.
            DR = sqrt(DispX.^2 + DispY.^2);

            [~, idx] = min(DR(:));
            [I, J] = ind2sub(size(DR), idx);

            % Distance from the minimum of displacement.
            [x,y] = meshgrid(1:size(DR,2),1:size(DR,1));
            distance = sqrt((x-J).^2+(y-I).^2);
            distance = distance*(1-Settings.Overlap)*Settings.Resolution;

            % Strain.
            strain = (DR./distance)*100;

            % Remove the singularity at the center.
            strain(distance==0) = NaN;

            strain = inpaint_nans(strain, 3);

        end

        writematrix(strain, [File.resultsname, filesep, 'Strains', filesep, 'Strain_', temp, '.xls']);

    end

    disp(' ');
    disp('          *** End strain finder ***        ');
    disp(' ');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

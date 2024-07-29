%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%             Function used to cross-correlate a timeseries of 2D images (beads) to a reference 2D image (trypsin)             %
%              and correct their relative shift so it is lower than 0.5 pixels. The 2D images are then cropped to              %
%           the desired size. Other fluorescence or brightfield channels can also be dedrifted and cropped, with the           %
%                             same crop parameters and drift values calculated for the beads image.                            %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [structure]: structure containing the different parameters for the analysis.                                  %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       Settings [structure]: structure containing the different parameters for the analysis.                                  %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Last Revison Date: 26/07/2024                                                                                              %
%   Based on a 2D cropper created by Prof. Xavier Trepat's group: Raimon Sunyer and Xavier Trepat.                             %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [File, Settings] = cell_cropper_2D_full(Settings, File)

    disp(' ');
    disp('         *** Start image cropping ***       ');
    disp(' ');

    % Paddings used to dedrift the images.
    padding = Settings.padding;                             % Padding for the edges of the the files.

    % Pre-allocate enough memory for the largest arrays.
    [Settings.SizeY, Settings.SizeX] = size(imread(File.Name(1).Trypsin));	        % Size of the Trypsin stack.

    % Define the ROI settings depending on the image size.
    % Settings.Size_ROI{x,y} indicate the size of the final cropped images.
    % Check that the cropped images will not span more than the original images.
    if ~isfield(Settings, 'Size_ROIx')
        Settings.Size_ROIx = Settings.SizeX;
    end

    if ~isfield(Settings, 'Size_ROIy')
        Settings.Size_ROIy = Settings.SizeY;
    end

    Settings.Size_ROIx = min(Settings.Size_ROIx, Settings.SizeX);
    Settings.Size_ROIy = min(Settings.Size_ROIy, Settings.SizeY);

    % Settings.CenterRoi{X,Y} indicate the center of the final cropped images.
    if ~isfield(Settings, 'CenterRoiX')
        Settings.CenterRoiX = round(Settings.SizeX/2);
    end

    if ~isfield(Settings, 'CenterRoiY')
        Settings.CenterRoiY = round(Settings.SizeY/2);
    end

    Settings.CenterRoiX = max(1, min(Settings.CenterRoiX, Settings.SizeX));
    Settings.CenterRoiY = max(1, min(Settings.CenterRoiY, Settings.SizeY));

    % Settings.Size_ROI{x,y}_Cropper indicate the size of the subimages used by the dedrifting algorithm.
    % Check that the de-drifted sub-images will not span more than the original image.
    if ~isfield(Settings, 'Size_ROIx_Cropper')
        Settings.Size_ROIx_Cropper = Settings.SizeX;
    end

    if ~isfield(Settings, 'Size_ROIy_Cropper')
        Settings.Size_ROIy_Cropper = Settings.SizeY;
    end

    Settings.Size_ROIx_Cropper = 2*min(Settings.CenterRoiX, round(Settings.Size_ROIx_Cropper/2));
    Settings.Size_ROIy_Cropper = 2*min(Settings.CenterRoiY, round(Settings.Size_ROIy_Cropper/2));

    % Settings.Roi.{i,j}{1,2} indicate the lower and upper bounds of the i and j indices of the sub-images used by the 
    % de-drifting algorithm.
    % Check that the sub-images will not span more than the original image.
    Settings.ROI.i1 = max(1, Settings.CenterRoiY-round(Settings.Size_ROIy_Cropper/2) + 1);
    Settings.ROI.i2 = min(Settings.SizeY, Settings.CenterRoiY+round(Settings.Size_ROIy_Cropper/2));
    Settings.ROI.j1 = max(1, Settings.CenterRoiX-round(Settings.Size_ROIx_Cropper/2) + 1);
    Settings.ROI.j2 = min(Settings.SizeX, Settings.CenterRoiX+round(Settings.Size_ROIx_Cropper/2));

    % Use f{i,j}{1,2} as short names for Settings.Roi.{i,j}{1,2}.
    fi1 = Settings.ROI.i1;
    fi2 = Settings.ROI.i2;
    fj1 = Settings.ROI.j1;
    fj2 = Settings.ROI.j2;

    % Define the ROI. This contains the region with the cells of interest.
    
    % {i,j}{1,2} indicate the lower and upper bounds of the i and j indices of the sub-images that constitute the final cropped 
    % images.
    % Check that the sub-images will not span more than the original image.
    i1 = max(1, Settings.CenterRoiY - round(Settings.Size_ROIy/2)+1);
    i2 = min(Settings.SizeY, Settings.CenterRoiY + round(Settings.Size_ROIy/2));
    j1 = max(1, Settings.CenterRoiX - round(Settings.Size_ROIx/2)+1);
    j2 = min(Settings.SizeX, Settings.CenterRoiX + round(Settings.Size_ROIx/2));

    % We calculate one timepoint each "stepbeads" times.
    stepbeads = Settings.StepBeads;

    % Load the TIFF Stack object of the Trypsin image. The Tripsin image is common for the whole experiment. Load it only once. 
    tsStack = imread(File.Name(1).Trypsin);

    im1 = double(tsStack);          % Array to store the Tripsin image.

    % Here, we are storing the cropped tripsin image, but we are not using this cropped stack any more. Thus, 
    % we don't need to store it in a variable. If we don't store it and instead we crop it on-the-fly when 
    % saving it, we save a lot of RAM for large stacks.
    File.Name(1).CropTrypsin = [File.CroppPath, filesep, 'crop', '_', 'trypsin', '_', num2str(1), '.tif'];

    % When the stacks are already cropped, we will skip them. However, there is some information in the "File.mat" and 
    % "Settings.mat", obtained when the images where actually cropped, that we need to access.
    if exist([File.resultsname, filesep, 'File.mat'], 'file') && exist([File.resultsname, filesep, 'Settings.mat'], 'file')

        File_old = load([File.resultsname, filesep, 'File.mat']);

        if isfield(File_old, 'File')
            File_old = File_old.File;
        end

        save([File.resultsname, filesep, 'File_old.mat'], 'File_old', '-mat');

        n_elem_old = numel(File_old.Name);

    end

    re_calculate = 1;

    % Major loop that goes through each timepoint.
    for ii=1:stepbeads:File.NFiles.Beads

        % If the stack has already been cropped, but a correct "File.mat" file doesn't exist, cropp the stack again.
        if ~exist([File.CroppPath, filesep, 'crop','_', 'Beads', '_', num2str(ii),'.tif'], 'file') || ...
                re_calculate || (ii > n_elem_old)

            [~, temp, ext] = fileparts(File.Name(ii).Beads);
            disp(['Processing ', [temp, ext], ' of ', num2str(File.NFiles.Beads), ' images']);

            % Array to store the padded Beads image. Since it is larger than the original Beads image, we will use it to store 
            % it here too. This way, we will save  a lot of RAM for large stacks.
            im2 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding);

            % Load the TIFF Stack object.
            tsStack = imread(File.Name(ii).Beads);
            im2((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding) = double(tsStack);

            if isfield(Settings, 'Correct_Rotation') && Settings.Correct_Rotation

                [im2, theta_deg] = correct_rotation_2D(im1, ...
                                                       im2, ...
                                                       1, ...   % Window for FFT calculation. 1 means hanning
                                                       1, ...   % Method to compute max of cross-correlation. 1 means 2DCentroid
                                                       2 ...    % Size of the window
                                                       );

                File.Drift(ii).ratation_deg = theta_deg;

            end

            % In this function, we input a cropped section of matrices im1 and im2. These sub-matrices will not be used again, 
            % so we don't need to store them, and we save a lot of RAM for large stacks.
            [~, ~, shx_temp, shy_temp, MaxCorr] = disp_on_blocks_fast_2D(...
                im1(fi1:fi2, fj1:fj2), ...                  % Reference image.
                im2((fi1:fi2)+padding, (fj1:fj2)+padding), ...                  % Measurement image.
                fi2-fi1+1, ...                              % Resolution of the PIV grid, in y.
                fj2-fj1+1, ...                              % Resolution of the PIV grid. in x.
                0, ...                                      % Overlap between blocks.
                padding, ...                                % Padding for the PIV, usually 0.
                1, ...                                      % Window for FFT calculation. 1 means hanning.
            	1, ...                                      % Method to compute max of cross-correlation. 1 means 2DCentroid.
                2, ...                                      % Size of the window.
                3 ...                                       % Maximum number of iteration.
                );

            shx = round(shx_temp);
            shy = round(shy_temp);

            File.Drift(ii).shx_centroid = shx;
            File.Drift(ii).shy_centroid = shy;
            File.Drift(ii).shx_centroid_subpix = shx_temp;
            File.Drift(ii).shy_centroid_subpix = shy_temp;
            File.Drift(ii).pkh = max(MaxCorr);

            disp(['     ... current image displacement is ', '(', num2str(File.Drift(ii).shx_centroid), ',', ...
                  num2str(File.Drift(ii).shy_centroid), ') pixels... ', ...
                  ' and maximum correlation is ', num2str(File.Drift(ii).pkh)]);

            % Crop and save:
            % If correlation is smaller than this, or displacement is NaN, it means that something wrong happened to your image.
            if  (File.Drift(ii).pkh>Settings.MinCorrelationTrheshold) && (~isnan(shx_temp)) && (~isnan(shy_temp))

                File.Mask(ii).Crop = 1;             % This marks this file as good.

                padval = 0;

                % Crop beads image.
                % Shift the matrix im2 and store it again in the same variable. This will save a considerable ammount of RAM 

                indi0_delta = 1 + max(0, shy);
                indj0_delta = 1 + max(0, shx);
                indif_delta = size(im2, 1) + min(0, shy);
                indjf_delta = size(im2, 2) + min(0, shx);

                indi0_delta_2 = 1 + max(0, -shy);
                indj0_delta_2 = 1 + max(0, -shx);
                indif_delta_2 = size(im2, 1) + min(0, -shy);
                indjf_delta_2 = size(im2, 2) + min(0, -shx);

                im2(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2) = ...
                    double(im2(indi0_delta:indif_delta, indj0_delta:indjf_delta));

                im2(1:(indi0_delta_2-1), :) = padval;
                im2((indif_delta_2+1):end, :) = padval;
                im2(:, 1:(indj0_delta_2-1)) = padval;
                im2(:, (indjf_delta_2+1):end) = padval;

                % Here, we are storing the cropped beads image, but we are not using this cropped stack any more. Thus, 
                % we don't need to store it in a variable. If we don't store it and instead we crop it on-the-fly when 
                % saving it, we save a lot of RAM for large stacks.
                File.Name(ii).CropBeads = [File.CroppPath, filesep, 'crop', '_', 'Beads', '_', num2str(ii), '.tif'];
                imwrite(uint16(im2), File.Name(ii).CropBeads, 'compression', 'none');

                % Crop phase contast image. It might be a stack of images.
                if isfield(Settings.Do.Open, 'BF') && Settings.Do.Open.BF

                    padval = 0;

                    File.Name(ii).CropBF = [File.CroppPath, filesep, 'crop','_', 'BF', '_', num2str(ii),'.tif'];

                    % Load the TIFF Stack object.
                    tsStack = imread(File.Name(ii).BF);

                    % Dedrift.
                    imPhase = ones(size(im2, 1), size(im2, 2))*padval;

                    imPhase((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding) = double(tsStack);

                    if Settings.Correct_Rotation
                        imPhase = imrotate(imPhase, theta_deg, 'bilinear', 'crop');
                    end

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indif_delta = size(imPhase, 1) + min(0, shy);
                    indjf_delta = size(imPhase, 2) + min(0, shx);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indif_delta_2 = size(imPhase, 1) + min(0, -shy);
                    indjf_delta_2 = size(imPhase, 2) + min(0, -shx);

                    imPhase(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2) = ...
                        double(imPhase(indi0_delta:indif_delta, indj0_delta:indjf_delta));

                    imPhase(1:(indi0_delta_2-1), :) = padval;
                    imPhase((indif_delta_2+1):end, :) = padval;
                    imPhase(:, 1:(indj0_delta_2-1)) = padval;
                    imPhase(:, (indjf_delta_2+1):end) = padval;

                    % Stored the cropped BF image.
                    imwrite(uint16(imPhase), File.Name(ii).CropBF, 'compression', 'none');

                end

                % Crop Fluo_1 image.
                if  isfield(Settings.Do.Open, 'Fluo_1') && Settings.Do.Open.Fluo_1
      
                    padval = 0;

                    % Load the TIFF Stack object.
                    tsStack = imread(File.Name(ii).Fluo_1);

                    File.Name(ii).CropFluo_1 = [File.CroppPath, filesep, 'crop', '_', 'Fluo_1', '_', num2str(ii), '.tif'];

                    % Dedrift.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding);
                    im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding) = double(tsStack);

                    if Settings.Correct_Rotation
                        im3 = imrotate(im3, theta_deg, 'bilinear', 'crop');
                    end

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indif_delta = size(im3, 1) + min(0, shy);
                    indjf_delta = size(im3, 2) + min(0, shx);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indif_delta_2 = size(im3, 1) + min(0, -shy);
                    indjf_delta_2 = size(im3, 2) + min(0, -shx);

                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2) = ...
                        double(im3(indi0_delta:indif_delta, indj0_delta:indjf_delta));

                    im3(1:(indi0_delta_2-1), :) = padval;
                    im3((indif_delta_2+1):end, :) = padval;
                    im3(:, 1:(indj0_delta_2-1)) = padval;
                    im3(:, (indjf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo image.
                    imwrite(uint16(im3), File.Name(ii).CropFluo_1, 'compression', 'none');

                end

                % Crop Fluo_2 image.
                if  isfield(Settings.Do.Open, 'Fluo_2') && Settings.Do.Open.Fluo_2

                    padval = 0;

                    % Load the TIFF Stack object.
                    tsStack = imread(File.Name(ii).Fluo_2);

                    File.Name(ii).CropFluo_2 = [File.CroppPath, filesep, 'crop', '_', 'Fluo_2', '_', num2str(ii), '.tif'];

                    % Dedrift.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding);
                    im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding) = double(tsStack);

                    if Settings.Correct_Rotation
                        im3 = imrotate(im3, theta_deg, 'bilinear', 'crop');
                    end

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indif_delta = size(im3, 1) + min(0, shy);
                    indjf_delta = size(im3, 2) + min(0, shx);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indif_delta_2 = size(im3, 1) + min(0, -shy);
                    indjf_delta_2 = size(im3, 2) + min(0, -shx);

                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2) = ...
                        double(im3(indi0_delta:indif_delta, indj0_delta:indjf_delta));

                    im3(1:(indi0_delta_2-1), :) = padval;
                    im3((indif_delta_2+1):end, :) = padval;
                    im3(:, 1:(indj0_delta_2-1)) = padval;
                    im3(:, (indjf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo image.
                    imwrite(uint16(im3), File.Name(ii).CropFluo_2, 'compression', 'none');

                end

                % Crop Fluo image.
                if isfield(Settings.Do.Open, 'Fluo') && Settings.Do.Open.Fluo

                    padval = 0;

                    % Load the TIFF Stack object.
                    tsStack = imread(File.Name(ii).Fluo);

                    File.Name(ii).CropFluo = [File.CroppPath, filesep, 'crop', '_', 'Fluo', '_', num2str(ii), '.tif'];

                    % Dedrift.
                    im3 = zeros(Settings.SizeY + 2*padding, Settings.SizeX + 2*padding);
                    im3((1:Settings.SizeY)+padding, (1:Settings.SizeX)+padding) = double(tsStack);

                    if Settings.Correct_Rotation
                        im3 = imrotate(im3, theta_deg, 'bilinear', 'crop');
                    end

                    indi0_delta = 1 + max(0, shy);
                    indj0_delta = 1 + max(0, shx);
                    indif_delta = size(im3, 1) + min(0, shy);
                    indjf_delta = size(im3, 2) + min(0, shx);

                    indi0_delta_2 = 1 + max(0, -shy);
                    indj0_delta_2 = 1 + max(0, -shx);
                    indif_delta_2 = size(im3, 1) + min(0, -shy);
                    indjf_delta_2 = size(im3, 2) + min(0, -shx);

                    im3(indi0_delta_2:indif_delta_2, indj0_delta_2:indjf_delta_2) = ...
                        double(im3(indi0_delta:indif_delta, indj0_delta:indjf_delta));

                    im3(1:(indi0_delta_2-1), :) = padval;
                    im3((indif_delta_2+1):end, :) = padval;
                    im3(:, 1:(indj0_delta_2-1)) = padval;
                    im3(:, (indjf_delta_2+1):end) = padval;

                    % Stored the cropped Fluo image.
                    imwrite(uint16(im3), File.Name(ii).CropFluo, 'compression', 'none');

                end

            else
                disp(' Correlation is too low. File discarded from further analysis! ');
                File.Mask(ii).Crop = 0;             % Mark it as a wrong file and do not save it.

            end

        end

    end

    % Re-crop so that there are no black borders.
    x_0_crop = 1;
    x_f_crop = Settings.SizeX;

    y_0_crop = 1;
    y_f_crop = Settings.SizeY;

    for ii=1:stepbeads:File.NFiles.Beads

        x_0_crop = max(x_0_crop, 1 - File.Drift(ii).shx_centroid);
        x_f_crop = min(x_f_crop, Settings.SizeX - File.Drift(ii).shx_centroid);
        
        y_0_crop = max(y_0_crop, 1 - File.Drift(ii).shy_centroid);
        y_f_crop = min(y_f_crop, Settings.SizeY - File.Drift(ii).shy_centroid);

    end

    Settings.Size_ROIx = min(j2, x_f_crop) - max(j1, x_0_crop) + 1;
    Settings.Size_ROIy = min(i2, y_f_crop) - max(i1, y_0_crop) + 1;

    % Re-save the stacks with the re-cropping dimensions.
    imwrite(uint16(im1(max(i1, y_0_crop):min(i2, y_f_crop), max(j1, x_0_crop):min(j2, x_f_crop))), ...
            File.Name(1).CropTrypsin, 'compression', 'none');

    % Repeat for each timepoint.
    for ii=1:stepbeads:File.NFiles.Beads

        % Skip the croppinf of the stacks only if: it has previously been cropped and "File.mat" and "Settings.mat" do exist.
        if exist([File.CroppPath, filesep, 'crop', '_', 'Beads', '_', num2str(ii), '.tif'], 'file')

            % Beads stack.
            % Load the TIFF Stack object.
            tsStack = imread(File.Name(ii).CropBeads);
            im2 = double(tsStack);

            imwrite(uint16(im2((max(i1, y_0_crop):min(i2, y_f_crop))+padding, ...
                               (max(j1, x_0_crop):min(j2, x_f_crop))+padding)), File.Name(ii).CropBeads, ...
                    'compression', 'none');

            % Crop phase contast image.
            if isfield(Settings.Do.Open, 'BF') && Settings.Do.Open.BF

                % Load the TIFF Stack object.
                tsStack = imread(File.Name(ii).CropBF);
                imPhase = double(tsStack);

                % Stored the cropped BF image.
                imwrite(uint16(imPhase((max(i1, y_0_crop):min(i2, y_f_crop))+padding, ...
                                       (max(j1, x_0_crop):min(j2, x_f_crop))+padding)), File.Name(ii).CropBF, ...
                        'compression', 'none');

            end

            % Crop the Fluo_1 image.
            if  isfield(Settings.Do.Open, 'Fluo_1') && Settings.Do.Open.Fluo_1

                % Load the TIFF Stack object.
                tsStack = imread(File.Name(ii).CropFluo_1);
                im3 = double(tsStack);

                % Stored the cropped Fluo image.
                imwrite(uint16(im3((max(i1, y_0_crop):min(i2, y_f_crop))+padding, ...
                                   (max(j1, x_0_crop):min(j2, x_f_crop))+padding)), File.Name(ii).CropFluo_1, ...
                        'compression', 'none');

            end

            % Crop Fluo_2 image.
            if  isfield(Settings.Do.Open, 'Fluo_2') && Settings.Do.Open.Fluo_2

                % Load the TIFF Stack object.
                tsStack = imread(File.Name(ii).CropFluo_2);
                im3 = double(tsStack);

                % Store the cropped Fluo image.
                imwrite(uint16(im3((max(i1, y_0_crop):min(i2, y_f_crop))+padding, ...
                                   (max(j1, x_0_crop):min(j2, x_f_crop))+padding)), File.Name(ii).CropFluo_2, ...
                        'compression', 'none');

            end

            % Crop phase contrast image.
            if isfield(Settings.Do.Open, 'Fluo') && Settings.Do.Open.Fluo

                % Load the TIFF Stack object.
                tsStack = imread(File.Name(ii).CropFluo);
                im3 = double(tsStack(:, :, :));

                % Stored the cropped phase contrast image.
                imwrite(uint16(im3((max(i1, y_0_crop):min(i2, y_f_crop))+padding, ...
                                   (max(j1, x_0_crop):min(j2, x_f_crop))+padding)), File.Name(ii).CropFluo, ...
                        'compression', 'none');

            end

        end

    end

    save([File.resultsname, filesep, 'Settings.mat'], 'Settings', '-mat');
    save([File.resultsname, filesep, 'File.mat'], 'File', '-mat');

    disp(' ');
    disp('          *** End image cropping ***        ');
    disp(' ');

end

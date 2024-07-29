%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                     Function used to set the general settings of the 2D PIV code and strain calculation.                     %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%                                                                                                                              %
%   Last Revison Date: 26/07/2024                                                                                              %
%   Based on codes created by Xavier Trepat's group.                                                                           %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Settings] = set_settings_stretch(Settings)

    % Open files settings:
%     Settings.Imgfmt = 'tif';                            % Image format of input images.
    Settings.Imgfmt = 'TIF';
%     Settings.Imgfmt = 'ome.tif';
    Settings.Do.Open.Trypsin = 1;                       % Open the Trypsin images.
    Settings.Do.Open.Beads = 1;                         % Open the images of the Fluorescent Beads.
    Settings.StepBeads = 1;                             % Step between beads timepoint images 
                                                        % (1 by default, meaning that all timepoints will be considered).

    % Cropping settings:
    Settings.Do.LimitNumberOfFrames = 0;                % Limits the number of frames to process. 0 means all timepoints.

    Settings.MinCorrelationTrheshold = 0.00;            % For drift correction purposes. If trypsin and fluorescence 
                                                        % image correlate less that this value, then don't use this data point.
%     Settings.Size_ROIx = 1280;                          % This is the size of the ROIs (less than Settings.SizeX and
%     Settings.Size_ROIy = 1280;                          % multiple of Settings.resolution).
%     Settings.Size_ROIx_Cropper = 1024;                  % Size of the image used by the cropper.
%     Settings.Size_ROIy_Cropper = 1024;

%     Settings.CenterRoiX = Settings.Size_ROIx/2;         % X- and Y-coordinates of the center of the ROI used for the cropper.
%     Settings.CenterRoiY = Settings.Size_ROIy/2;         % If they are not defined, will use the center of the image.
    
    % Displacements settings:
    Settings.Resolution = 128;                          % XY resolution. Changed from 8,16,32,64. (PIV box size).

    Settings.Overlap = 0.75;                            % Overlap for PIV analysis. 0 means no overlap.
                                                        % 0.5 means overlap of half window (recommended).

    Settings.padding = 0;                               % Padding in the X and Y directions.
    Settings.Blocksizedec = 16;                         % Amount by which each side of the interrogation window shrinks at 
                                                        % every iteration, for longdisp algorithm.

    Settings.fft_window = 1;                            % Window for FFT calculation. 1 means hanning, 2 means cosine.
    Settings.corr_method = 2;                           % Method to compute max of cross-correlation. 1 means 2DCentroid, 
                                                        % 2 means Gaussian fit.
    Settings.cm_th = 0;                                 % Threshold of center of mass calculation, if 2DCentroid is used.
    Settings.pk_window = 1;                             % Half-size of the window, around the centroid peak, used to locate the 
                                                        % peak of the cross-correlation.
    Settings.iter_max = 100;                            % Maximum number of iterations after decrease.
    Settings.iter_epsilon = 0.01;                       % epsilon, maximum acceptable relative error of the iteration.
    Settings.N_conver = 10;                             % Stop after N_conver iterations without improving the convergence.
    Settings.longdisp_neigh = 10;                       % del: for each box where longdisp is required, 
                                                        % apply longdisp to the neighbours del boxes around it.
    Settings.dynamic_window = 0;                        % When determining the initial window size, we could use different 
                                                        % criteria for the different boxes. 1 if the sizes of the PIV boxes are 
                                                        % determined dynamically; 2 if all the boxes larger than 
                                                        % Settings.Resolution have the maximum initial window size; 0 if all 
                                                        % the windows have the maximum window size.
    Settings.filt_size = 2;                             % Size of the Predictor filter, in PIV resolution units,
                                                        % i.e. filter_size = filt_size/(1 - overlap).
    Settings.filt_size_strain = 2;                      % Size of the Corrector filter used for the traction calculation,
                                                        % in PIV resolution units. Same as above.
    Settings.filt_pow = 0.1;                            % Power of the Predictor filter.
    Settings.filt_pow_strain = 0.1;                     % Power of the Corrector filter.

    Settings.PixelSizeXY = 0.1625;                      % Pixel size in the x and y directions.

    Settings.maxstrain = 20;                            % For the strain plots, maximum value of the strain scale.
    Settings.maxDR = 20;                                % For the radial displacement plots, maximum value of the displacement.

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

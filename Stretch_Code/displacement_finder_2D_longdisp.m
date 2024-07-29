%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                 Function used to apply a 2D PIV to a pair of 2D stack images.                                %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [structure]: structure containing the different parameters for the analysis.                                  %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [structure]: structure containing the parameters related to paths and file names.                                 %
%                                                                                                                              %
%   Last Revison Date: 29/07/2024                                                                                              %
%   Based on a 2D PIV created by Prof. Xavier Trepat's group: Raimon Sunyer and Xavier Trepat.                                 %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%   References:                                                                                                                %
%       N/A.                                                                                                                   %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function File = displacement_finder_2D_longdisp(Settings, File)

    disp(' ');
    disp('       *** Start displacement finder ***       ');
    disp(' ');

   % Import local variables from File and Settings structures.

    nBeadsFile = File.NFiles.Beads;             % Number of timepoints to analyze.
    CroppPath = File.CroppPath;                 % Location of the cropped files.
    DispPath = File.DispPath;                   % Location of the displacement files.
    namebeads = {File.Name.Beads};              % Base of the beads image file names.
    resolution = Settings.Resolution;           % XY-resolution of the PIV.
    blocksizedec = Settings.Blocksizedec;       % Amount by which each side of the interrogation window shrinks at every 
                                                % iteration, in xy.
    overlap = Settings.Overlap;                 % Overlap of the PIV boxes. 0 means no overlap. 0.5 means overlap of half box.
                                                % In xy.

    padding = Settings.padding;                 % Padding in the X and Y directions for the PIV boxes.

    fft_window = Settings.fft_window;           % Window for FFT calculation. 1 means hanning.
    corr_method = Settings.corr_method;         % Method to compute max of cross-correlation: 1 for Centroid, 2 for Gaussian fit.
    pk_window = Settings.pk_window;             % Size of the sub-window for the peak of cross-correlation subpixel location.
    iter_max = Settings.iter_max;               % Maximum number of iterations after decrease of box size.
    iter_epsilon = Settings.iter_epsilon;       % epsilon, maximum acceptable relative error of the iteration.
    iter_N_conver = Settings.N_conver;          % N_conver: stop after N_conver iterations without improving the convergence.
    longdisp_del = Settings.longdisp_neigh;	    % del: for each box where longdisp is required, apply longdisp to the 
                                                % neighbouring "del" boxes around it. In XY.
    dynamic_window = Settings.dynamic_window;   % 1 if the sizes of the PIV windows are determined dynamically; 0 otherwise.
    filt_size = Settings.filt_size;             % Size of the Predictor filter, in window lengths, i.e. 
                                                % filter_size = filt_size*window_size. In XY.
    filt_pow = Settings.filt_pow;               % Power of the Predictor filter.

    % Load the reference (trypsin) 2D image. It is common for the whole experiment. Load it only once.
    tt = imread([CroppPath, filesep, 'crop_trypsin_', num2str(1), '.tif']);

    im1 = double(tt);                           % Array to store the Tripsin image.

    % When the stacks are already analyzed, we will skip them. However, there is some information in the "File.mat", 
    % obtained when the images where actually analyzed, that we need to access.
    if exist([File.resultsname, filesep, 'File_old.mat'], 'file')

        File_old = load([File.resultsname, filesep, 'File_old.mat']);

        if isfield(File_old, 'File_old')
            File_old = File_old.File_old;
        end

        n_elem_old = numel(File_old.Name);

    end

    re_calculate = ~exist([File.resultsname, filesep, 'File_old.mat'], 'file');

    % Loop over all the timepoints.
    for ifile = 1:nBeadsFile

        [~, temp, ext] = fileparts(namebeads{ifile});

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        temp = [temp, ext];
        temp = temp(1:end-length(Settings.Imgfmt)-1);

        beads_file = [DispPath, filesep, 'GelDisp_', temp, '.csv'];

        if exist([CroppPath, filesep, 'crop_Beads_', num2str(ifile), '.tif'], 'file') && (~exist(beads_file, 'file') || ...
           re_calculate || (ifile > n_elem_old))

        	disp( [ 'Processing crop_Beads_', num2str(ifile),'.tif', ' of ', num2str(nBeadsFile), ' images' ] ) ;

            % Image to process.
            tt = imread([CroppPath, filesep, 'crop_Beads_', num2str(ifile), '.tif']);
            im2 = double(tt);                                   % Array to store the Beads image.

            [x, y, dx, dy, pkh, ~, ~, ws0, Convergence_Mask] = disp_on_blocks_longdisp_2D(...
            	im1, ...                    % Reference image.
                im2, ...                    % Measurement image.
                resolution, ...             % Resolution, in the Y direction (rows), of the PIV grid.
                resolution, ...             % Resolution, in the X direction (columns), of the PIV grid.
                overlap, ...                % Overlap between blocks, in XY.
                padding, ...                % Padding for the PIV, in XY, usually 0.
                blocksizedec, ...           % Amount by which each side of the interrogation box shrinks at every iteration. 
                ...                         % In XY.
                fft_window, ...             % Window for FFT calculation. 1 means hanning.
                corr_method, ...            % Method to compute max of cross-correlation. 1 for 3DCentroid, 2 for Gaussian fit.
                pk_window, ...              % Size of the sub-window for the peak of cross-correlation subpixel location.
                iter_max, ...               % Maximum number of iterations after decrease of box size.
                iter_epsilon, ...           % Maximum acceptable relative error of the iteration.
                iter_N_conver, ...          % N_conver: stop after N_conver iterations without improving the convergence.
                longdisp_del, ...           % For each box where longdisp is required, apply longdisp to the neighbouring "del" 
                ...                         % boxes around it. In XY.
                dynamic_window, ...         % 1 or 2 if the sizes of the PIV windows are determined dynamically; 0 otherwise.
                filt_size, ...              % Size of the Predictor filter, in window lengths, i.e. 
                ...                         % filter_size = filt_size*window_size. 
                ...                         % In Y, X and Z respectively.
                filt_pow ...                % Power of the Predictor filter.
                );

            dx = inpaint_nans(dx);
            dy = inpaint_nans(dy);

            % Save displacements:
            tmp_header = [{'x (px)'}, {'y (px)'}, ...
                          {'u_x (px) unfiltered'}, {'u_y (px) unfiltered'}, ...
                          {'correlation (non-dim)'}, {'Initial Window Size (px)'}, {'Convergence (binary)'}];
            tmp = [x(:), y(:), dx(:), dy(:), pkh(:), ws0(:), Convergence_Mask(:)];

            vers = version('-release') ;

            if str2double( vers(1:end-1) )<2020
                writetable(array2table(tmp_header), beads_file, 'WriteVariableNames', false);

                f3 = fopen(beads_file, 'a');

                for i = 1:numel(x)
                    fprintf(f3, '%15.5e,%15.5e,%15.5e,%15.5e,%15.5e,%15.5e,%15.5e\n', ...
                            x(i), y(i), dx(i), dy(i), pkh(i), ws0(i), Convergence_Mask(i));
                end

                fclose(f3);

            else
                writetable(array2table(tmp, 'VariableNames', tmp_header), beads_file);
            end

            File.TractionSize(ifile).i = size(x,1);
            File.TractionSize(ifile).j = size(x,2);

            File.StrainSize(ifile).i = size(x,1);
            File.StrainSize(ifile).j = size(x,2);

        end

    end

    save([File.resultsname, filesep, 'File.mat'], 'File', '-mat');

    if exist([File.resultsname, filesep, 'File_old.mat'], 'file')

        delete([File.resultsname, filesep, 'File_old.mat']);

    end

    disp(' ');
    disp('        *** End displacement finder ***        ');
    disp(' ');

 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                     Script used to calculate the 2D displacement and 2D strain field of a soft substrate.                    %
%                                                                                                                              %
%   Last Revison Date: 26/07/2024                                                                                              %
%   Based on codes created by Xavier Trepat's group.                                                                           %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%------------------------------------------------------------------------------------------------------------------------------%
%                                 Parameters regarding the file names, threads and compilation.                                %
%------------------------------------------------------------------------------------------------------------------------------%

fold = '../Examples';                                 % Parent folder.
exp_fold = 'M15b';                                                                        % Experiment folder.
pathname = [fold, filesep, exp_fold];                                                   % Data folder.
resultsname = [pathname, '_Results'];                                                   % Results folder.

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                                          Setup fftw.                                                         %
%------------------------------------------------------------------------------------------------------------------------------%

% Optimize a little bit the fft algorithm (it will not matter much when all dimensions are powers of 2).
fftw( 'planner', 'exhaustive' ) ;

if exist( 'fftw_wisdom.mat', 'file' )
    load( 'fftw_wisdom.mat' ) ;
    fftw( 'wisdom', fftw_wisdom ) ;
end

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                             Calculate displacements and strains.                                             %
%------------------------------------------------------------------------------------------------------------------------------%

Settings = {};
File.pathname = pathname;
File.resultsname = resultsname;

% When the stacks are already cropped, we will skip them. However, there is some information in the "File.mat" and 
% "Settings.mat", obtained when the images where actually cropped, that we need to access.
if exist([File.resultsname, filesep, 'File.mat'], 'file') && exist([File.resultsname, filesep, 'Settings.mat'], 'file')

    Settings = load([File.resultsname, filesep, 'Settings.mat']);
    File = load([File.resultsname, filesep, 'File.mat']);

    if isfield(File, 'File')
        File = File.File;
    end

    if isfield(Settings, 'Settings')
        Settings = Settings.Settings;
    end

end

Settings = set_settings_stretch(Settings);

% Open files:
File = xOpenFiles_stretch(Settings,File);

% Crop Images:
[File, Settings] = cell_cropper_2D_full(Settings,File);

% Find Displacement Fields:
File = displacement_finder_2D_longdisp(Settings, File);

% Find Strain map.
File = strain_finder_2D(Settings, File);

% Plot strains and displacements.
File = Plot_strains_2D(Settings, File) ;

%------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------%
%                                                       Save parameters.                                                       %
%------------------------------------------------------------------------------------------------------------------------------%

% Save the fftw wisdom.
fftw_wisdom = fftw( 'wisdom' ) ;
save( 'fftw_wisdom.mat', 'fftw_wisdom' ) ;

%------------------------------------------------------------------------------------------------------------------------------%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                      Function used to set the file and folder names of the images used for the analysis.                     %
%                                                                                                                              %
%   Inputs:                                                                                                                    %
%       Settings [struct]: structure containing the different parameters for the analysis.                                     %
%       File [struct]: structure containing the file name parameters for the analysis.                                         %
%                                                                                                                              %
%   Outputs:                                                                                                                   %
%       File [struct]: structure containing the file name parameters for the analysis.                                         %
%                                                                                                                              %
%   Last Revison Date: 03/06/2024                                                                                              %
%   Author: Xavier Trepat's group.                                                                                             %
%   Modified by Manuel Gomez Gonzalez and Ernest Latorre Ibars.                                                                %
%                                                                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [File] = xOpenFiles_stretch(Settings, File)

    % Create necessary folders.
    if ~exist(File.pathname, 'dir')
        mkdir(File.pathname);
    end

    File.CroppPath = [File.resultsname, filesep, 'Croppeddata'];       % Folder where the cropped images will be stored.
    File.StrainPath = [File.resultsname, filesep, 'Strains'];          % Folder where the strain results will be stored.
    File.DispPath = [File.resultsname, filesep, 'Displacements'];      % Folder where the 2D PIV results will be stored.
    File.ImPath = [File.resultsname, filesep, 'Images'];               % Folder where the images will be stored.

    if ~exist(File.CroppPath, 'dir')
        mkdir(File.CroppPath);
    end

    if ~exist(File.StrainPath, 'dir')
        mkdir(File.StrainPath);
    end

    if ~exist(File.DispPath, 'dir')
        mkdir(File.DispPath);
    end    

    if ~exist(File.ImPath, 'dir')
        mkdir(File.ImPath);
    end    

    % Open the beads file:
    if Settings.Do.Open.Beads && (~isfield(File, 'BeadsName') || (exist(File.BeadsName, 'file') ~= 2))
        [Temp.BeadsName, File.BeadsPath] = uigetfile([File.pathname, filesep, '*.', Settings.Imgfmt], 'Open stretched image');

        % Establish the bases of the names.
        File.Base.Beads = Temp.BeadsName(1:end-length(Settings.Imgfmt)-1);
        File.Base.Beads = File.Base.Beads(1:find( isletter(File.Base.Beads), 1, 'last'));

    elseif Settings.Do.Open.Beads && isfield(File, 'BeadsName') && (exist(File.BeadsName, 'file') == 2)
        [File.BeadsPath, Temp.BeadsName, Temp.ext] = fileparts(File.BeadsName);

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        Temp.BeadsName = [Temp.BeadsName, Temp.ext];
        Temp.BeadsName = Temp.BeadsName(1:end-length(Settings.Imgfmt)-1);

        % uigetfile ends PathName with filesep, while fileparts does not.
        File.BeadsPath = [File.BeadsPath, filesep];

        % Establish the bases of the names.
        File.Base.Beads = Temp.BeadsName(1:find(isletter(Temp.BeadsName), 1, 'last'));

    end

    % Open the first trypsin file:
    if Settings.Do.Open.Trypsin && (~isfield(File, 'TrypsinName') || (exist(File.TrypsinName, 'file') ~= 2))
        [Temp.TrypsinName, File.TrypsinPath] = uigetfile([File.pathname, filesep, '*.', Settings.Imgfmt], ...
                                                         'Open reference image');

        % Establish the bases of the names.
        Temp.TrypsinName = Temp.TrypsinName(1:end-length(Settings.Imgfmt)-1);
        File.Base.Trypsin = Temp.TrypsinName(1:find(isletter(Temp.TrypsinName), 1, 'last'));

    elseif Settings.Do.Open.Trypsin && isfield(File, 'TrypsinName') && (exist(File.TrypsinName, 'file') == 2)
        [File.TrypsinPath, Temp.TrypsinName, Temp.ext] = fileparts(File.TrypsinName);

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        Temp.TrypsinName = [Temp.TrypsinName, Temp.ext];
        Temp.TrypsinName = Temp.TrypsinName(1:end-length(Settings.Imgfmt)-1);

        % uigetfile ends PathName with filesep, while fileparts does not.
        File.TrypsinPath = [File.TrypsinPath, filesep];

        % Establish the bases of the names.
        File.Base.Trypsin = Temp.TrypsinName(1:find(isletter(Temp.TrypsinName), 1, 'last'));

    end

    % How many images form the Beads stack?
    if Settings.Do.Open.Beads

        % Find all the files that start with File.Base.Beads, have the extension Settings.Imgfmt, 
        % and sort them in natural order.
        dirTrac = dir([File.BeadsPath, File.Base.Beads, '*.', Settings.Imgfmt]);
        dirTrac = natsortfiles({dirTrac.name});

        imNum = length(dirTrac);

        if Settings.Do.LimitNumberOfFrames
            imNum = min(imNum, Settings.Do.LimitNumberOfFrames);
        end

        File.NFiles.Beads = imNum;

        % Now build the names of the files:
        for i = 1:File.NFiles.Beads
            File.Name(i).Beads = [File.BeadsPath, dirTrac{i}];
        end

    end

    % How many images form the trypsin stack?
    if Settings.Do.Open.Trypsin

        File.NFiles.Trypsin = 1;

        % Now build the names of the files:
        File.Name(1).Trypsin = [File.TrypsinPath, Temp.TrypsinName, '.', Settings.Imgfmt];

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

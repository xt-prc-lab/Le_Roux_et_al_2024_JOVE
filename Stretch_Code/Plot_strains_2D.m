%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                                              %
%                                          Function used to plot the 2D strain field.                                          %
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

function File = Plot_strains_2D(Settings, File)

    disp(' ');
    disp('         *** Start strain image plotting ***       ');
    disp(' ');

    frames = File.NFiles.Beads;

    % Loop over all the timepoints.
    for k=1:frames

        [~, temp, ext] = fileparts(File.Name(k).Beads);

        % There might be problems with non-standard extensions, such as ".ome.tif". This is fixed with these two lines.
        temp = [temp, ext];
        temp = temp(1:end-length(Settings.Imgfmt)-1);

        tempFile = [File.resultsname filesep 'Displacements' filesep 'GelDisp_', temp, '.csv' ];

        if exist(tempFile, 'file') == 2

            TempData = readmatrix(tempFile);

            X = reshape(TempData(:,1), File.TractionSize(k).i, File.TractionSize(k).j);
            Y = reshape(TempData(:,2), File.TractionSize(k).i, File.TractionSize(k).j);
            DispX = reshape(TempData(:,3), File.TractionSize(k).i, File.TractionSize(k).j);
            DispY = reshape(TempData(:,4), File.TractionSize(k).i, File.TractionSize(k).j);
%             c = reshape(TempData(:,5), File.TractionSize(k).i, File.TractionSize(k).j);

            x = X(1, :);
            x = x(:);
            y = Y(:, 1);
            y = y(:);

            DR = sqrt(DispX.^2 + DispY.^2);

            scale_DR = Settings.maxDR*Settings.PixelSizeXY;

            figh1=figure;
            cla(gca);
            clf(figh1);
            set(figh1, 'visible', 'off');
            imagesc(x*Settings.PixelSizeXY, y*Settings.PixelSizeXY, DR*Settings.PixelSizeXY);
            cmap = jet(256);
            colormap(cmap);
            cb = colorbar;
            xlabel('x [\mum]');
            ylabel('y [\mum]');
            ylabel(cb, 'u_r [\mum]', 'FontSize', 14);
            caxis([0 scale_DR]);
            axis image;
            print(figh1, '-djpeg90', '-r150', [File.resultsname filesep 'Images' filesep 'DR',num2str(k),'.jpg']);
            close;

            strain = readmatrix([File.resultsname, filesep, 'Strains', filesep, 'Strain_', temp, '.xls']);

            figh1=figure;
            cla(gca);
            clf(figh1);
            set(figh1, 'visible', 'off');
            imagesc(x*Settings.PixelSizeXY, y*Settings.PixelSizeXY, strain);
            cmap = jet(256);
            colormap(cmap);
            cb = colorbar;
            xlabel('x [\mum]');
            ylabel('y [\mum]');
            ylabel(cb, 'radial strain [%]', 'FontSize', 14);
            caxis([0 Settings.maxstrain]);
            axis image;
            print(figh1, '-djpeg90', '-r150', [File.resultsname filesep 'Images' filesep 'strain',num2str(k),'.jpg']);
            close;

        end

    end

    disp(' ');
    disp('          *** strain image plotting ***        ');
    disp(' ');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

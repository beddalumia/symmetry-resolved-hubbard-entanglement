function import_colorlab()
%% Imports ALL functionality from the COLORLAB (github:bellomia/colorlab)
%
%  >> QcmP.plot.import_colorlab()
%  makes available the +rgb and +palette namespaces, plut all colortools
%
%  See also: rgb, palette, paletteshow, paletteditor, str2rgb, set_palette
%            get_palette, set_colororder

    % Relative path building...
    path_to_plot = fileparts(mfilename('fullpath'));
    path_to_repo = erase(path_to_plot,'+QcmP/+plot');
    path_to_clab = [path_to_repo, 'lib/colorlab/'];

    % Add colorlab library folders to path
    addpath(path_to_clab); colorlab.enter('quiet');
    
end
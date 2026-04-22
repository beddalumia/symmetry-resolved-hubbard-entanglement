function spectral_gif(filename,style,dt,list,var,varargin)
%% SPECTRAL_GIF: Builds a GIF for the line-evolution of the spectral functions
%
%   >> QcmP.plot.spectral_gif(filename,style,dt,list,var,varargin)
%
%  filename : filename of the complex spectral function to be plotted
%  style    : plotting style to be passed to QcmP.plot.spectral_frame()
%  dt       : delay-time for the GIF frames, in seconds
%  list     : an array of values for main line variable (could be empty!)
%  var      : an optional charvec, defining the name of the line variable [default: 'U']
%  varargin : additional options passed to QcmP.plot.spectral_frame()
%  ------------------------------------------------------------------------
    if nargin < 1
        help QcmP.plot.spectral_gif
        return
    end
    if nargin < 5
        var = 'U';
    end
    if nargin < 4 || isempty(list)
        [list, ~] = QcmP.post.get_list(var); 
    else
        list = sort(list);
    end

    Nt = length(list);

    fprintf('Start GIF building...\n\n');

    for i = 1:Nt

        % Check for <U=%f> directory
        t = list(i);
        DIR = sprintf('%s=%f',var,t);
        if ~isfolder(DIR)
            errstr = sprintf('%s_list appears to be inconsistent: ',var);
            errstr = [errstr,DIR]; %#ok
            errstr = [errstr,' folder has not been found.']; %#ok
            error(errstr);
        end

        % Enter directory
        cd(DIR); 

        % Pick and plot the requested filename
        QcmP.plot.spectral_frame(filename,style,varargin{:});
        
        % Adjust title to highlight U value
        title(DIR);

        % Exit directory
        cd('..');

        % Push frame to gif file
        [~,body,~] = fileparts(filename);
        QcmP.plot.push_frame([body,'.gif'],i,Nt,dt); close(gcf);

    end
    
    fprintf('...GIF has been built.\n\n');

end 
     
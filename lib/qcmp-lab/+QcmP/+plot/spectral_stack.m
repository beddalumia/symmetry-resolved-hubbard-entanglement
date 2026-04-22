function spectral_stack(filename,dx,dy,cmap_name,list,var,which,varargin)
%% SPECTRAL_STACK: Builds a classic stacked plot for the requested spectral
%                  tensors, reading along a computed U-driven line
%
%   >> QcmP.plot.spectral_stack(filename,dx,dy,cmap_name,list,var,which,varargin)
%
%  filename : filename of the complex spectral function to be plotted
%  dx       : horizontal step for stacking, in units of U, for proper scaling
%  dy       : vertical step for stacking, in units of U, for proper scaling
%  cmap_name: name of the desired colormap as a string (optional, see colorlab)
%  list     : an array of values for main line variable (could be empty!)
%  var      : an optional charvec, defining the name of the line variable [default: 'U']
%  which    : which function to plot? ['real' or 'imag', if not given both]
%  varargin : additional options to be passed to plotter
%
% NOTE: this plot would convey misleading information if the list is not
%       evenly spaced, in such case please consider spectral_gif instead.
%
% See also get_palette palette paletteshow
%  ------------------------------------------------------------------------
    if nargin < 1
       help QcmP.plot.spectral_stack
       return
    end

    if nargin < 4 || isempty(cmap_name)
        cmap_name = 'berlin';
    end

    if nargin < 6
        var = 'U';
    end

    if nargin < 5 || isempty(list)
        [list, ~] = QcmP.post.get_list(var); 
    else
        list = sort(list);
    end
    
    if nargin < 7 || isempty(which)
        figure("Name",'real')
        QcmP.plot.spectral_stack(filename,dx,dy,cmap_name,list,var,'real',varargin{:})
        figure("Name",'imag')
        QcmP.plot.spectral_stack(filename,dx,dy,cmap_name,list,var,'imag',varargin{:})
        return
    end

    Ny = length(list);

    fprintf('Start stacking %s spectra...\n\n',which);

    QcmP.plot.import_colorlab();
    colorlist = get_palette(cmap_name,Ny);

    for i = 1:Ny

        % Check for <U=%f> directory
        y = list(i);
        DIR = sprintf('%s=%f',var,y);
        if ~isfolder(DIR)
            errstr = sprintf('%s_list appears to be inconsistent: ',var);
            errstr = [errstr,DIR]; %#ok
            errstr = [errstr,' folder has not been found.']; %#ok
            error(errstr);
        end

        % Enter directory
        cd(DIR); 

        % Pick and plot the requested filename
        f = QcmP.plot.spectral_load(filename);

        % Plot f(x+dx)+dy curves
        switch which
            case 'real'
            plot(f.zeta+y*dx,f.real+y*dy,'Color',colorlist(i,:),varargin{:});
            case 'imag'
            plot(f.zeta+y*dx,f.imag+y*dy,'Color',colorlist(i,:),varargin{:});
        end
        % We could directly use waterfall(), but it requires Z to be a
        % meshgrid or something (a matrix) and the plotting to be done
        % outside of the loop: it would be also faster... and the plot
        % actually flexible regarding the 'viewpoint'. >>>> TODO? <<<<

        % Hold the figure handle
        hold on

        % Exit directory
        cd('..');

    end
    
    % Add label on x and y axis
    [~,body,~] = fileparts(filename);
    if      any(strfind(body,'realw'))
        xlabel('$\omega$','Interpreter','latex');
    elseif  any(strfind(body,'iw'))
        xlabel('$i\omega$','Interpreter','latex');
    end
    ylabel(body,'Interpreter','none');
    
    % Adjust title to highlight U values
    title([upper(which),' PART']);
    
    % Add legend (as a colorbar... which makes sense for an evenly spaced
    %             list, but that's already the case for the whole plot..)
    try 
        clim([min(list),max(list)]);
    catch % ver < R2022a
        caxis([min(list),max(list)]);
    end
    colormap(colorlist);
    cbar = colorbar('Location','eastoutside');
    cbar.Label.String = 'Hubbard interaction';
    
    fprintf('...DONE.\n\n');

end 
     
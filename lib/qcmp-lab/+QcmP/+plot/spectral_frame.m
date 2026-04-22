function spectral_frame(filename,style,varargin)
%% SPECTRAL_FRAME: Nicely plots a QcmPlab generated spectral tensor
%
%   >> QcmP.plot.spectral_frame(filename,style,varargin)
%
%  filename : filename of the complex spectral function to be plotted
%  style    : plotting style to be passed to internal 1d-plotter ['line','area','scatter']
%  varargin : additional options passed to internal 1d-plotter (e.g. 'FaceColor', etc...)
%  ------------------------------------------------------------------------

    % Load complex function of real/imag frequency z, according to QcmPlab conventions
    F = QcmP.plot.spectral_load(filename); % F.zeta, F.real, F.imag

    % Plot both real and imag parts of F
    if ~exist('style','var'); style = 'area'; end
    cplot(F.zeta,F.real,style,varargin{:}); hold on
    cplot(F.zeta,F.imag,style,varargin{:});

    % Legend, Labels
    legend('real','imag');
    if      any(strfind(filename,'realw'))
                xlabel('$\omega$','Interpreter','latex');
    elseif  any(strfind(filename,'iw'))
                xlabel('$i\omega$','Interpreter','latex');
    end
    [~,body,~] = fileparts(filename);
    ylabel(body,'Interpreter','none');

    % Limits
    xlim([F.zeta(1),F.zeta(end)]);
    try % To deal with self-energies
    ylim([min(min(F.real),min(F.imag)),max(max(F.real),max(F.imag))]);
    catch
    ylim([-1,1]);
    end

end

%% Custom overloaded 1d-plotter
function cplot(x,y,style,varargin)

    if     strcmp(style,'area')

                area(x,y,'FaceAlpha',0.5,varargin{:});

    elseif strcmp(style,'scatter')

                scatter(x,y,'filled',varargin{:});

    else   % default: line-plot

                plot(x,y,varargin{:});

    end

end

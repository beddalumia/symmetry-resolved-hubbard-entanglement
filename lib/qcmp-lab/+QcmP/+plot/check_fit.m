function check_fit(filter,layout)
%% CHECK_FIT loads all fit_*.ed files in current directory and makes a plot
%
%  >> QcmP.plot.check_fit(<filter>,<layout>)
%
%     filter: an optional string, giving a "hint" to filter the files
%     • 'ineq0001' --> all files in the "fit_*ineq0001.ed" form
%     • 's1'       --> all files in the "fit_*s1*.ed" form
%     • <glob>     --> any glob or regexp, but keep in mind that would be
%                      sandwiched as "fit-*" + <glob> + "*.ed"
%
%     layout: an optional string, to allow choosing between...
%     • 'vertical'   -->  vertical subplot concatenation (default)
%     • 'horizontal' -->  horizontal subplot concatenation
%     • 'matrix'     -->  matrix-like subplot concatenation 

    if nargin < 1
        filter = "_";
    end

    if nargin < 2
        layout = 'vertical';
    end

    if isa(filter,"char")
        glob = ['fit_*',filter,'*.ed'];
    elseif isa(filter,"string")
        glob = "fit_*" + filter + "*.ed";
    end

    fit_list = dir(glob);
    n_fit = length(fit_list);
    n_plt = n_fit * 2;
    quadL = ceil(sqrt(n_plt));
    if n_fit==0
        fprintf(2,"ERROR: No fit files found.\n")
        return
    end

    figure_handle = figure();

    switch layout
        case 'matrix'
            handle = tiledlayout(quadL,quadL);
            set(figure_handle, ...
                'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
        case 'horizontal'
            handle = tiledlayout(2,n_fit);
            set(figure_handle, ...
                'Units', 'Normalized', 'OuterPosition', [0, 0.25, 1, 0.5]);
        case 'vertical'
            handle = tiledlayout(n_fit,2);
            set(figure_handle, ...
                'Units', 'Normalized', 'OuterPosition', [0.25, 0, 0.5, 1]);
        otherwise
            fprintf(2,'\nERROR: invalid form-factor.\n\n')
            close(figure_handle); help QcmP.plot.check_fit;
            return
    end

    QcmP.plot.import_colorlab

    for i = 1:n_fit
        % Collect fit data
        fit_name = fit_list(i).name;
        fit_data = load(fit_name);
        Freq = fit_data(:,1);
        ImFG = fit_data(:,2);
        ImFT = fit_data(:,3);
        ReFG = fit_data(:,4);
        ReFT = fit_data(:,5);
        % Pre-process file name
        fit_name = erase(fit_name,'fit_');
        fit_name = erase(fit_name,'.ed');
        % Plot real in a new tile
        nexttile
        plot(Freq,ReFG,'-','LineWidth',2,'Color',str2rgb('deep purple'))
        hold on
        plot(Freq,ReFT,'-.','LineWidth',2,'Color',str2rgb('light purple'))
        hold off
        title(fit_name,"Interpreter","none")
        xlabel("i$\omega$","Interpreter","latex")
        ylabel("REAL")
        legend(["FG","Fit"])
        % Plot imag in a new tile
        nexttile
        plot(Freq,ImFG,'-','LineWidth',2,'Color',str2rgb('deep blue'))
        hold on
        plot(Freq,ImFT,'-.','LineWidth',2,'Color',str2rgb('sky blue'))
        hold off
        title(fit_name,"Interpreter","none")
        xlabel("i$\omega$","Interpreter","latex")
        ylabel("IMAG")
        legend(["FG","Fit"])
    end
end
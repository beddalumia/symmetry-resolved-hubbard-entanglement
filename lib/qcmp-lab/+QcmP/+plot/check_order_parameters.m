function check_order_parameters(var,suffix)
%% Loads all order parameters in a U-line and builds a stacked-plot
%  >> QcmP.plot.check_order_parameters(<var>,<suffix>)
%     var    : an optional charvec, defining the name of the line variable [default: 'U']
%     suffix : an optional charvec, handling inequivalent filename endings
    if nargin < 1
        var = 'U';
    end
    if nargin < 2
        suffix = [];
    end
    [ids,O_cell,X_list] = QcmP.post.order_parameter_line(var,suffix);
    O_matrix = cell2mat(O_cell);
    s = stackedplot(X_list,O_matrix);
    xlabel(var)
    s.DisplayLabels = ids;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0, 0.2, 1]);
end
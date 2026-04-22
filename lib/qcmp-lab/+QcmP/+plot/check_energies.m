function check_energies(var)
%% Loads all impurity energy values in a U-line and builds a stacked-plot
%  >> QcmP.plot.check_energies(<var>)
%     var    : an optional charvec, defining the name of the line variable [default: 'U']
    if nargin < 1
        var = 'U';
    end
    [ids,E_cell,X_list] = QcmP.post.energy_line(var);
    E_matrix = [cell2mat(E_cell)];
    try
        Kin_vect = QcmP.post.kinetic_line(X_list);
        E_matrix = [E_matrix,Kin_vect];
    catch
        disp("No kinetic energy available, skipping it.")
    end
    s = stackedplot(X_list,E_matrix);
    xlabel(var)
    try
        s.DisplayLabels = [ids,"Kin"];
    catch
        s.DisplayLabels = ids;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.4, 0, 0.2, 1]);
end
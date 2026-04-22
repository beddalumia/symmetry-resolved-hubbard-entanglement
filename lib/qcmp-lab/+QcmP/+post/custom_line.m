function [C_list,X_list] = custom_line(filename,var,X_list)
%% Getting a list of values, from custom files containing a single real value.
%
%       [C_list,X_list] = QcmP.post.custom_line(filename,var,X_list)
%
%  filename : a *required* str/char, to provide the custom filename to look for
%  var      : an optional char for the name of the line variable [default: 'U'] 
%  C_list   : a double-precision array, forall X, giving all the retrieved values
%  X_list   : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<2)
        var = 'U';
    end
    if ~exist('X_list','var') || isempty(X_list)
       [X_list, ~] = QcmP.post.get_list(var); 
    else
       X_list = sort(X_list);
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    C_list = zeros(Nx,1);
    for i = 1:length(X_list)
        x = X_list(i);
        DIR= sprintf('%s=%f',var,x);
        if ~isfolder(DIR)
            errstr = sprintf('%s_list appears to be inconsistent: ',var);
            errstr = [errstr,DIR]; %#ok
            errstr = [errstr,' folder has not been found.']; %#ok
            error(errstr);
        end
        cd(DIR);
        if not(isfile(filename))
            C_list(i) = NaN;
        else
            C_list(i) = load(filename);
        end
        cd('..');
    end
    [~,filename] = fileparts(filename);
    filename = filename + ".txt"; % Always use txt extension
    QcmP.post.writematrix(C_list,filename,'Delimiter','tab');
end


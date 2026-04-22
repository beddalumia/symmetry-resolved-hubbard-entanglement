function [L_list,X_list] = luttinger_line(var,suffix,X_list)
%% Getting a list of luttinger invariants, from directories.
%
%       [L_list,X_list] = QcmP.post.luttinger_line(var,suffix,X_list)
%
%  var    : an optional char for the name of the line variable [default: 'U'] 
%  suffix : an optional charvec, handling inequivalent filename endings
%  L_list : a float-array, forall X, giving all the luttinger integral vals
%  X_list : an optional array of values for the line variable: where to search
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(nargin<1)
        var = 'U';
    end
    if(nargin<2)
        suffix = [];
    end
    if(nargin<3 || isempty(X_list))
        [X_list, ~] = QcmP.post.get_list(var); 
    else
        X_list = sort(X_list);
    end
    % Then we can proceed spanning all the U-values
    Nx = length(X_list);
    L_list = zeros(Nx,1);
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
        if(~isempty(suffix))
            filename = ['luttinger_',suffix,'.dat'];
        else
            filename = 'luttinger.dat';
        end
        if not(isfile(filename))
            L_list(i) = NaN;
        else
            L_list(i) = load(filename);
        end
        cd('..');
    end
    if(~isempty(suffix))
        filename = ['luttinger_line_',suffix,'.txt'];
    else
        filename = 'luttinger_line.txt';
    end
    QcmP.post.writematrix(L_list,filename,'Delimiter','tab');
end
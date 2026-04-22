function [Z_list,X_list] = zeta_line(var,suffix,X_list)
%% Getting a list of Z-weight values, from directories.
%
%       [Z_list,X_list] = QcmP.post.zeta_line(var,suffix,X_list)
%
%  var    : an optional char for the name of the line variable [default: 'U'] 
%  suffix : an optional charvec, handling inequivalent filename endings
%  Z_list : a float-array, forall X, giving all the quasiparticle weight values
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
    Z_list = zeros(Nx,1);
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
            filename = ['zeta_last_',suffix,'.ed'];
        else
            filename = 'zeta_last.ed';
        end
        if not(isfile(filename))
            Z_list(i) = NaN;
        else
            Z_list(i) = load(filename);
        end
        cd('..');
    end
    if(~isempty(suffix))
        filename = ['zeta_',suffix,'.txt'];
    else
        filename = 'zeta.txt';
    end
    QcmP.post.writematrix(Z_list,filename,'Delimiter','tab');
end
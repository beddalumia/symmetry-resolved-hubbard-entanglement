function refresh_line(EXE,doMPI,Nnew,var,where,varargin)
%% Refreshes a pre-existent calculation by Nnew loops [whole-line]
%
%   QcmP.run.refresh_line(EXE,doMPI,Nnew,var,where,varargin)
%
%   EXE                 : Executable driver [mandatory]
%   doMPI               : Flag to activate OpenMPI [default: true]
%   Nnew                : How much loops to add in the refresh [default: 1]
%   var                 : Main variable of the line [$(var)=%f directories]
%   where               : Optional input defining whichˆ points to refresh
%   varargin            : Set of fixed control parameters ['name',value]
%
% ˆ['all': all subfolders, 'conv': read $(var)_conv.txt, 'nonconv': all - conv]
%  >> default string is 'nonconv'
%  >> you shall pass even an array of floats instead of a string, in that
%     case those would be the values to be refreshed (at your own risk).

%% Print docstring if no input is provided
if nargin < 1
   help QcmP.run.refresh_line 
   return
end

%% Other optional input wrangling
if nargin < 2
   doMPI = true;
end
if nargin < 3
   Nnew = 1;
end
if nargin < 4
   var = 'U';
end

%% Default behavior if no provided where
if nargin < 5
   mode = 'nonconv';
elseif isnumeric(where)
   mode = 'custom';
else
   mode = where; 
end

%% Open file to write/update converged $(var)-values
fconv_name = sprintf('%s_conv.txt',var);
if(not(isfile(fconv_name)))
   warning('You may have passed a wrong main variable name!')
end
fileID_conv = fopen(fconv_name,'a');

%% Retrieve the list of *converged* U-values
x_conv = unique(sort(load(fconv_name)));

%% Retrieve the list of all available U-values
x_list = QcmP.post.get_list(var);

%% Define which list of U-values refresh
switch mode
    
    case 'all'
        
        X = x_list;
        
    case 'conv'
        
        X = x_conv;
        if isempty(X) 
           warning('There are no converged points!')
           return
        end
        
    case 'nonconv'
        
        X = setdiff(x_list,x_conv);
        if isempty(X) 
           warning('There are no unconverged points!')
           return
        end
        
    case 'custom'
        
        X = where;
        
    otherwise

        help QcmP.run.refresh_line
        error('Invalid set of points (what you passed as where?).')
        
end

%% Tell the user what we are refreshing...
fprintf('Refreshing the following points:')
disp(X)

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(X)

    DIR = sprintf('%s=%f',var,X(i));   % Define the $(var)-folder name.

    if isfolder(DIR)
        cd(DIR);                  % Enter the $(var)-folder (if it exists)
    else
        errstr = sprintf('%s_list appears to be inconsistent: ',var);
        errstr = [errstr,DIR]; %#ok
        errstr = [errstr,' folder has not been found.']; %#ok
        error(errstr);
    end

    unconverged = QcmP.run.refresh_point(EXE,doMPI,Nnew,varargin{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
    if not(unconverged) && isempty(find(x_conv == X(i),1))
        fprintf(fileID_conv,'%f\n', X(i));      % Update x_conv, only
    end                                         % if *newly* converged
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cd ..                           % Exit the $(var)-folder

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fileID_conv);

end




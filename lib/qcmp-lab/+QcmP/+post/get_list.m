function [flist, strlist] = get_list(VARNAME)
%% Getting a list of variable values, from directories.
%
%       [flist, strlist] = QcmP.post.get_list(VARNAME)
%
%  VARNAME: a string, identifying the listed variable (e.g. 'U')
%  flist: an array of float_values (e.g. U=[:] )
%  strlist: an array of dir_name strings (e.g. ['U=%f'] )
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    subentries = dir('.'); % Retrieves every folder and file inside current
    subfolders = subentries([subentries(:).isdir]); % Keeps only subfolders
    subfolders = subfolders(~ismember({subfolders(:).name},{'.','..'}));
    N = length(subfolders); flist = zeros(N,1); 
    try
        strlist = strings(N,1);
    catch % GNU Octave
        strlist = cell(N,1);
    end
    for i = 1:N
        DIR = subfolders(i).name; % Let's get the indexed string...
      try
        flist(i) = sscanf(DIR, [VARNAME,'=%f']); %...and extract the value!
        strlist(i) = DIR;
      catch
        disp(' ')
        disp(DIR)
        disp('> spurious directory: will be discarded')
        flist(i) = NaN;
      end
    end
    % Prune away all the NaNs
    strlist = strlist(~isnan(flist));
    flist = flist(~isnan(flist));
    % We need to sort the lists by floats (not strings, as it is now)
    [flist, sortedIDX] = sort(flist); strlist = strlist(sortedIDX);
end

function [conv_ratio,tot_points,unconv_pts] = inspect_line(var,print_list)
%% To monitor convergence ratio and total number of points at runtime.
%
%   [conv_ratio,tot_points,unconv_pts] = QcmP.run.inspect_line(print_list)
%
%   conv_ratio   : Percentage of converged points
%   tot_points   : Total number of computed points
%   unconv_pts   : List of unconverged points [$(var)-values]
%   var          : Name of the line variable [default: 'U']
%   print_list   : Flag to activate printing of unconv_pts [default: don't]
%
%   It also logs the info to stdout, if nargout == 0 [interactive usage]

% Define the default $(var)
if nargin < 1
   var = 'U';
end

% Load files to read (un)converged U-values
X_conv = load(sprintf('%s_conv.txt',var));
X_list = load(sprintf('%s_list.txt',var));

% Compute output variables
tot_points = length(X_list); 
unconv_pts = setdiff(X_list,X_conv);
conv_ratio = length(X_conv)/tot_points*100;

% Print to standard output
if nargout == 0
   fprintf('\n')
   fprintf('Total number of points: %d (so far)\n', tot_points)
   fprintf('Converged points ratio: %d%% (so far)\n', round(conv_ratio))
   fprintf('\n')
   if nargin > 1 && print_list
      disp Unconverged:
      disp ' '
      disp(unconv_pts)
      if isempty(unconv_pts)
         disp('> None')
      end
   end
end

end
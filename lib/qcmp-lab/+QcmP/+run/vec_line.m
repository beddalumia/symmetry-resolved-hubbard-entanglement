function vec_line(EXE,doMPI,old,values,var,varargin)
%% Runs a $(var)-line with no other convergence control than writing $(var)_conv.txt
%
%   QcmP.run.dry_line(EXE,doMPI,old,values,var,varargin)
%
%   EXE                 : Executable driver
%   doMPI               : Flag to activate OpenMPI
%   old                 : Restart point [NaN or empty -> no restart]
%   values              : Input values for $(var) [a generic array]
%   var                 : Main variable of the line [default: 'U']
%   varargin            : Set of fixed control parameters ['name',value]

%% Print docstring if no input is provided
if nargin < 1
   help QcmP.run.dry_line
   return
end

%% Handle default value for var
if nargin < 4
   var = 'U';
end

xlist = fopen(sprintf('%s_list.txt',var),'a');
xconv = fopen(sprintf('%s_conv.txt',var),'a');

%% Phase-Line: single loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(values)

   x = values(i); 
   if i > 1
      old = values(i-1);   % update the restarting point
   end
   unconverged = QcmP.run.single_point(EXE,doMPI,var,x,old,varargin{:});

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %% HERE WE CATCH A FAILED (unconverged) DMFT LOOP
   if unconverged
      fprintf(xlist,'%f\n', NaN); % Update xlist with a null value
   else
      fprintf(xlist,'%f\n', x);   % Update xlist with actual value
      fprintf(xconv,'%f\n', x);   % Update xconv only if converged
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(xlist); fclose(xconv);

end



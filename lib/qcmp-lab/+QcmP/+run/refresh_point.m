function unconverged = refresh_point(EXE,doMPI,Nnew,varargin)
%% Refreshes a pre-existent calculation by Nnew loops [single-point]
%
%   unconverged = QcmP.run.refresh_point(EXE,doMPI,Nnew,varargin)
%
%   EXE                 : Executable driver
%   doMPI               : Flag to activate OpenMPI
%   Nnew                : How much loops to add in the refresh
%   varargin            : Set of fixed control parameters ['name',value]

%% Print docstring if no input is provided
if nargin < 1
   help QcmP.run.refresh_point 
   return
end

%% Move used.input into input [crucial: we do not know about touched ctrl-vars]
usedinput  = dir([pwd, '/used.*']);
parts = strsplit(usedinput.name,'.');
newname = cell2mat(join(parts(2:end),'.'));
copyfile(usedinput.name,newname); % [NOT MOVEFILE PLEASE: corruption prone!]

%% Move *.used into *.restart [desirable: first new loop == last old loop]
usedpack = dir([pwd, '/*.used']);
for i = 1:numel(usedpack)
    file = fullfile(usedpack(i).folder, usedpack(i).name);
    [folder, name, extension] = fileparts(file); %#ok (mlint suppression)
    extension = '.restart';
    movefile(file, fullfile(folder, [name, extension]));
end

%% Reset the unconverged status
delete('ERROR.README')

%% Run FORTRAN code (already compiled and added to PATH!) %%%%%%%%%%%%%%%%%
if doMPI
    mpi = 'mpirun ';                        % Control of MPI
else                                        % boolean flag...
    mpi = [];
end
NLOOP = sprintf(' nloop=%d',Nnew);          % OVERRIDE of #{loops}
VAR = [];                                   % OVERRIDE of
for i = 1:2:(length(varargin)-1)            % ALL
    VARname = varargin{i};                  % the
    VARsval = char(string(varargin{i+1}));  % OTHER
    VAR = [VAR,' ',VARname,'=',VARsval];    % PARAMETERS
end	
out = ' | tee -a LOG_refresh.out';          % STDOUT destination
sys_call = [mpi,EXE,NLOOP,VAR,out]          %#ok (we want to print)
tic
system(sys_call);                           % Fortran-call
chrono = toc;
file_id = fopen('LOG_refresh.time','w');    % Write on time-log
fprintf(file_id,'%f\n', chrono);
fclose(file_id);

%% HERE WE CATCH A FAILED (unconverged) DMFT LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
errorfile = 'ERROR.README';
unconverged = isfile(errorfile);
if(unconverged)
   fprintf(2,['\n< ',sys_call,' > point has not converged!\n\n']); beep;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end




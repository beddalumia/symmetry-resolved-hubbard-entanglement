function push_frame(gifname,iframe,nframe,dt,infig)
%% PUSH_FRAME: appends a frame to a .gif file [high quality]
%
%   >> QcmP.plot.push_frame(gifname,iframe,nframe,dt,infig)
%
% gifname   :   filename for the gif, should end with .gif
% iframe    :   id number of the given frame, as an integer
% nframe    :   total number of expected frames, as an integer
% dt        :   duration for the given frame, in seconds
% infig     :   optional handle for the figure to print
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                                                        if nargin < 5
                                                            infig = gcf;
                                                        end    
    % Write info to stdout
    info = sprintf('Appending %d-th frame of %d to ',iframe,nframe);
    info = [info,'<',gifname,'>','\n'];
    fprintf(info);
    % Capture the plot as an image
    im = print(infig,'-RGBImage');
    % Suitable conversion
    [ind,cm] = rgb2ind(im,256,'nodither');
    % Write to the GIF file
    if iframe == 1
        imwrite(ind,cm,gifname,'gif','Loopcount',inf,'DelayTime',dt);
    else
        imwrite(ind,cm,gifname,'gif','WriteMode','append','DelayTime',dt);
    end
    
end

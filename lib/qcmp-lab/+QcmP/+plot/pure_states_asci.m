function pure_states_asci(file_rdm,file_states,indices)
%% Visualizing all information about pure states, for a given point.
%
%       >> pure_states_asci(file_rdm,file_states,indices)
%
%  ➜ it draws a treemap with patches-area proportional to p(i)*|c{i}(j)|^2,
%    where i spans the pure states, and p(i) is the associated probability
%    and j spans the Fock-space basis, as defined by the ASCI solver.
%    The c{i}(j) are basis expansion coefficients for the pure states, so
%    we have each patch with an area that intuitively represents how much
%    a configuration state (from the basis) contributes to the groundstate.
%    Nevertheless we retain information about which are the pure states, by
%    grouping such patches in bigger rectangles, each one corresponding to
%    a single pure state, with area proportional to p(i). Suitable labels
%    are assigned to the patches, when big enough to accomodate text.
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % Input checking
    if(nargin<1)
        help QcmP.plot.pure_states_asci
        return
    end
    
    % Retrieve pure state data
    %[c,p] = QcmP.post.get_pure_states(suffix,check);
    [c,p] = get_pure_states(file_rdm);
    
    % We set Nlat to the passed input (there's no way to infer from the files, yet)
    Nlat = length(indices);
    % Then we can determine number of local orbitals, too :)
    Nrdm = length(p);

    % Build the nested treemap: outer patches ~ p(i), inner ~ |c{i}(j)|^2
    outer_tree = treemap(p);
    QcmP.plot.import_colorlab
    colormap = get_palette('Spectral',Nrdm);
    rng(7); % fix the seed!
    colormap = colormap(randperm(Nrdm),:); % avoid too similar colors for similar states
    abs_weight = zeros(Nrdm,1);
    labels = cell(Nrdm^2,1);
    for i = 1:Nrdm
        disp(i)
        expansion = c{i}.*conj(c{i});
        inner_tree = treemap(expansion,outer_tree(3,i),outer_tree(4,i));
        inner_tree(1,:) = inner_tree(1,:) + outer_tree(1,i);
        inner_tree(2,:) = inner_tree(2,:) + outer_tree(2,i);
        % Assign labels, when appropriate (p(i)*|c{i}(j)|^2 > 1%)
        labels = build_ket(file_states,c{i}(:));
        fprintf('p_%d = %f\n\n',i,p(i))
        for j = 1:Nrdm
            abs_weight(j) = p(i) * norm(c{i}(j))^2;
            fprintf('  w_%d = %f\n',j,norm(c{i}(j))^2)
            fprintf('  p_%d * w_%d = %f\n\n',i,j,abs_weight(j))
            if abs_weight(j) < 0.01
                labels{j} = '';
            end
            % Uncomment to debug treemap areas:
            % labels{j} = num2str(abs_weight(j));
        end
        plotRectangles(inner_tree,labels,colormap);
    end
    outline(outer_tree,'-',1.5);

    function kets = build_ket(file_state,coefficients)
        %% BUILD_KET : Puts together a pretty label for a pure state component
        %              (following conventions adopted by ASCI solver)
        file_stream = fopen(file_state,'r');
        data = textscan(file_stream,'%s');
        dets = data{:}(3:2:end);
        for k = 1:Nrdm % How to vectorize this crap remains to be found...
            dets{k} = [dets{k}(indices(1:2))];%,char(10),dets{k}(indices(3:4))];
        end
        kets = string(dets);
        kets = strrep(kets,'u','\uparrow');
        kets = strrep(kets,'d','\downarrow');
        kets = strrep(kets,'2','\leftrightarrow');
        kets = strrep(kets,'0','\bullet');
        %kets = "\mid " + kets + " \rangle";
        kets = string(coefficients) + " \times " + kets;%char(10) + kets;
        fclose(file_stream);
    end
end

function rectangles = treemap(data,w,h)
    %%  TREEMAP : Partition a rectangular area according to data
    %   rectangles = treemap(data,w,h) partitions the rectangle [0,0,w,h]
    %   into an array of rectangles, one for each element of the vector "data"
    %   The areas of the rectangles are proportional to the values in data and
    %   the aspect ratio of the rectangles are close to one. If w and h are
    %   not given we would assume w = h = 1, if neither data is provided some
    %   random sample would be generated to provide a demonstrative example.
    %   Finally if no output is explicitly requested we'll automatically make
    %   a plot for you.
    %
    %   The algorithm used is as follows. The data are sorted from largest to
    %   smallest. We then lay out a row or column along the shortest side of
    %   the remaining rectangle. Blocks are added to this new row or column
    %   until adding a new block would worsen the average aspect ratio of the
    %   blocks in the row or column. Once this row or column is laid out, we
    %   recurse on the remaining blocks and the remaining rectangle.
    %
    %   Examples:
    %    treemap(rand(20,1))
    %
    %    r = treemap(1:15,1.6,1);
    %    plotrectangles(r)
    %    outline(r)
    %
    %   Copyright 2007-2013 The MathWorks, Inc. and 2022 Gabriele Bellomia
    %
    if (nargin == 0)
        % make up some random data
        data = abs(randn(1,80)) .^ 2;
    end
    % if you don't specify a rectangle, we use the unit square
    if(nargin < 2)
        w = 1;
        h = 1;
    end
    % this is the ratio of rectangle area to data values
    areaRatio = (w * h) / sum(data);
    % we lay out the date from largest to smallest
    [v,sortIndex] = sort(data,'descend');
    % this chooses the number of blocks in each row or column
    p = partitionArea(v,w,h);
    % this leaves us with the task of assigning rectangles to each element.
    % storage for the result
    rectangles = zeros(4,length(data));
    % these hold the origin for each row
    oX = 0;
    oY = 0;
    % the index of the first value in v for the current row
    first = 1;
    % where to put the resulting rectangles in the rectangles array
    resultIndex = 1;
    % for each row to layout...
    for i = 1:length(p)
        % which values are we placing
        last = first + p(i) - 1;
        chunk = v(first:last);
        % for the next iteration
        first = last + 1;
        % how much area should each block have?
        blockArea = chunk * areaRatio;
        % how much area must the entire column have?
        columnArea = sum(blockArea);
        % the origin for the blocks starts as the origin for the column
        blockX = oX;
        blockY = oY;
        % are we laying out a row or a column?
        if((w < h)) % a row
            % so how thick must the row be?
            columnHeight = columnArea / w;
            % lets place each value
            for j = 1:p(i)
                % so how wide should it be?
                blockWidth = blockArea(j) / columnHeight;
                % remember the rectangle
                rectangles(:,sortIndex(resultIndex)) = [blockX,blockY,blockWidth,columnHeight];
                resultIndex = resultIndex + 1;
                % move the origin for the nextblock
                blockX = blockX + blockWidth;
            end
            % move the origin for the next row
            oY = oY + columnHeight;
            h = h - columnHeight;
        else % a column
            columnWidth = columnArea / h;
            % lets place each value
            for j = 1:p(i)
                % so how high should it be?
                blockHeight = blockArea(j) / columnWidth;
                % if one corner is at oX,oY, where is the opposite corner?
                rectangles(:,sortIndex(resultIndex)) = [blockX,blockY,columnWidth,blockHeight];
                resultIndex = resultIndex + 1;
                % move the origin for the nextblock
                blockY = blockY  +  blockHeight;
            end
            % move the origin for the next row
            oX = oX + columnWidth;
            w = w - columnWidth;
        end
    end
    % show the results if they are not returned
    if(nargout == 0)
        plotRectangles(rectangles)
    end
end
%
% contains
%
function partition = partitionArea(v,w,h)
    % Returns a vector telling us how many blocks go in each column or
    % row. Normalization: sum(partition) == length(v).
    % The rest of the code only wands to deal with long side and short
    % side. It is not interested in orientation (w or h)
    %
    % Copyright 2007-2013 The MathWorks, Inc. and 2022 Gabriele Bellomia
    %
    if(w > h)
        longSide = w;
        shortSide = h;
    else
        longSide = h;
        shortSide = w;
    end
    % this is the ratio of value units to associated area.
    areaRatio = (w * h) / sum(v);
    % we want to minimize cost
    bestAverageAspectRatio = inf;
    % we will return an array of integers that tell how many blocks to
    % place in each row (or column)
    partition = [];
    % How many blocks should go in the next column of rectangles?
    % i is our current guess. We'll keep incrementing it until aspect ratio
    % (cost) gets worse.
    for i = 1:length(v)
        columnTotal = sum(v(1:i));
        columnArea = columnTotal * areaRatio;
        columnWidth = columnArea / shortSide;
        % this is the cost associated with this value of i.
        sumOfAspectRatio = 0;
        % for each block in the column
        for j = 1:i
            % what is the aspect ratio
            blockArea = v(j) * areaRatio;
            blockHeight = blockArea / columnWidth;
            aspectRatio = blockHeight / columnWidth;
            if(aspectRatio < 1)
                aspectRatio = 1 / aspectRatio;
            end
            sumOfAspectRatio = sumOfAspectRatio + aspectRatio;
        end
        averageAspectRatio = sumOfAspectRatio/i;
        % if the cost of this value of i worse than for (i-1) we are done
        % and we will use i-i blocks in this column.
        if(averageAspectRatio >= bestAverageAspectRatio)
            if(partition < i) % if we are not done yet, we recurse
                p = partitionArea(v(i:end),shortSide,longSide - columnWidth);
                partition = [partition,p];
            end
            return
        end
        bestAverageAspectRatio = averageAspectRatio;
        partition = i;
    end
end

function plotRectangles(rectangles,labels,colors)
    %%  PLOTRECTANGLES : Create patches to display rectangles
    %
    %   >> plotRectangles(rectangles,labels,colors)
    %
    %   rectangles is a 4-by-N matrix where each column is a rectangle in
    %   the form [x0,y0,w,h]. x0 and y0 are the coordinates of the lower
    %   left-most point in the rectangle.
    %
    %   Example:
    %    cla
    %    plotRectangles([0 0 1 1; 1 0 1 1; 0 1 2 1]')
    %
    %   Copyright 2007-2013 The MathWorks, Inc. and 2022 Gabriele Bellomia
    %
    %   See also OUTLINE.
    
    if(nargin < 2)
        labels = [];
    end
    
    if(nargin < 3)
        colors = rand(size(rectangles,2),3).^0.5;
    end
    
    % make a patch for each rectangle
    
    for i = 1:size(rectangles,2)
        r = rectangles(:,i);
        if isnan(r(3)) || isnan(r(4)) || r(3)<1e-12 || r(4)<1e-12
           continue 
        end
        xPoints = [r(1), r(1), r(1) + r(3),r(1) + r(3)];
        yPoints = [r(2), r(2)+ r(4), r(2)+ r(4),r(2)];
        if r(3) >= r(4)
           angle = 0;
        else
           angle = 90;
        end
        patch(xPoints,yPoints,colors(i,:),'EdgeColor','none');
        if(~isempty(labels))
            text(r(1) + r(3)/2,r(2) + r(4)/2, 1, labels{i}, ...
                'VerticalAlignment','middle', ...
                'HorizontalAlignment','center', ...
                'Rotation',angle, ...
                'Interpreter','tex', ...
                'FontSize',30*max(r(3),r(4)))
                % Paste above and uncomment to have "dark-theme" labels
                %'Color', [1 1 1], ...
                %'BackgroundColor', colors(i,:)/3, ...
        end
    end
    
    axis equal
    axis tight
    axis off
    
end

function outline(rectangles,style,thickness)
    %%  OUTLINE : Outline in black all the rectangles
    %
    %   >> outline(rectangles)
    %
    %   rectangles is a 4-by-N matrix where each column is a rectangle in
    %   the form [x0,y0,w,h]. x0 and y0 are the coordinates of the lower
    %   left-most point in the rectangle.
    %
    %   Example:
    %    cla
    %    outline([0 0 1 1; 1 0 1 1; 0 1 2 1]')
    %
    %   Copyright 2007-2013 The MathWorks, Inc. and 2022 Gabriele Bellomia
    %
    %   See also PLOTRECTANGLES. 
    
    if nargin < 2
       style = '-';     % matlab's default
    end
    
    if nargin < 3
       thickness = 0.5; % matlab's default 
    end
    
    for i = 1:size(rectangles,2)
        r = rectangles(:,i);
        xPoints = [r(1),      r(1), r(1)+r(3), r(1)+r(3)];
        yPoints = [r(2), r(2)+r(4), r(2)+r(4), r(2)     ];
        patch(xPoints,yPoints,[0 0 0],...
            'FaceColor','none','LineWidth',thickness,'LineStyle',style)
    end

end

function [c,p] = get_pure_states(file_rdm)

    asci_matrix = load(file_rdm);
    full_matrix = full(spconvert(asci_matrix(2:end,:)));
    
    [V,D] = eig(full_matrix);
    
    p = diag(D);
    
    c = cell(size(D));
    for i = 1:size(D,1)
        c{i} = V(:,i);
    end

end
    

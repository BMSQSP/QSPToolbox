function quadrantCounts = quadrantCount(x1, x2)
% This function is roughly a 2D equivalent to creating the comparison grid for the
% 1D KS test.  Here, all of the quadrant comparisons are run for each point in
% x1 as a reference against x2
%
%  Arguments
%   x1:       2xN matrix of reference points
%   x2:       2xM matrix of reference points
%
%  Returns
%   quadrantCounts

[~,nsample1] = size(x1);
[~,nsample2] = size(x2);

% - A function handle to perform comparisons in all possible directions
pointCompare = @(x, point)([((x(1,:) > point(1)) & (x(2,:) > point(2))),...
   ((x(1, :) <= point(1)) & (x(2, :) > point(2))),...
   ((x(1, :) <= point(1)) & (x(2, :) <= point(2))),...
   ((x(1, :) > point(1)) & (x(2, :) <= point(2)))]);

quadrantCounts = nan(nsample1,4*nsample2);
for j = 1:(nsample1)
	point = x1(:,j);
	quadrantCounts(j,:) = pointCompare(x2, point);
end



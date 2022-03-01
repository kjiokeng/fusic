%% Find the path that has been used for this communcation among all the resolved paths
% paths            -- the identified most dominant paths

% Return:
% used_path        -- the path that has been used
%                     It is an array in the form [aoa, tof, power]

function [used_path] = find_used_path(paths)
	% Consider that the path with the greatest power is the one that has been used
	[~, path_ind] = max(paths(:,3));
	used_path = paths(path_ind,:);
end
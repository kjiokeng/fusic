function [y, x] = cdf_helper(vals)
	% Compute the histogram of A and B.
	[counts_vals, bins_vals] = hist(vals);

	% Compute the cumulative distribution function of A and B.
	cdf_vals = cumsum(counts_vals) / sum(counts_vals);

	% Result
	x = bins_vals;
	y = cdf_vals;

	% Matlab built-in function
	[y, x] = ecdf(vals);
end
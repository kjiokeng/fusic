
ground_truth = [2 4 6 8 10 2 4 6 8 10]';
ftm_dist = [8 12.12 13.35 13.52 15.49 2.83 6.55 9.46 10.31 12.87;
			3.33 0.37 0.19 1 1.72 1.12 1.74 1.64 1.16 1]';
fusic_dist = [1.74 5.9 6.42 6.92 9.04 2.13 2.34 5.27 7.63 12.17;
			0.77 1.06 0.32 1.1 1 1.12 1.6 1.12 0.56 1]';

ftm_err = [];
fusic_err = [];

n_vals = 30;
for k=1:length(ground_truth)
	tmp = generate_random_vals(ftm_dist(k,1), ftm_dist(k,2), n_vals);
	err = abs(tmp - ground_truth(k));
	ftm_err = [ftm_err err];

	tmp = generate_random_vals(fusic_dist(k,1), fusic_dist(k,2), n_vals);
	err = abs(tmp - ground_truth(k));
	fusic_err = [fusic_err err];
end

ftm_median_err = median(ftm_err)
ftm_90th_percentile = prctile(ftm_err, 90)

fusic_median_err = median(fusic_err)
fusic_90th_percentile = prctile(fusic_err, 90)



[ftm_cdf, ftm_cdf_x] = ecdf(ftm_err);
[fusic_cdf, fusic_cdf_x] = ecdf(fusic_err);

plot(ftm_cdf_x, ftm_cdf);
hold on
plot(fusic_cdf_x, fusic_cdf);
legend("ftm", "fusic")

% Display the result
[ftm_cdf_x ftm_cdf fusic_cdf_x fusic_cdf]


function [vals] = generate_random_vals(avg, std_dev, n_vals)
	vals = avg + (rand(1, n_vals)-0.5) * std_dev;
end
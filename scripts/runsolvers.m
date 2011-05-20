function [out] = runsolvers(N, app, Dmax, Vmax, K)

	T = size(app.D, 1); 				% number of steps in a track

	tic;
	fprintf('Sampling tracks. Initializing solvers.\n');
	out.train = RandomTracks.sampleApp(N, app);
	out.test = RandomTracks.sampleApp(N, app);

	out.msol1 = MSDSolver(false, 10);
	out.msol2 = MSDSolver(false, 50);
	out.ssol = SimpleHMMSolver(Dmax, Vmax);
	out.lsol = HMMSolver(K);
	toc;

	tic;
	fprintf('MSD1.\n');
	[out.mres1, out.errs1] = out.msol1.test(out.test);
	toc

	tic;
	fprintf('MSD2.\n');
	[out.mres2, out.errs2] = out.msol2.test(out.test);
	toc

	tic;
	fprintf('fHMM.\n');
	[out.fres, out.errs3] = out.ssol.test(out.test);
	toc

	tic;
	fprintf('cHMM training.\n');
	out.LL = out.lsol.train(out.train);
	fprintf('cHMM testing.\n');
	[out.cres, out.errs4] = out.lsol.test(out.test);
	toc

	% Compare performance on first track
	m1 = out.mres1(1);
	m2 = out.mres2(1);
	f = out.fres(1);
	c = out.cres(1);

	% Diffusion
	figure;
	plot(1:T, app.D, '-+', ...
			 1:T, [m1.D, m2.D, f.D(:,1), c.D(:,1)]);
	xlabel('n (# steps)', 'FontSize', 16);
	ylabel('Diffusion Coefficient ([dist]^2/s)', 'FontSize', 16);
	labs = {'Actual', 'MSD (W=10)', 'MSD (W=50)', ...
					sprintf('fHMM (D_{max}=%d)',Dmax), ...
					sprintf('cHMM (K=[%d,%d,%d])',K)};
	legend(labs, 'FontSize', 12, 'Location', 'Best');
	title(sprintf('Diffusion Estimates for %s', app.name));

	% Velocity
	figure;
	plot(1:T, app.V(:,1), '-+', ...
			 1:T, [m1.V, m2.V, f.V(:,1), c.V(:,1)]);
	xlabel('n (# steps)', 'FontSize', 16);
	ylabel('Velocity v_x ([dist]/s)', 'FontSize', 16);
	labs{4} = sprintf('fHMM (V_{max}=%d)',Vmax);
	legend(labs, 'FontSize', 12, 'Location', 'Best');
	title(sprintf('Velocity Estimates for %s', app.name));
end

function [out] = runsolvers(app, Dmax, Vmax, K)

	N = 100;
	T = 100;

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

end

function [out] = chmmtrain(N, T, app, K, opts)
	tic;	
	fprintf('Sampling tracks. Initializing solver.\n');
	out.train = RandomTracks.sampleApp(N, app);
	out.test = RandomTracks.sampleApp(N, app);
	out.sol = HMMSolver(K, opts);
	out.Nrestarts = 2;
	out.Niters = 50;
	toc


	tic;
	fprintf('cHMM training.\n');
	out.LL = out.sol.train(out.train, out.Nrestarts, out.Niters);
	fprintf('cHMM testing.\n');
	[out.res, out.errs] = out.sol.test(out.test);
	fprintf('Mean AE error...');
	mean(out.errs.AE)
	fprintf('Mean RE error...');
	mean(out.errs.RE)
	toc
end

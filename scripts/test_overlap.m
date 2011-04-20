% OVERLAP test: how does EM's ability to distinguish between two states vary as
% the degree to which the states' emission distributions overlap increase?

N = 100;
T = 100;
S = 2;
maxIter = 10;

mc = MarkovChain([0.5; 0.5], [0.8 0.2; 0.2 0.8]);

trials = cell(50, 1);
i=1;
for mu=0.5:0.5:2.5
	for sig=0.2:0.2:2
		disp(sprintf('Trial %d', i));
		tic;
		means = [0; mu];
		stddevs = [1; sig];
		hmm = CHMM(mc, means, stddevs);
		X = hmm.sample(90, T);
		hmm0 = CHMM.random(S, 2.5, 2);
		hmm1 = hmm0.em(X, maxIter);
		
		[X, Z] = hmm.sample(10, T);
		Z1 = hmm1.infer(X);

		trials{i}.hmm = hmm;
		trials{i}.hmmEst = hmm1;
		trials{i}.Z = Z;
		trials{i}.Zest = Z1;
		i = i + 1;
		toc
	end
end

% test_em: Build a simple true model, sample a training and test set from it.
% Starting from a random initial model, try to learn the true model using EM.

%% Build a true model
S = 3;
means = [-3 0 3]';
stddevs = [0.5 1 2]';
mc = MarkovChain.random(S);
hmm = CHMM(mc, means, stddevs);

figure;
subplot(311);
hmm.plotNormals();
title('True HMM');

%% Training
Ntrain = 10;
Nrestarts = 9;
Niters = 500;
T = 100; 	% length of sequences

[X] = hmm.sample(Ntrain, T);

disp('Training with EM...');
tic;

hmm0 = CHMM.random(S, 5, 3);
hmm_orig = hmm0;
logliks = nan(Niters, Nrestarts+1);

[hmm_est, L, logliks(:,1)] = hmm0.em(X, Niters);

% do several restarts of HMM
for i=1:Nrestarts
	disp(sprintf('EM restart %d', i));
	hmm0 = CHMM.random(S, 5, 3);
	[hmm1, L1, logliks(:,i+1)] = hmm0.em(X, Niters);
	if (L1 > L)
		hmm_est = hmm1;
		hmm_orig = hmm0;
		L = L1;
	end	
end
toc

subplot(312);
hmm_orig.plotNormals();
title('Initial HMM');

subplot(313);
hmm_est.plotNormals();
title('Estimated HMM');

%% Testing
Ntest = 10;
[X, Z] = hmm.sample(Ntest, T);
disp('Testing...');
tic;
Zest = hmm_est.infer(X);
toc

disp('Computing errors...');
errs = zeros(Ntest, 1);
for i=1:Ntest
	zest = Zest(:,i); 	% states already sorted by mean
	ztru = Z(:,i);
	errs(i) = sum(zest ~= ztru);
end

figure;
plot(errs/T);
title('Percent Error for Test Sequences');

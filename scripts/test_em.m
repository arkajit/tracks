S = 3;
means = [-3 0 3];
stddevs = [0.5 1 2];
mc = MarkovChain.random(S);
hmm = CHMM(mc, means, stddevs);

%% Training
Ntrain = 15;
Nrestarts = 5;
T = 100; 	% length of sequences

[X] = hmm.sample(Ntrain, T);

disp('Training with EM...');
hmm0 = CHMM.random(S, 5, 3);
tic;
[hmm_est, L] = hmm0.em(X);
toc

% do several restarts of HMM
for i=1:Nrestarts
	disp(sprintf('EM restart %d', i));
	tic;
	hmm0 = CHMM.random(S, 5, 3);
	[hmm1, L1] = hmm0.em(X);
	if (L1 > L)
		hmm_est = hmm1;
		L = L1;
	end	
	toc
end

%% Testing
Ntest = 10;
[X, Z] = hmm.sample(Ntest, T);
disp('Testing...');
tic;
Zest = hmm_est.infer(X);
toc

disp('Computing errors...');
[~, idx] = sort(hmm_est.means);
errs = zeros(Ntest, 1);
for i=1:Ntest
	zest = idx(Zest{i}); 	% interpret states in sorted order
	ztru = Z{i};
	errs(i) = sum(zest ~= ztru);
end

plot(errs);

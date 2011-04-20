S = 3;
means = [-3 0 3]';
stddevs = [0.5 1 2]';
mc = MarkovChain.random(S);
hmm = CHMM(mc, means, stddevs);

%% Training
Ntrain = 15;
Nrestarts = 1;
T = 200; 	% length of sequences

[X] = hmm.sample(Ntrain, T);

disp('Training with EM...');
tic;
hmm0 = CHMM.random(S, 5, 3);
[hmm_est, L] = hmm0.em(X);

% do several restarts of HMM
for i=1:Nrestarts
	disp(sprintf('EM restart %d', i));
	hmm0 = CHMM.random(S, 5, 3);
	[hmm1, L1] = hmm0.em(X);
	if (L1 > L)
		hmm_est = hmm1;
		L = L1;
	end	
end
toc

%% Testing
Ntest = 10;
[X, Z] = hmm.sample(Ntest, T);
disp('Testing...');
tic;
Zest = hmm_est.infer(X);
toc

disp('Computing errors...');
[~, idx] = sort(hmm_est.means);
% invert idx function
iidx = zeros(size(idx));
for i=1:length(idx)
	iidx(idx(i)) = i;
end
errs = zeros(Ntest, 1);
for i=1:Ntest
	zest = iidx(Zest{i}); 	% interpret states in sorted order
	ztru = Z{i};
	errs(i) = sum(zest ~= ztru);
end

plot(errs);

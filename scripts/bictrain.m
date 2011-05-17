%% Train several cHMM models with different number of states, compare BIC scores
[LL, bic, train] = function bictrain(app)

N = 10;
T = 100;
n = N*T; 	% number of observations
K = 5; 		% maximum number of states

LL = zeros(K, 3);
bic = zeros(K, 3);

tic;
fprintf('Generating training data.\n');
train = RandomTracks.sampleApp(N, app);

for k=1:K
	fprintf('Training Solver %d\n', k);
	sol(k) = HMMSolver(k);
	LL(k,:) = sol(k).train(train);
	bic(k,:) = LL(k,:) - 0.5 * (k^2 + 2*k - 1) * log(n);
end

%plot results
plot(1:K, [LL/n, bic/n]);
legend('LL_x', 'LL_y', 'LL_z', 'BIC_x', 'BIC_y', 'BIC_z');
xlabel('K (# states)');
ylabel('Normalized Max. Log Likelihood (LL) and BIC');
toc

end

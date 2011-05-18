%% Train several cHMM models with different number of states, compare BIC scores
function [out] = bictrain(app)

N = 10;
T = 100;
n = N*T; 	% number of observations
K = 7; 		% maximum number of states

out.LL = zeros(K, 3);
out.bic = zeros(K, 3);
out.aic = zeros(K, 3);

tic;
fprintf('Generating training data.\n');
train = RandomTracks.sampleApp(N, app);

for k=1:K
	fprintf('Training Solver %d\n', k);
	sol(k) = HMMSolver(k);
	out.LL(k,:) = sol(k).train(train);
	deg = k^2 + 2*k - 1;
	out.bic(k,:) = out.LL(k,:) - 0.5 * deg * log(n);
	out.aic(k,:) = out.LL(k,:) - deg;
end

%plot results
figure;
plot(1:K, [out.LL/n, out.bic/n]);
legend('LL_x', 'LL_y', 'LL_z', 'BIC_x', 'BIC_y', 'BIC_z');
xlabel('K (# states)');
ylabel('Normalized Max. Log Likelihood (LL) and BIC');
toc

end

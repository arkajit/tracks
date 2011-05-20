%% Train several cHMM models with different number of states, compare BIC scores
%		Sample N tracks to measure the LL required to compute BIC.
function [out] = bictrain(N, app)

T = size(app.D, 1); 				% number of steps in a track
n = N*T; 										% number of observations
K = 7; 											% maximum number of states

Nrestarts = 2;							% number of times to restart EM
Niters = 30;								% number of times to run EM

out.LL = zeros(K, 3);				% log-likelihood scores
out.bic = zeros(K, 3);			% BIC scores
out.aic = zeros(K, 3);			% AIC scores

tic;
fprintf('Generating training data.\n');
out.train = RandomTracks.sampleApp(N, app);

for k=1:K
	fprintf('Training Solver %d\n', k);
	out.sol(k) = HMMSolver(k);
	out.LL(k,:) = out.sol(k).train(out.train, Nrestarts, Niters);
	deg = k^2 + 2*k - 1;
	out.bic(k,:) = out.LL(k,:) - 0.5 * deg * log(n);
	out.aic(k,:) = out.LL(k,:) - deg;
end

%plot results
figure;
plot(1:K, [out.LL/n, out.bic/n]);
labs = {'LL_x', 'LL_y', 'LL_z', 'BIC_x', 'BIC_y', 'BIC_z'};
legend(labs, 'FontSize', 12, 'Location', 'Best');
xlabel('K (# states)', 'FontSize', 16);
ylabel('Normalized Max. Log Likelihood (LL), BIC', 'FontSize', 16);

toc

end

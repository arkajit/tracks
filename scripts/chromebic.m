%% Train several cHMM models with different number of states, compare BIC scores
function [out] = chromebic(train)

N = length(train);
T = length(train(1).steps);
n = N*T; 	% number of observations
K = 7; 		% maximum number of states

Nrestarts = 5;
Niters = 50;

out.LL = zeros(K, 3);
out.bic = zeros(K, 3);
out.aic = zeros(K, 3);

tic;
for k=1:K
	fprintf('Training Solver %d\n', k);
	out.sol(k) = HMMSolver(k);
	out.LL(k,:) = out.sol(k).train(train, Nrestarts, Niters);
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

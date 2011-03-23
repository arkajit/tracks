function [hmm] = make_chmm(means, stddevs)

hmm.means = means;
hmm.stddevs = stddevs;
hmm.k = length(hmm.means);
if (hmm.k ~= length(hmm.stddevs)) 
	disp('Error: length mismatch');
	return;
end

hmm.t = ones(hmm.k, 1) / hmm.k;
hmm.log_t = log(hmm.t);

hmm.T = ones(hmm.k, hmm.k);
hmm.T = hmm.T ./ repmat(sum(hmm.T, 2), 1, hmm.k);
hmm.log_T = log(hmm.T);

hmm.E = @emit;
hmm.log_E = @(s, x) log(hmm.E(s, x));
hmm.logEall = @(x) log(emitAll(x));

function [prob] = emit(s, x)
	if stddevs(s)
		prob = normpdf(x, means(s), stddevs(s));
	else
		prob = double(x == means(s));
	end

function [probs] = emitAll(x)
	probs = zeros(hmm.k, 1);
	for s=1:hmm.k
		probs(s) = emit(s, x);
	end

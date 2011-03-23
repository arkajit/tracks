function [hmm] = fact_hmm(means, stddevs, t1, t2, o1, o2)
% FACT_HMM creates a Factorial HMM over two Markov Chains.
% Two MCs, one each for means and stddevs, are represented
% as a column vector of parameter values and a second
% column vector of initial probabilities
%

% o1 and o2 are odds that represent the odds against switching 
% states: i.e. odds = p/(1-p) where we
% model transitions as unlikely with a probability p
% of sticking to current state. And with probability
% (1-p) state transitions to a different state with
% equal probability
% Example: odds=10 means it is 10 times more likely
% that we'll stick to one state rather than switch.

hmm.means = means;
hmm.stddevs = stddevs;

k1 = length(means); % number of means
k2 = length(stddevs); % number of stddevs

% if initial state distributions not provided, assume uniform
if isempty(t1)
  t1 = ones(k1,1) / k1;
end;

if isempty(t2)
  t2 = ones(k2,1) / k2;
end;

t = reshape(t1 * t2', k1*k2, 1);
hmm.t = t/sum(t);
hmm.log_t = log(hmm.t);
k = length(t);
hmm.k = k;

if k1 == 1
  T1 = [1];
else
  T1 = ones(k1,k1) + ((k1-1)*o1 - 1)*eye(k1); % diag entries = (k1-1)*o1
end
hmm.T1 = T1 ./ repmat(sum(T1,2),1,k1);
log_T1 = log(hmm.T1);

if k2 == 1
  T2 = [1];
else
  T2 = ones(k2,k2) + ((k2-1)*o2 - 1)*eye(k2); % diag entries = (k2-1)*o1
end
hmm.T2 = T2 ./ repmat(sum(T2,2),1,k2);
log_T2 = log(hmm.T2);

hmm.log_T = zeros(k,k);
for s1=1:k
  for s2=1:k
    [i1, j1] = s2ind(s1);
    [i2, j2] = s2ind(s2);
    hmm.log_T(s1,s2) = log_T1(i1,i2) + log_T2(j1,j2); 
  end;
end;

hmm.T = exp(hmm.log_T);

hmm.E = @emit;
hmm.log_E = @(s,x) log(hmm.E(s,x));
hmm.logEall = @x log(emitAll(x));
hmm.s2ind = @s2ind;
hmm.ind2s = @ind2s;
hmm.to_pairs = @states2pairs;

function [i,j] = s2ind(s)
  i = mod(s, k1);
  if i == 0, i = k1; end;
  j = ceil(s / k1);
end

function [s] = ind2s(i, j)
	s = (j-1) + i;
end

function [pairs, indices] = states2pairs(states)
  n = length(states);
  pairs = zeros(n,2);
  indices = zeros(n,2);
  for k=1:n
    [i,j] = s2ind(states(k));
    pairs(k,:) = [hmm.means(i), hmm.stddevs(j)];
    indices(k,:) = [i,j];
  end
end

function [prob] = emit(s, x)
  [i, j] = s2ind(s);
  prob = normpdf(x, hmm.means(i), hmm.stddevs(j));
  if isnan(prob)
    prob = double(x == hmm.means(i));
  end
end

function [probs] = emitAll(x)
	probs = zeros(hmm.k, 1);
	for s=1:hmm.k
		probs(i) = emit(s, x);
	end	
end

end

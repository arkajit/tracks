
% ----------------------------------------------
% Viterbi assuming a continuous (Gaussian)
% emission distribution

function [state, log_v] = viterbi(x,hmm) 
    
k = size(hmm.log_t,1); % number of states
m = length(x);         % length of sequence

log_v = zeros(k,m); 
for i=1:k
  log_v(i,1) = hmm.log_t(i) + hmm.log_E(i,x(1)); 
end;
for j=2:m,
    log_V = repmat( log_v(:,j-1),1,k ); % makes a k-by-k matrix with each
                                        % col being the prev Viterbi col
    best = max(log_V+hmm.log_T,[],1)';
    for i=1:k
      log_v(i,j) = best(i) + hmm.log_E(i,x(j));
    end;
end;

% backtrack

state = zeros(1,m);
[tmp,ind] = max(log_v(:,m));
state(m) = ind(1); 
for j=m-1:-1:1,
    [tmp,ind] = max(log_v(:,j)+hmm.log_T(:,state(j+1))); 
    state(j) = ind(1); 
end;

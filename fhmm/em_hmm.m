% TODO(arkajit): update this to work with continuous Gaussian emissions

function [hmm,bic] = em_hmm(X,hmm,verbose)

if (nargin<3), verbose = 0; end; 

% hmm.log_t = log-initial state distribution
% hmm.log_T = log-trasition probabilities 
% hmm.log_E = log-emission probabilities

log_lik = 0; cont = 1; it = 1; 
while (cont), 
    log_lik_tmp = log_lik;
    [Nt,NT,NE,log_lik] = Estep(X,hmm); 
    if (verbose), fprintf('%4d %5.4f \n',it,log_lik); end;
    hmm = Mstep(Nt,NT,NE);
    
    cont = (it<5) | (log_lik-log_lik_tmp)>1e-5*abs(log_lik); 
    it = it + 1; 
end;

ndata = 0; for i=1:length(X), ndata = ndata + length(X{i}); end; %summing up the number of observation in all training instances.
[k,l] = size(hmm.E);

%counting the number of independent parameter d
%initial probalilities t has k parameters, but they have to sum to 1 => k-1 independent parameters
%transition matrix T has k x k entries, but each row has to sum to 1 => k*(k-1) independent parameters
%emission probabilities E is a k x l matrix, but each row has to sum to 1 => k*(l-1) independent parameters

bic = log_lik - (k-1+k*(k-1)+k*(l-1))/2*log(ndata); 

% ----------------------------------------------
% M-step of the EM algorithm

function [hmm] = Mstep(Nt,NT,NE)

k = size(Nt,1); % number of states
l = size(NE,2); % number of output symbol

hmm.t = Nt/sum(Nt);
hmm.log_t = log(hmm.t);

hmm.T = NT./repmat(sum(NT,2),1,k);
hmm.log_T = log(hmm.T);

hmm.E = NE./repmat(sum(NE,2),1,l);
hmm.log_E = log(hmm.E);

% ----------------------------------------------
% E-step of the EM algorithm

function [Nt,NT,NE,log_lik] = Estep(X,hmm)

% Nt = expected counts for the initial state
% NT(s,s') = expected number of transitions from s to s'
% NE(s,z) = expected number of of times we were in s and generated z
Nt = zeros(size(hmm.log_t)); 
NT = zeros(size(hmm.log_T)); 
%NE = zeros(size(hmm.log_E)); % can't take the size of a function 

k = size(hmm.log_t,1); % number of states

log_lik = 0; 
for i=1:length(X), 	% for each data example
    x = X{i}; m = length(x); 
    
    log_a = forward(x,hmm); % forward probabilities on a log-scale
    log_b = backward(x,hmm);% backward probabilities on a log-scale
    log_lik = log_lik + log_sum_exp(log_a(:,end),1); % add Pr(Data)

    % accumulate posterior counts 
    Nt = Nt + norm_exp(log_a(:,1)+log_b(:,1)); 
    for j=1:m-1,	% for each timepoint in an example
    %    NE(:,x(j)) = NE(:,x(j)) + norm_exp(log_a(:,j)+log_b(:,j)); 
    
        log_A = repmat(log_a(:,j),1,k);
        log_B = repmat(log_b(:,j+1),1,k);

				log_e = zeros(k,1);
				for s=1:k
					log_e(s) = hmm.log_E(s, x(j+1));
				end

				log_E = repmat(log_e, 1, k);
        NT = NT + norm_exp(log_A+hmm.log_T+log_E'+log_B'); 
    end;
    %NE(:,x(m)) = NE(:,x(m)) + norm_exp(log_a(:,m)+log_b(:,m));

end;

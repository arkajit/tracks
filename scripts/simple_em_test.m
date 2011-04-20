% A simple test of EM on a 1-state HMM

hmm = CHMM(MarkovChain([1], [1]), [0], [1]);
X = hmm.sample(10, 1000);

% let's start with a really bad model
hmm0 = CHMM.random(1, 5, 10);
tic; hmmEst = hmm0.em(X, 10); toc

% should be pretty close to the standard normal we sampled from
hmmEst.plotNormals()

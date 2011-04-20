% Simple test to see how well EM does when there is very little overlap between
% states in a 2-state HMM

N = 10;
T = 1000;
hmm = CHMM(MarkovChain.random(2), [0; 5], [1; 0.5]);
[X, S] = hmm.sample(N, T);

hmm0 = CHMM.random(2, 10, 5);
tic; hmmEst = hmm0.em(X, 10); toc

hmmEst.plotNormals();

tic; Sest = hmmEst.infer(X); toc

errs = HMMSolver.errors(2, S, Sest);
disp(sprintf('Error: %4.2f%%', 100*errs / (N*T)))

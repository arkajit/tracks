ta = TrackAnalyzer(30, ... % number of trials
                   20, ... % number of bins
                    5);    % number of transitions used in random trajectories

nSteps = 100;
errS = ta.analyzeLength(nSteps, 5, 1);
errTau = ta.analyzeTau(nSteps, 2, 0.1:0.1:5);
errRat = ta.analyzeRatio(nSteps, 0:0.1:5, 1);

ta = TrackAnalyzer(30, ... % number of trials
                   20, ... % number of bins
                    5);    % number of transitions used in random trajectories

errS = ta.analyzeLength(100, 5, 1);
errTau = ta.analyzeTau(100, 2, 0.1:0.1:5);
errRat = ta.analyzeRatio(100, 0:0.1:5, 1);

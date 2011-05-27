S = 100; 		% number of steps
T = S+1;		% number of positions
N = 10;			% number of tracks to sample

%% CREATE TRACKS

% From raw position data
positions = rand(T, 3);		% matrix with each column containing x,y,z positions
tau = 1;										% timestep between positions
t1 = Track.fromTrajectory(positions, tau);

% Or generate random tracks from known motion parameters
D = rand(S, 1);
V = rand(S, 3);
t2 = RandomTracks.from(D, V, tau);

% Several commonly occuring pairs of D and Vs are stored as "applications" in
% classes/TrackApps.m. Use the following to generate many tracks of the same
% application type

ta = TrackApps(S, tau);
app = ta.diffstep(3);		% a DIFFSTEP app with a step of 3
testTracks = RandomTracks.sampleApp(N, app);	% N tracks of the DIFFSTEP app
trainTracks = RandomTracks.sampleApp(N, app);	% N more tracks of the DIFFSTEP app
% can access all tracks like testTracks(1), testTracks(2), ..., testTracks(N)

% Tracks are analyzed using Solvers. All solvers implement a test method that
% can be used to solve a set of tracks at once.
% [res, errs] = sol.test(testTracks);
% res is an array of structures with the results for each input track. Each
% structure has D, V fields with inferred estimates.
% and errs is a struct with absolute (AE) and relative (RE) error fields

% HMMSolver also needs to be trained first.

%% CREATE SOLVERS

% MSD
W = 10; 	% MSD Window size
msdSol = MSDSolver(false, W);

% fHMM (SimpleHMMSolver)
Dmax = 1;
Vmax = 3;
nBins = 20;		% optional
odds = 10;		% optional
fhmmSol = SimpleHMMSolver(Dmax, Vmax); % use default nBins, odds
% or uncomment the next line to specify all the parameters
%fhmmSol = SimpleHMMSolver(Dmax, Vmax, nBins, odds);

% cHMM (HMMSolver) parameters
K = [1;2;3]; 		% states in x,y,z
hmmSol = HMMSolver(K);

%% RUN SOLVERS

% MSD
[mres, merrs] = msdSol.test(testTracks);

% fHMM
[fres, ferrs] = fhmmSol.test(testTracks);

% cHMM

% Train
[loglik] = hmmSol.train(trainTracks);
% Instead of using default number of EM restarts and iterations, can also
% specify additional parameters as below (uncomment to use).

%nRestarts = 2;
%nIterations = 15;
%[loglik] = hmmSol.train(trainTracks, nRestarts, nIterations);

% Test
[cres, cerrs] = hmmSol.test(testTracks);

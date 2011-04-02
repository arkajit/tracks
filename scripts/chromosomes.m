% TODO(arkajit): update this to run the new algorithm
% specifically, should use a train/test approach as in the test_tracks script

data = importdata('../data/chromosome_tracks.mat');

Dmax = 0.1;
Vmax = 0.2;
nBins = 10;
odds = 25;
tau = data.info.timestep;
tracknos = data.tracks(:,1);
n = max(tracknos);

trajectories = cell(n, 1);
tracks = cell(n, 1);
D = cell(n, 1);
V = cell(n, 1);
hmmSolver = HMMSolver(Dmax, Vmax, nBins, odds);

for i=1:n
	disp(sprintf('Analyzing Track %d...', i));
	savedir = sprintf('../ctfigs/%d', i);
	trajectories{i} = data.tracks(find(tracknos==i), 3:5);
	tracks{i} = Track.fromTrajectory(trajectories{i}, tau);
	[D{i}, V{i}] = hmmSolver.solve(tracks{i});
	plotter = TrackPlotter(tracks{i}, savedir);
	plotter.plotDiffusion(D{i});
	plotter.plotVelocity(V{i});	
end

save('../data/chromosome_analysis.mat', 'trajectories', 'tracks', 'D', 'V');

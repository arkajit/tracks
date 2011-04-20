TRACKS_HOME = '~/Dropbox/tracks/code';
cd(TRACKS_HOME);
data = importdata('data/chromosome_tracks.mat');

Dmax = 0.1;
Vmax = 0.2;
S = 4;
maxIter = 10;
tau = data.info.timestep;

tracknos = data.tracks(:,1);
N = max(tracknos);

trajectories = cell(N, 1);
tracks = cell(N, 1);

for i=1:N
	disp(sprintf('Creating Track %d...', i));
	trajectories{i} = data.tracks(find(tracknos==i), 3:5);
	tracks{i} = Track.fromTrajectory(trajectories{i}, tau);
end

split = floor(.8*N);
trainTracks = tracks(1:split);
%testTracks = tracks(split+1:end);

disp('Training model...');
tic;
solver = HMMSolver(Dmax, Vmax, S);
solver = solver.train(trainTracks, maxIter);
toc

disp('Testing model...');
% test on all the tracks (we need all the answers anyway...)
tic;
D = cell(N, 1);
V = cell(size(D));
[D, V] = solver.solveAll(tracks);
toc

disp('Plotting results...');
tic;
for i=1:N
	disp(sprintf('Plotting track %d...', i));
	savedir = sprintf('chromem/%d', i);
	plotter = TrackPlotter(tracks{i}, savedir);
	plotter.plotDiffusion(D{i}(:,1));
	plotter.plotVelocity(V{i});	
end
toc

save('data/chromem.mat', 'solver', 'trajectories', 'trainTracks', 'tracks', 'D', 'V');

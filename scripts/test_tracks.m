N = 100;			%		number of tracks
T = 100; 		% 	number of steps
S = 4;			% 	expected number of states (vel or diff)

% values specific for the chromosome application
Dmax = 0.2;		% microns^2/s
Vmax = 0.1;		% microns/s
tau = 15;			% s

disp('Generating data...');
tic;
tracks = cell(N, 1);
for i=1:N
	tracks{i} = RandomTracks.fixedTransitions(T, S-1, Dmax, Vmax, tau);
end

split = floor(.8*N);
trainTracks = tracks(1:split);
testTracks = tracks(split+1:end);
toc

disp('Training model...');
tic;
solver = HMMSolver(Dmax, Vmax, S);
solver = solver.train(trainTracks);
toc

disp('Testing model...');
tic;
D = cell(length(testTracks), 1);
V = cell(size(D));
[D, V] = solver.solveAll(testTracks);
toc

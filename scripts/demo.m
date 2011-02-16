nSteps = 100;
Dmax = 1;
Vmax = 5;
nBins = 20;
tau = 1;

track = RandomTracks.fromParams(nSteps, Dmax, Vmax, tau);
plotter = TrackPlotter(track);
hmmSolver = HMMSolver(Dmax, Vmax, nBins);
msdSolver = MSDSolver();
rollSolver = MSDSolver(true);

%solvers = [hmmSolver, msdSolver, rollSolver];
%TrackSolver.solveAll(track, solvers);

disp('HMM');
tic; [D1, V1] = hmmSolver.solve(track); toc
plotter.plotDiffusion(D1);
plotter.plotVelocity(V1);

disp('MSD');
tic; [D2, V2] = msdSolver.solve(track); toc
plotter.plotDiffusion(D2);
plotter.plotVelocity(V2, true);

disp('RollingMSD');
tic; [D3, V3] = rollSolver.solve(track); toc
plotter.plotDiffusion(D3);
plotter.plotVelocity(V3, true);

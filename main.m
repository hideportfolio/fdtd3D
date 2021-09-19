clc

nonref = MATERIAL(0, 0, 1000);
air = MATERIAL(340, 0, 1293);
water = MATERIAL(1500, 0, 1000);

bone = MATERIAL(4250, 1500, 2000);
%骨中縦波音速 1方向       これを4250にしたらバグる
velocity_longitudinal_bone1 = 3460;

vessel = MATERIAL(160, 481, 1072);
blood = MATERIAL(1570, 0, 1025);
purizum = MATERIAL(6000, 3670, 2500);

material = [water, bone];

% res.dx = 0.00002857/2;
% res.dy = 0.00002857/2;
% dx=305e-6;
res.dx = 0.00005;
res.dy = 0.00005;
res.dz = 0.00005;
dim = 3;
res.dt = res.dx / purizum.velocity / sqrt(dim);

% courant = 1 / (velocity_longitudinal_water * sqrt((1 / dx^2) + (1 / dy^2))); %courantの安定条件式
% dt=dx/velocity_longitudinal_bone3/sqrt(dimension);

model = SimulationSpace(res);
[space.nx, space.ny, space.nz] = size(model);
space.nt = 4000;

f = 2.0e+6;
senLen = round(25e-3 / res.dx);
inputwave = IncidenceWave(f, res, space, water, senLen);

simu = FDTD(res, space, material, model);
waveform = zeros(space.nt, 1);

for s = 1:space.nt

    tic
    waitbar(s / space.nt);
    pause(0.1);

    simu = simu.SendWave(inputwave, senLen, s);
    simu = simu.ElasticFDTD();
    simu.Image();
    waveform(s) = simu.ObservedWave();

    clc
    T = toc;
    time = T * (space.nt - s) / 3600
    s

end

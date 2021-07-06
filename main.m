clc

nonref = MATERIAL(0, 0, 1000);
air = MATERIAL(340, 0, 1293);
water = MATERIAL(1500, 0, 1000);

bone = MATERIAL(3000, 1500, 2000);
vessel = MATERIAL(160, 481, 1072);
blood = MATERIAL(1570, 0, 1025);
purizum = MATERIAL(6000, 3670, 2500);

material = [water, purizum];

% res.dx = 0.00002857/2;
% res.dy = 0.00002857/2;
res.dx = 0.00005;
res.dy = 0.00005;

res.dt = res.dx / purizum.velocity / sqrt(2);

model = SimulationSpace(res);
[space.nx, space.ny] = size(model);
space.nt = 20000;

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

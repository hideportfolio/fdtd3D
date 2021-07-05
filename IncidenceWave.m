function [inputwave] = IncidenceWave(f, res, space, water, senLen)
    %UNTITLED この関数の概要をここに記述
    %   詳細説明をここに記述
    % clc
    % water = MATERIAL(1500, 0, 1000);
    % res.dx = 0.00002857/2;
    % res.dy = 0.00002857/2;
    % res.dt = res.dx / water.velocity / sqrt(2);

    % model = SimulationSpace(res);
    % [space.nx, space.ny] = size(model);
    % space.nt = 10000;

    f = 2.0e+6;
    n = round(1 / f / res.dt);
    dis = (40 + 60/11) * 1e-3 / res.dx;

    wave(1:n) = sin(2 * pi * f * res.dt * (1:n));
    han(1:n) = 0.5 - 0.5 * cos(2 * pi * f * res.dt * (1:n));
    wave_han = wave .* han;

    X = dis;
    Y = ((1:senLen) - senLen / 2);
    arrayTimeDiff = (sqrt(X^2 + (Y .* Y)) - dis) .* res.dx ./ water.velocity ./ res.dt;
    arrayTimeDiff = abs(max(arrayTimeDiff) - arrayTimeDiff);

    inputwave = zeros(round(max(arrayTimeDiff)) + size(wave_han, 2), round(senLen));

    % inputwave(round(L(1:sl)) + 1:round(L(1:sl)) + size(wave_han, 2), 1:sl) ...
    %     = repmat(transpose(wave_han), 1, sl);
    for l = 1:round(senLen)
        inputwave(round(arrayTimeDiff(l)) + 1:round(arrayTimeDiff(l) + size(wave_han, 2)), l) = wave_han / max(wave_han);
    end

end

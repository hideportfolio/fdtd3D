function model = SimulationSpace(res)

    %     x = round(57e-3 / res.dx);
    x = round(20e-3 / res.dx);
    y = round(27e-3 / res.dx);
    z = round(27e-3 / res.dx);

    model = ones(x, y, z);
    %     [nx, ny] = size(model);

    %     puri_len = round(14.1e-3 / res.dx);
    %     purizum = ones(puri_len, puri_len);
    %     purizum(round(puri_len / 2), round(puri_len / 2)) = 2;
    %     for l = 1:round(14.1e-3 * sqrt(3) / 2 / res.dx)
    %         purizum(l, round(141/244 * l):round(puri_len - 141/244 * l + 1)) = 2;
    %     end

    % %     model(round(41e-3 / res.dx):round(41e-3 / res.dx) + puri_len - 1, round(ny / 2 - puri_len / 2):round(ny / 2 + puri_len / 2) - 1) = purizum;
    %     model(round(2e-3 / res.dx):round(2e-3 / res.dx) + puri_len - 1, round(ny / 2 - puri_len / 2):round(ny / 2 + puri_len / 2) - 1) = purizum;

end

% % bvtv=bone_data(loop_count);
% bvtv = 30;
% bone = csvread(sprintf('../bone_%d.csv', bvtv), 0, 0, [0, 0, 2350, 300]);
% % bone = csvread('bone_001.csv',0,0,[0,0,2350,300]);
% bone = imresize(bone, 2);
% bone = imrotate(bone, 90, 'bilinear');
% bone = int16(bone);
% [nx, ny] = size(bone);
% bone_point = round(10e-3 / dx);
% bone_model = zeros(round(78e-3 / dx), round(67e-3 / dy));
% % bone_model(bone_point:bone_point+round(4e-3/dx),1+5:ny+5)=1;
% bone_model(bone_point:bone_point + nx - 1, 1 + 5:ny + 5) = bone; %

% [nx, ny] = size(bone_model);

% bone_model(:, ny - 5:ny) = 0;

% inputpointx = round(8e-3 / dx);
% bone_model(inputpointx + round(60e-3 / dx):inputpointx + round(60e-3 / dx) + round(2.5e-3 / dx), :) = 2;

% bone_model(inputpointx + round(60e-3 / dx) + round(0.3e-3 / dx):inputpointx + round(60e-3 / dx) + round(2.5e-3 / dx) - round(0.3e-3 / dx), :) = 3;
% % bone_model(1:inputpointx-1,:)= 4;
% % bone_model(1003:1003+20,1:85)=3;
% % bone_model(1003:1003+20,335:50+ny-1)=3;

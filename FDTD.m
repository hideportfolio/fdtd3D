classdef FDTD
    %UNTITLED このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        res
        space
        C11
        C12
        C22
        C66
        Txx
        Tyy
        Txy
        Ux
        Uy
        rhox
        rhoy
        hig1
        hig2
        hig3
        hig4
        hig5
        hig6
        hig7
        hig8
        b
    end

    methods

        function obj = FDTD(res, space, material, model)
            c11 = zeros(1, length(material)); %map関数
            c66 = zeros(1, length(material));
            density = zeros(1, length(material));

            for m = 1:length(material)
                c11(m) = material(m).density * material(m).velocity^2;
                c66(m) = material(m).density * material(m).side_velocity^2;
                density(m) = material(m).density;
            end

            c12 = c11 - 2 * c66;

            obj.space = space;
            obj.res = res;
            obj.C11 = c11(model(:, :));
            obj.C22 = obj.C11;
            obj.C12 = c12(model(:, :));
            obj.C66 = (c66(model(1:space.nx - 1, 1:space.ny - 1)) + c66(model(1:space.nx - 1, 2:space.ny)) + c66(model(2:space.nx, 1:space.ny - 1)) + c66(model(2:space.nx, 2:space.ny))) / 4;

            obj.Txx(:, :) = zeros(space.nx, space.ny);
            obj.Tyy(:, :) = zeros(space.nx, space.ny);
            obj.Txy(:, :) = zeros(space.nx + 1, space.ny + 1);
            obj.Ux(:, :) = zeros(space.nx + 1, space.ny);
            obj.Uy(:, :) = zeros(space.nx, space.ny + 1);
            obj.rhox = (density(model(2:space.nx, :)) + density(model(1:space.nx - 1, :))) / 2;
            obj.rhoy = (density(model(:, 2:space.ny)) + density(model(:, 1:space.ny - 1))) / 2;

            obj.hig1 = zeros(6, space.ny);
            obj.hig2 = zeros(6, space.ny);
            obj.hig3 = zeros(space.nx, 6);
            obj.hig4 = zeros(space.nx, 6);
            obj.hig5 = zeros(6, space.ny);
            obj.hig6 = zeros(6, space.ny);
            obj.hig7 = zeros(space.nx, 6);
            obj.hig8 = zeros(space.nx, 6);

            obj.b = (material(1).velocity * res.dt - res.dx) / (material(1).velocity * res.dt + res.dx);
        end

        function obj = SendWave(obj, inputwave, senLen, s)
            posi.x = round(1e-3 / obj.res.dx);
            posi.y = round(obj.space.ny / 2);

            if s <= size(inputwave, 1)
                % obj.Txx(posi.x, posi.y - round(senLen / 2) + 1:posi.y + round(senLen / 2) - 1) = 1;
                % obj.Tyy(posi.x, posi.y - round(senLen / 2) + 1:posi.y + round(senLen / 2) - 1) = 1;
                obj.Txx(posi.x, posi.y - round(senLen / 2) + 1:posi.y + round(senLen / 2)) = inputwave(s, 1:senLen);
                obj.Tyy(posi.x, posi.y - round(senLen / 2) + 1:posi.y + round(senLen / 2)) = inputwave(s, 1:senLen);
                % else
                %     obj.Txx(posi.x, posi.y - round(senLen / 2):posi.y + round(senLen / 2)) = 0;
                %     obj.Tyy(posi.x, posi.y - round(senLen / 2):posi.y + round(senLen / 2)) = 0;
            end

        end

        function obj = ElasticFDTD(obj)
            a = 1;
            b1 = obj.b;
            b2 = b1;
            b3 = b1;
            b4 = b1;
            d1 = 0.00001;
            d2 = d1;
            D = obj.res.dt / obj.res.dx;
            nx = obj.space.nx;
            ny = obj.space.ny;

            obj.Ux(2:nx, :) = obj.Ux(2:nx, :) + D ./ obj.rhox .* ((obj.Txx(2:nx, :) - obj.Txx(1:nx - 1, :)) + (obj.Txy(2:nx, 2:ny + 1) - obj.Txy(2:nx, 1:ny)));
            obj.Uy(:, 2:ny) = obj.Uy(:, 2:ny) + D ./ obj.rhoy .* ((obj.Txy(2:nx + 1, 2:ny) - obj.Txy(1:nx, 2:ny)) + (obj.Tyy(:, 2:ny) - obj.Tyy(:, 1:ny - 1)));
            obj.Txx(:, :) = obj.Txx(:, :) + obj.C11 * D .* (obj.Ux(2:nx + 1, :) - obj.Ux(1:nx, :)) + obj.C12 .* D .* (obj.Uy(:, 2:ny + 1) - obj.Uy(:, 1:ny));
            obj.Tyy(:, :) = obj.Tyy(:, :) + obj.C12 * D .* (obj.Ux(2:nx + 1, :) - obj.Ux(1:nx, :)) + obj.C22 .* D .* (obj.Uy(:, 2:ny + 1) - obj.Uy(:, 1:ny));
            obj.Txy(2:nx, 2:ny) = obj.Txy(2:nx, 2:ny) + obj.C66 * D .* ((obj.Uy(2:nx, 2:ny) - obj.Uy(1:nx - 1, 2:ny)) + (obj.Ux(2:nx, 2:ny) - obj.Ux(2:nx, 1:ny - 1)));

            obj.Txx(1, :) = (b1 + b2) * (obj.Txx(2, :) - obj.hig1(1, :)) - b1 * b2 * (obj.Txx(3, :) - 2 * obj.hig1(2, :) + obj.hig1(4, :)) - (b1 * (1 - d2) + b2 * (1 - d1)) * (obj.hig1(3, :) - obj.hig1(5, :)) + ((1 - d1) + (1 - d2)) * obj.hig1(2, :) - (1 - d1) * (1 - d2) * obj.hig1(6, :); %��
            obj.Txx(obj.space.nx, :) = (b1 + b2) * (obj.Txx(obj.space.nx - 1, :) - obj.hig2(3, :)) - b1 * b2 * (obj.Txx(obj.space.nx - 2, :) - 2 * obj.hig2(2, :) + obj.hig2(6, :)) - (b1 * (1 - d2) + b2 * (1 - d1)) * (obj.hig2(1, :) - obj.hig2(5, :)) + ((1 - d1) + (1 - d2)) * obj.hig2(2, :) - (1 - d1) * (1 - d2) * obj.hig2(4, :); %��
            obj.Txx(:, 1) = (b3 + b4) * (obj.Txx(:, 2) - obj.hig3(:, 1)) - b3 * b4 * (obj.Txx(:, 2) - 2 * obj.hig3(:, 2) + obj.hig3(:, 4)) - (b3 * (1 - d2) + b4 * (1 - d1)) * (obj.hig3(:, 3) - obj.hig3(:, 5)) + ((1 - d1) + (1 - d2)) * obj.hig3(:, 2) - (1 - d1) * (1 - d2) * obj.hig3(:, 6); %��
            obj.Txx(:, obj.space.ny) = (b3 + b4) * (obj.Txx(:, obj.space.ny - 1) - obj.hig4(:, 3)) - b3 * b4 * (obj.Txx(:, obj.space.ny - 2) - 2 * obj.hig4(:, 2) + obj.hig4(:, 6)) - (b3 * (1 - d2) + b4 * (1 - d1)) * (obj.hig4(:, 1) - obj.hig4(:, 5)) + ((1 - d1) + (1 - d2)) * obj.hig4(:, 2) - (1 - d1) * (1 - d2) * obj.hig4(:, 4);
            obj.Tyy(1, :) = (b1 + b2) * (obj.Tyy(2, :) - obj.hig5(1, :)) - b1 * b2 * (obj.Tyy(3, :) - 2 * obj.hig5(2, :) + obj.hig5(4, :)) - (b1 * (1 - d2) + b2 * (1 - d1)) * (obj.hig5(3, :) - obj.hig5(5, :)) + ((1 - d1) + (1 - d2)) * obj.hig5(2, :) - (1 - d1) * (1 - d2) * obj.hig5(6, :); %��
            obj.Tyy(obj.space.nx, :) = (b1 + b2) * (obj.Tyy(obj.space.nx - 1, :) - obj.hig6(3, :)) - b1 * b2 * (obj.Tyy(obj.space.nx - 2, :) - 2 * obj.hig6(2, :) + obj.hig6(6, :)) - (b1 * (1 - d2) + b2 * (1 - d1)) * (obj.hig6(1, :) - obj.hig6(5, :)) + ((1 - d1) + (1 - d2)) * obj.hig6(2, :) - (1 - d1) * (1 - d2) * obj.hig6(4, :); %��
            obj.Tyy(:, 1) = (b3 + b4) * (obj.Tyy(:, 2) - obj.hig7(:, 1)) - b3 * b4 * (obj.Tyy(:, 3) - 2 * obj.hig7(:, 2) + obj.hig7(:, 4)) - (b3 * (1 - d2) + b4 * (1 - d1)) * (obj.hig7(:, 3) - obj.hig7(:, 5)) + ((1 - d1) + (1 - d2)) * obj.hig7(:, 2) - (1 - d1) * (1 - d2) * obj.hig7(:, 6); %��
            obj.Tyy(:, obj.space.ny) = (b3 + b4) * (obj.Tyy(:, obj.space.ny - 1) - obj.hig8(:, 3)) - b3 * b4 * (obj.Tyy(:, obj.space.ny - 2) - 2 * obj.hig8(:, 2) + obj.hig8(:, 6)) - (b3 * (1 - d2) + b4 * (1 - d1)) * (obj.hig8(:, 1) - obj.hig8(:, 5)) + ((1 - d1) + (1 - d2)) * obj.hig8(:, 2) - (1 - d1) * (1 - d2) * obj.hig8(:, 4); %�E

            obj.hig1(4:6, :) = obj.hig1(1:3, :);
            obj.hig1(1:3, :) = obj.Txx(1:3, :);
            obj.hig2(4:6, :) = obj.hig2(1:3, :);
            obj.hig2(1:3, :) = obj.Txx(obj.space.nx - 2:obj.space.nx, :);
            obj.hig3(:, 4:6) = obj.hig3(:, 1:3);
            obj.hig3(:, 1:3) = obj.Txx(:, 1:3);
            obj.hig4(:, 4:6) = obj.hig4(:, 1:3);
            obj.hig4(:, 1:3) = obj.Txx(:, obj.space.ny - 2:obj.space.ny);

            obj.hig5(4:6, :) = obj.hig5(1:3, :);
            obj.hig5(1:3, :) = obj.Tyy(1:3, :);
            obj.hig6(4:6, :) = obj.hig6(1:3, :);
            obj.hig6(1:3, :) = obj.Tyy(obj.space.nx - 2:obj.space.nx, :);
            obj.hig7(:, 4:6) = obj.hig7(:, 1:3);
            obj.hig7(:, 1:3) = obj.Tyy(:, 1:3);
            obj.hig8(:, 4:6) = obj.hig8(:, 1:3);
            obj.hig8(:, 1:3) = obj.Tyy(:, obj.space.ny - 2:obj.space.ny);
        end

        function Image(obj)
            % time(s) = s * obj.res.dt;
            T = sqrt(obj.Txx.^2 + obj.Tyy.^2);
            imagesc(T);
            axis equal
            caxis([0 0.1]);
            M = getframe;
        end

        function waveform = ObservedWave(obj)
            center = round(obj.space.ny / 2);
            sen = 3e-3 / obj.res.dx;
            posix = round(41e-3 / obj.res.dx) - 1;
            waveform = sum(obj.Txx(posix, center - round(sen / 2):center + round(sen / 2)))

            % focalx = inputpointx + round(60e-3 / dx);
            % focaly = 2637;
            % if mod(s, 100) == 0
            %     csvwrite(sprintf('D:/movie/movie_%04d.csv', s), T);
            % end
            %             posi.x = 1e-3
            %             waveform=obj.Txx(posi.x,)

            % sen_length = round(0.05e-3 / obj.res.dy);
            % senLen = round(1e-3 / res.dx);
            % receivecount = round(senLen / sen_length);
            % receive_start = posi.y - round(sen / 2);
            % waveform1 = zeros(space.nt, receivecount + 1);
            % res_senLen = 1e-3 / obj.res.dt;

            % for num = 1:receivecount - 1

            %     obj.waveform1(s, num) = obj.waveform1(s, num) + sum(obj.Txx(posi.x, receive_start + nn + sen_length * (num - 1):posi.x, receive_start + nn + sen_length * (num - 1) + res_senLen)) / res_senLen;

            % end

            % waveform2(s) = sum(Txx(inputpointx + 10, inputpointy - round(senLen / 2):inputpointy + round(senLen / 2))) / round(senLen / 2);
        end

    end

end

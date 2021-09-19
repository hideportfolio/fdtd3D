classdef FDTD

    properties
        res
        space
        C11
        C22
        C33
        C12
        C13

        C44
        C55
        C66
        Txx
        Tyy
        Tzz

        Txy
        Tyz
        Tzx

        Ux
        Uy
        Uz

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
        hig9
        hig10
        hig11
        hig12
        hig13
        hig14
        hig15
        hig16
        hig17
        hig18
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

            %%%%%  密度パラメータ
            % density_x=(Model3d(2:nx,:,:)+Model3d(1:nx-1,:,:))/2;  %密度パラメータ
            % density_y=(Model3d(:,2:ny,:)+Model3d(:,1:ny-1,:))/2;
            % density_z=(Model3d(:,:,2:nz)+Model3d(:,:,1:nz-1))/2;
            % %%%%%%%%%%%%     密度情報を入れていく
            % for n=1:nx-1
            %     for m=1:ny
            %         for l=1:nz
            %             if density_x(n,m,l)==0
            %                 density_x(n,m,l)=rhow;

            %             else
            %                 density_x(n,m,l)=rhob;

            %             end
            %         end
            %     end
            % end
            %%%%%%%%%%%%%%弾性定数 C13,C23 の設定
            C13 = zeros(nx, ny, nz);
            ratio_elastic_C12C13 = 1.0212; %C12の1.0212倍がC13

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

            obj.hig1 = zeros(6, ny, nz); %上端 Txx
            obj.hig2 = zeros(6, ny, nz); %下端 Txx
            obj.hig3 = zeros(nx, 6, nz); %左端 Txx
            obj.hig4 = zeros(nx, 6, nz); %右端 Txx
            obj.hig5 = zeros(nx, ny, 6); %奥端 Txx
            obj.hig6 = zeros(nx, ny, 6); %手前端 Txx

            obj.hig7 = zeros(6, ny, nz); %上端 Tyy
            obj.hig8 = zeros(6, ny, nz); %下端 Tyy
            obj.hig9 = zeros(nx, 6, nz); %左端 Tyy
            obj.hig10 = zeros(nx, 6, nz); %右端 Tyy
            obj.hig11 = zeros(nx, ny, 6); %奥端 Tyy
            obj.hig12 = zeros(nx, ny, 6); %手前端 Tyy

            obj.hig13 = zeros(6, ny, nz); %上端 Tzz
            obj.hig14 = zeros(6, ny, nz); %下端 Tzz
            obj.hig15 = zeros(nx, 6, nz); %左端 Tzz
            obj.hig16 = zeros(nx, 6, nz); %右端 Tzz
            obj.hig17 = zeros(nx, ny, 6); %奥端 Tzz
            obj.hig18 = zeros(nx, ny, 6); %手前端 Tzz

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

            %FDTDの式 
            Ux(2:nx,:,:)=Ux(2:nx,:,:)+D/density_x.*((Txx(2:nx,:,:)-Txx(1:nx-1,:,:))+...
                                                        (Txy(2:nx,2:ny+1,1:nz)-Txy(2:nx,1:ny,1:nz))+...
                                                        (Tzx(2:nx,1:ny,2:nz+1)-Tzx(2:nx,1:ny,1:nz)));

            Uy(:,2:ny,:)=Uy(:,2:ny,:)+D/density_y.*((Txy(2:nx+1,2:ny,1:nz)-Txy(1:nx,2:ny,1:nz))+...
                                                        (Tyy(:,2:ny,:)-Tyy(:,1:ny-1,:))+...
                                                        (Tyz(1:nx,2:ny,2:nz+1)-Tyz(1:nx,2:ny,1:nz)));

            Uz(:,:,2:nz)=Uz(:,:,2:nz)+D/density_z.*((Tzx(2:nx+1,1:ny,2:nz)-Tzx(1:nx,1:ny,2:nz))+...
                                                        (Tyz(1:nx,2:ny+1,2:nz)-Tyz(1:nx,1:ny,2:nz))...
                                                        +(Tzz(:,:,2:nz)-Tzz(:,:,1:nz-1)));


            Txx(:,:,:)=Txx(:,:,:)+C11*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                                C12*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                                C13*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

            Tyy(:,:,:)=Tyy(:,:,:)+C12*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                                C22*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                                C23*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

            Tzz(:,:,:)=Tzz(:,:,:)+C13*X.*(Ux(2:nx+1,:,:)-Ux(1:nx,:,:))+...
                                                C23*Y.*(Uy(:,2:ny+1,:)-Uy(:,1:ny,:))+...
                                                C33*Z.*(Uz(:,:,2:nz+1)-Uz(:,:,1:nz));

            Tyz(2:nx,2:ny,2:nz)=Tyz(2:nx,2:ny,2:nz)+C44*D.*((Uz(2:nx,2:ny,2:nz)-Uz(2:nx,1:ny-1,2:nz))+...
                (Uy(2:nx,2:ny,2:nz)-Uy(2:nx,2:ny,1:nz-1)));

            Tzx(2:nx,2:ny,2:nz)=Tzx(2:nx,2:ny,2:nz)+C55*D.*((Uz(2:nx,2:ny,2:nz)-Uz(1:nx-1,2:ny,2:nz))+...
                (Ux(2:nx,2:ny,2:nz)-Ux(2:nx,2:ny,1:nz-1)));

            Txy(2:nx,2:ny,2:nz)=Txy(2:nx,2:ny,2:nz)+C66*D.*((Ux(2:nx,2:ny,2:nz)-Ux(2:nx,1:ny-1,2:nz))+...    
                (Uy(2:nx,2:ny,2:nz)-Uy(1:nx-1,2:ny,2:nz)));


        %吸収境界条件式
            Txx(:,ny,:)=(b3+b4)*(Txx(:,ny-1,:)-hig4(:,3,:))-b3*b4*(Txx(:,ny-2,:)-2*hig4(:,2,:)+hig4(:,6,:))...
                -(b3*(1-d2)+b4*(1-d1))*(hig4(:,1,:)-hig4(:,5,:))+((1-d1)+(1-d2))*hig4(:,2,:)...
                -(1-d1)*(1-d2)*hig4(:,4,:); %右端
            Tyy(:,ny,:)=(b1+b2)*(Tyy(:,ny-1,:)-hig10(:,3,:))-b1*b2*(Tyy(:,ny-2,:)-2*hig10(:,2,:)+hig10(:,6,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig10(:,1,:)-hig10(:,5,:))+((1-d1)+(1-d2))*hig10(:,2,:)...
                -(1-d1)*(1-d2)*hig10(:,4,:); %右端
            Tzz(:,ny,:)=(b1+b2)*(Tzz(:,ny-1,:)-hig16(:,3,:))-b1*b2*(Tzz(:,ny-2,:)-2*hig16(:,2,:)+hig16(:,6,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig16(:,1,:)-hig16(:,5,:))+((1-d1)+(1-d2))*hig16(:,2,:)...
                -(1-d1)*(1-d2)*hig16(:,4,:); %右端
            
            Txx(:,1,:)=(b3+b4)*(Txx(:,2,:)-hig3(:,1,:))-b3*b4*(Txx(:,3,:)-2*hig3(:,2,:)+hig3(:,4,:))...
                -(b3*(1-d2)+b4*(1-d1))*(hig3(:,3,:)-hig3(:,5,:))+((1-d1)+(1-d2))*hig3(:,2,:)...
                -(1-d1)*(1-d2)*hig3(:,6,:); %左端
            Tyy(:,1,:)=(b1+b2)*(Tyy(:,2,:)-hig9(:,1,:))-b1*b2*(Tyy(:,3,:)-2*hig9(:,2,:)+hig9(:,4,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig9(:,3,:)-hig9(:,5,:))+((1-d1)+(1-d2))*hig9(:,2,:)...
                -(1-d1)*(1-d2)*hig9(:,6,:); %左端
            Tzz(:,1,:)=(b1+b2)*(Tzz(:,2,:)-hig15(:,1,:))-b1*b2*(Tzz(:,3,:)-2*hig15(:,2,:)+hig15(:,4,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig15(:,3,:)-hig15(:,5,:))+((1-d1)+(1-d2))*hig15(:,2,:)...
                -(1-d1)*(1-d2)*hig15(:,6,:); %左端
            
            Txx(1,:,:)=(b1+b2)*(Txx(2,:,:)-hig1(1,:,:))-b1*b2*(Txx(3,:,:)-2*hig1(2,:,:)+hig1(4,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig1(3,:,:)-hig1(5,:,:))+((1-d1)+(1-d2))*hig1(2,:,:)...
                -(1-d1)*(1-d2)*hig1(6,:,:); %上端
            Tyy(1,:,:)=(b1+b2)*(Tyy(2,:,:)-hig7(1,:,:))-b1*b2*(Tyy(3,:,:)-2*hig7(2,:,:)+hig7(4,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig7(3,:,:)-hig7(5,:,:))+((1-d1)+(1-d2))*hig7(2,:,:)...
                -(1-d1)*(1-d2)*hig7(6,:,:); %上端
            Tzz(1,:,:)=(b1+b2)*(Tzz(2,:,:)-hig13(1,:,:))-b1*b2*(Tzz(3,:,:)-2*hig13(2,:,:)+hig13(4,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig13(3,:,:)-hig13(5,:,:))+((1-d1)+(1-d2))*hig13(2,:,:)...
                -(1-d1)*(1-d2)*hig13(6,:,:); %上端
            
            Txx(nx,:,:)=(b1+b2)*(Txx(nx-1,:,:)-hig2(3,:,:))-b1*b2*(Txx(nx-2,:,:)-2*hig2(2,:,:)+hig2(6,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig2(1,:,:)-hig2(5,:,:))+((1-d1)+(1-d2))*hig2(2,:,:)...
                -(1-d1)*(1-d2)*hig2(4,:,:); %下端
            Tyy(nx,:,:)=(b1+b2)*(Tyy(nx-1,:,:)-hig8(3,:,:))-b1*b2*(Tyy(nx-2,:,:)-2*hig8(2,:,:)+hig8(6,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig8(1,:,:)-hig8(5,:,:))+((1-d1)+(1-d2))*hig8(2,:,:)...
                -(1-d1)*(1-d2)*hig8(4,:,:); %下端
            Tzz(nx,:,:)=(b1+b2)*(Tzz(nx-1,:,:)-hig14(3,:,:))-b1*b2*(Tzz(nx-2,:,:)-2*hig14(2,:,:)+hig14(6,:,:))...
                -(b1*(1-d2)+b2*(1-d1))*(hig14(1,:,:)-hig14(5,:,:))+((1-d1)+(1-d2))*hig14(2,:,:)...
                -(1-d1)*(1-d2)*hig14(4,:,:); %下端

            Txx(:,:,1)=(b1+b2)*(Txx(:,:,2)-hig6(:,:,1))-b1*b2*(Txx(:,:,3)-2*hig6(:,:,2)+hig6(:,:,4))...
                -(b1*(1-d2)+b2*(1-d1))*(hig6(:,:,3)-hig6(:,:,5))+((1-d1)+(1-d2))*hig6(:,:,2)...
                -(1-d1)*(1-d2)*hig6(:,:,6); %手前
            Tyy(1:nx,1:ny,1)=(b1+b2)*(Tyy(1:nx,1:ny,2)-hig12(:,:,1))-b1*b2*(Tyy(1:nx,1:ny,3)-2*hig12(:,:,2)+hig12(:,:,4))...
                -(b1*(1-d2)+b2*(1-d1))*(hig12(:,:,3)-hig12(:,:,5))+((1-d1)+(1-d2))*hig12(:,:,2)...
                -(1-d1)*(1-d2)*hig12(:,:,6); %手前
            Tzz(1:nx,1:ny,1)=(b1+b2)*(Tzz(1:nx,1:ny,2)-hig18(:,:,1))-b1*b2*(Tzz(1:nx,1:ny,3)-2*hig18(:,:,2)+hig18(:,:,4))...
                -(b1*(1-d2)+b2*(1-d1))*(hig18(:,:,3)-hig18(:,:,5))+((1-d1)+(1-d2))*hig18(:,:,2)...
                -(1-d1)*(1-d2)*hig18(:,:,6); %手前
            
            Txx(:,:,nz)=(b1+b2)*(Txx(:,:,nz-1)-hig5(:,:,3))-b1*b2*(Txx(:,:,nz-2)-2*hig5(:,:,2)+hig5(:,:,6))...
                -(b1*(1-d2)+b2*(1-d1))*(hig5(:,:,1)-hig5(:,:,5))+((1-d1)+(1-d2))*hig5(:,:,2)...
                -(1-d1)*(1-d2)*hig5(:,:,4); %奥
            Tyy(:,:,nz)=(b1+b2)*(Tyy(:,:,nz-1)-hig11(:,:,3))-b1*b2*(Tyy(:,:,nz-2)-2*hig11(:,:,2)+hig11(:,:,6))...
                -(b1*(1-d2)+b2*(1-d1))*(hig11(:,:,1)-hig11(:,:,5))+((1-d1)+(1-d2))*hig11(:,:,2)...
                -(1-d1)*(1-d2)*hig11(:,:,4); %奥
            Tzz(:,:,nz)=(b1+b2)*(Tzz(:,:,nz-1)-hig17(:,:,3))-b1*b2*(Tzz(:,:,nz-2)-2*hig17(:,:,2)+hig17(:,:,6))...
                -(b1*(1-d2)+b2*(1-d1))*(hig17(:,:,1)-hig17(:,:,5))+((1-d1)+(1-d2))*hig17(:,:,2)...
                -(1-d1)*(1-d2)*hig17(:,:,4); %奥
            

            
            
            hig1(4:6,:,:)=hig1(1:3,:,:); %2秒前上端
            hig1(1:3,:,:)=Txx(1:3,:,:);  %1秒前上端
            
            hig2(4:6,:,:)=hig2(1:3,:,:); %2秒前下端
            hig2(1:3,:,:)=Txx(nx-2:nx,:,:);  %1秒前下端
            
            hig3(:,4:6,:)=hig3(:,1:3,:); %2秒前左端
            hig3(:,1:3,:)=Txx(:,1:3,:);  %1秒前左端
            
            hig4(:,4:6,:)=hig4(:,1:3,:); %2秒前右端
            hig4(:,1:3,:)=Txx(:,ny-2:ny,:);  %1秒前右端
            
            hig6(:,:,4:6)=hig6(:,:,1:3); %2秒前手前
            hig6(:,:,1:3)=Txx(:,:,1:3);  %1秒前手前
            
            hig5(:,:,4:6)=hig5(:,:,1:3); %2秒前奥
            hig5(:,:,1:3)=Txx(:,:,nz-2:nz);  %1秒前奥

            hig7(4:6,:,:)=hig7(1:3,:,:); %2秒前上端
            hig7(1:3,:,:)=Tyy(1:3,:,:);  %1秒前上端
            
            hig8(4:6,:,:)=hig8(1:3,:,:); %2秒前下端
            hig8(1:3,:,:)=Tyy(nx-2:nx,:,:);  %1秒前下端
            
            hig9(:,4:6,:)=hig9(:,1:3,:); %2秒前左端
            hig9(:,1:3,:)=Tyy(:,1:3,:);  %1秒前左端
            
            hig10(:,4:6,:)=hig10(:,1:3,:); %2秒前右端
            hig10(:,1:3,:)=Tyy(:,ny-2:ny,:);  %1秒前右端
            
            hig12(:,:,4:6)=hig12(:,:,1:3); %2秒前手前
            hig12(:,:,1:3)=Tyy(:,:,1:3);  %1秒前手前
            
            hig11(:,:,4:6)=hig11(:,:,1:3); %2秒前奥
            hig11(:,:,1:3)=Tyy(:,:,nz-2:nz);  %1秒前奥

            hig13(4:6,:,:)=hig13(1:3,:,:); %2秒前上端
            hig13(1:3,:,:)=Tzz(1:3,:,:);  %1秒前上端
            
            hig14(4:6,:,:)=hig14(1:3,:,:); %2秒前下端
            hig14(1:3,:,:)=Tzz(nx-2:nx,:,:);  %1秒前下端
            
            hig15(:,4:6,:)=hig15(:,1:3,:); %2秒前左端
            hig15(:,1:3,:)=Tzz(:,1:3,:);  %1秒前左端
            
            hig16(:,4:6,:)=hig16(:,1:3,:); %2秒前右端
            hig16(:,1:3,:)=Tzz(:,ny-2:ny,:);  %1秒前右端
            
            hig18(:,:,4:6)=hig18(:,:,1:3); %2秒前手前
            hig18(:,:,1:3)=Tzz(:,:,1:3);  %1秒前手前
            
            hig17(:,:,4:6)=hig17(:,:,1:3); %2秒前奥
            hig17(:,:,1:3)=Tzz(:,:,nz-2:nz);  %1秒前奥
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
            posix = round(2e-3 / obj.res.dx) - 1;
            waveform = sum(obj.Txx(posix, center - round(sen / 2):center + round(sen / 2)));

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

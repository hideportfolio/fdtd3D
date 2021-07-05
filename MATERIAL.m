classdef MATERIAL
    %MATERIAL このクラスの概要をここに記述
    %   詳細説明をここに記述

    properties
        velocity
        side_velocity
        density
    end

    methods

        function obj = MATERIAL(velocity, side_velocity, density)
            %MATERIAL このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj.velocity = velocity;
            obj.side_velocity = side_velocity;
            obj.density = density;
        end

    end

end

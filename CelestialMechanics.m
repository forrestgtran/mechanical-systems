% Now with animation!
function res = celestial_mechanics_1 ()
%% initialize values 
% Star 1
m1  =  2e30                    ;%mass of Sun (kg)
P1i =  [0, 0]                  ;%star location vector
V1i =  [0, 0]                  ;%initial orbital velocity of star 2 (m/s)

% Body 2
m2  =  1.8986e27               ;%mass of body 2 (kg)
P2i =  [778500000000, 0]       ;%planet location vector
V2i =  [0, 13100]              ;%initial orbital velocity of star 2 (m/s)

initial_values = [P1i, V1i, P2i,V2i];

%% universe variables
G = 6.67e-11 ; %Nm^2/kg^2
tdays = 4332;
tseconds = tdays * 24 *60 * 60;
step = tseconds/200;
time_span = (0 : step : tseconds);

%% ode 45
%options = odeset ('Events', @events);
[T, Output] = ode45 (@jupitergoesvroom, time_span, initial_values );

%% plot
% hold on
% plot (Output(:,1), Output(:,2), 'y', 'linewidth', 4)
% plot (Output(:,5), Output(:,6), 'r', 'linewidth', 4)

animate_func (T, Output)

%% functions
    function res = jupitergoesvroom (~, W)
        %unpack positions
        P1 = [W(1), W(2)];
        
        P2 = [W(5), W(6)];
        
        %unpack velocities
        V1x = W(3);
        V1y = W(4);
        
        V2x = W(7);
        V2y = W(8);
        
        %Pack Velocitties
        dP1dt = [V1x;
                V1y];
        
        dP2dt = [V2x;
                V2y];
            
        %Pack Accelerations    
        dV = acceleration (P1, P2);
        dV1dt = dV (1:2);
        dV2dt = dV (3:4);
        
        res = [dP1dt;
            dV1dt;
            dP2dt;
            dV2dt];
    end
    function res = acceleration (P1, P2)
        
        %gravitational_force = gravity_force_func (m1, P1, m2, P2);
        f1 =  F12(m1, P1, m2, P2);
        fx1 = f1 (1);
        fy1 = f1 (2);
        
        fx2 = fx1;
        fy2 = fy1;
        
        xacceleration1 = fx1/m1;
        yacceleration1 = fy1/m1;
        
        xacceleration2 = -fx2/m2;
        yacceleration2 = -fy2/m2;
        
        res =   [xacceleration1;
            yacceleration1;
            xacceleration2;
            yacceleration2];
    end

 
    function res = F12 (m1, P1, m2, P2)
        % define r terms
        R = P2 - P1;
        r = norm (R);
        rhat = R/r;
        
        % force equation
        fg = G * (m1 * m2)/ (r^2) * rhat;
        res = fg;
    end


%% animation
    function animate_func(T,M)
        
        % animate the positions of the planets, assuming that the
        % columns of M are x1, y1, x2, y2.
        X1 = M(:,1);
        Y1 = M(:,2);
        X2 = M(:,5);
        Y2 = M(:,6);
        
        minmax = [min([X1;X2]), max([X1;X2]), min([Y1;Y2]), max([Y1;Y2])];
        
        for i=1:length(T)
            clf;
            axis(minmax);
            hold on;
            draw_func(X1(i), Y1(i), X2(i), Y2(i));
            drawnow;
        end
        
        function draw_func(x1, y1, x2, y2)
            plot(x1, y1, 'r.', 'MarkerSize', 50);
            plot(x2, y2, 'b.', 'MarkerSize', 20);
        end
    end

end

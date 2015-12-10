% Now with animation!
function res = rocket1(vx, vy)
%% initialize values
%up to 64 bodies
Number_of_Bodies = 5;

%% User - Define Initial Values for each Body

%Matrix of 4 wide and number of bodies down
%initialize values
masses = zeros(Number_of_Bodies, 1);             % 1 value for each star
initial_values = zeros(4 * Number_of_Bodies, 1); % 4 values for each body

P = 12e10; 
V = 15e4;

%sun
m0  =  5e31                       ;%mass of Sun (kg)
P0i =  [P, 0]                       ;%star location vector
V0i =  [0, -V]                       ;%initial orbital velocity of star 2 (m/s)

% mercury
m1  =  m0                       ;%mass of mercury (kg)
P1i =  [-P, 0]                ;%star location vector
V1i =  [0, V]                ;%initial orbital velocity of mercury (m/s)

% venus
m2  =  m0                     ;%mass of body venus (kg)
P2i =  [0, -P]               ;%planet location vector
V2i =  [-V, 0]                   ;%initial orbital velocity of venus (m/s)

% rocket
m3  =  2e6                     ;%mass of body venus (kg)
P3i =  [-9e11, 0]               ;%planet location vector
V3i =  [vx, vy];%[vx, vy]                   ;%initial orbital velocity of venus (m/s)

% venus
m4  =  m0                     ;%mass of body venus (kg)
P4i =  [0, P]               ;%planet location vector
V4i =  [V, 0]                   ;%initial orbital velocity of venus (m/s)

initial_values = [P0i, V0i, P1i, V1i, P2i, V2i, P3i, V3i, P4i, V4i];
masses = [m0, m1, m2, m3, m4];

%% universe and time variables
G = 6.67e-11 ; %Nm^2/kg^2
tyears = 1;
tdays = tyears*365;
tseconds = tdays * 24 *60 * 60;
step = tseconds/1000;
time_span = (0 : step : tseconds);

%% ode 45
options = odeset ('Events', @events);
[T, Output, TE, YE, IE] = ode45 (@jupitergoesvroom, time_span, initial_values, options);

%% plot
%   uncomment the following code to display plots

 hold on
 animate_func (T, Output);

%% return value
res = TE;

%% ODE 45 call Function
    function res = jupitergoesvroom (~, W)
        %W comes from ODE
        
        %unpack positions
        %create vector of positions
        %and initialize it so MATLAB doesn't die
        Position_Vector = zeros (2 * Number_of_Bodies, 1);
        for i = 1:Number_of_Bodies
            Position_Vector (2*i-1) = W(4*i-3); % x position of star n
            Position_Vector (2*i)   = W(4*i-2); % y position of star n
        end
        
        % Generate Acceleration values from positions
        Accelerations = acceleration (Position_Vector);
        
        %Rpack into a column vector
        dW = zeros (4 * Number_of_Bodies, 1); %initialize dW vector
        for m = 1 : Number_of_Bodies
            dW(4*m-3) = W(4*m-1);               % change in position x
            dW(4*m-2) = W(4*m);                 % change in position y
            dW(4*m-1) = Accelerations (2*m-1);  % change in velocity x
            dW(4*m)   = Accelerations (2*m);    % change in velocity y
        end
        
        res = (dW) ; %transpose into column vector
    end

%% Acceleration Function
    function res = acceleration (Pos_Vec)
        %% calculate forces on each body from all other bodies
        %gravitational_force = gravity_force_func (m1, P1, m2, P2);
        %initialize force vector so MATLAB doesn't die
        Force = zeros (Number_of_Bodies * 2, 1);
        
        for i = 1:Number_of_Bodies
            %repeat the following force calculations for each body
            %find position of body being acted on
            P1x = Pos_Vec (2*i-1);
            P1y = Pos_Vec (2*i);
            P1 = [P1x; P1y];
            
            %initialize two vectors so MATLAB doesn't die
            forcex = zeros (Number_of_Bodies - 1, 1);
            forcey = zeros (Number_of_Bodies - 1, 1);
            
            for j = 1:Number_of_Bodies
                %make sure it doesn't calculate the force the body exerts on
                %itself
                if (j == i)
                    continue
                end
                %calculate all the forces on a given body
                %find position of body causing the force
                P2x = Pos_Vec (2*j-1);
                P2y = Pos_Vec (2*j);
                P2 = [P2x; P2y];
                
                gravforce = GravForce (masses(i), P1 , masses(j), P2);
                %create vector of all the forces in each direction
                forcex(j) = gravforce(1);
                forcey(j) = gravforce(2);
            end
            
            %sum all values in each direction, then pack into vector
            Force(2*i-1) =  sum(forcex);
            Force(2*i)   =  sum(forcey);
        end
        
        %% derive acceleration from forces (divide by mass)
        %initialize Acceleration vector so MATLAB doesn't die
        Acceleration = zeros(2 * Number_of_Bodies, 1);
        for k = 1:Number_of_Bodies
            Acceleration(2*k-1) = Force(2*k-1) / masses(k);
            Acceleration(2*k)   = Force(2*k)   / masses(k);
        end
        
        %% pack accelerations back into column vector
        res =   (Acceleration);
    end

    function res = GravForce(m1, P1, m2, P2)
        % Finding unit vector from P1 to P2
        % this gives the force ON P1, use this for P1
        R = (P2 - P1);
        R_mag = norm(R);
        r_hat = R/R_mag;
        % Gravitational Force Equation in direction of r_hat
        F_grav = (G * (m1 * m2) / (R_mag^2)) * r_hat;
        res = F_grav;
    end

%% odeset

function [value, isterminal, direction] = events (~,W)
        x = W(13);
        safezone = -1*P3i(1);
        value = [x-safezone];
        isterminal = [1];
        direction = [1];
    end

%% animation
    function animate_func(T,M)
        %% unpack positions
        % animate the positions of the planets, assuming that the
        % columns of M are x1, y1, x2, y2.
        xPositions = zeros(length(M(:,1)), Number_of_Bodies);
        yPositions = zeros(length(M(:,2)), Number_of_Bodies);
        
        for i= 1:Number_of_Bodies
            xPositions(:,i) = M (:, (4*i-3));
            yPositions(:,i) = M (:, (4*i-2));
        end
        
        %% find minimum and maximums
        minmax = [min(min([xPositions])), max(max([xPositions])), min(min([yPositions])), max(max([yPositions]))];
        %minmax = [-5e11, 5e11, -15e11, 15e11];
        %% do the animations
        for l=1:length(T)
            clf;
            axis(minmax);
            hold on;
            draw_func(xPositions, yPositions, l);
            drawnow;
        end
        
        %% plot each body
        function draw_func(xPos, yPos, i) 
            Colors = planetcolors();
            MarkerSizes = [50, 50, 50, 20, 50] ; %GenerateMarkerSizes ();
            for j= 1:Number_of_Bodies
               plot(xPos(i,j), yPos(i,j), '.', 'MarkerSize', MarkerSizes(j), 'MarkerEdgeColor', Colors (j,:));
            end
        end
    end

    function res = planetcolors()
        sun     = [.8 .8 0]; %yellow
        mercury = [0.9 0.4 0.1]; %brown
        venus   = [.9 0.5 0]; %orange
        earth   = [0 0 1]; %blue
        mars    = [1 0 0]; %red
        jupiter = [1 .5 0]; %orange
        saturn  = [1 1 0]; %yellow
        uranus  = [0.5 0.5 1]; %pale blue
        neptune = [0 0 0.5]; %dark blue
        res     = [sun; mercury; venus; earth; mars; 
                   jupiter; saturn; uranus; neptune];
    end

    function res = randomcolors ()
       randcolors = zeros (Number_of_Bodies, 3);
        for i = 1 : Number_of_Bodies
            randcolors(i, 1:3) = [rand, rand, rand];
        end
        res = randcolors;
    end

    function res = GenerateMarkerSizes ()
        %((mass value/ max mass) *40) + 10
        res = ((masses ./ 2e30) .* 40) + 10;
        %mas mass is the sun's mass for the solar system scenario
    end

%% Random Generation
    function res = randMass
        %Range = 1e27 to 1e31
        res = GenerateRandom (1e27, 5e31);
    end
    function res = randPosX
        %Range = 1e10 to 1e12
        res = GenerateRandom (1e10, 1e12);
    end
    function res = randPosY
        %Range = 1e10 to 1e12
        res = GenerateRandom (1e10, 1e12);
    end
    function res = randVelX
        %Range = 1e4 to 1e5
        g = GenerateRandom (-1,1);
        res = GenerateRandom (1e1, 1e4) * (g/ abs(g));
    end
    function res = randVelY
        %Range = 1e4 to 1e5
        g = GenerateRandom (-1,1);
        res = GenerateRandom (1e4, 1e5) * (g/ abs(g));
    end

    function res = GenerateRandom (x1, x2)
        %generating random value in a range
        %starting value + delta * rand
        res = x1 + (rand * abs(x1-x2));
    end

end

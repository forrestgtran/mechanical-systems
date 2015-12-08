% Now in 3d!
function res = Solar_System_3D()
%% initialize values
%up to 64 bodies
Number_of_Bodies = 9;

%% User - Define Initial Values for each Body

%Matrix of 4 wide and number of bodies down
%initialize values
masses = zeros(Number_of_Bodies, 1);             % 1 value for each star
initial_values = zeros(6 * Number_of_Bodies, 1); % 6 values for each body

%sun
m0  =  2e30                    ;%mass of Sun (kg)
P0i =  [0, 0, -1]                  ;%star location vector
V0i =  [0, 0, -10000]                  ;%initial orbital velocity of star 2 (m/s)

% mercury
m1  =  3.3e23                    ;%mass of mercury (kg)
P1i =  [-4.6e10, 0, 0]                  ;%star location vector
V1i =  [-100, 58980, 0]                  ;%initial orbital velocity of mercury (m/s)

% venus
m2  =  4.87e24               ;%mass of body venus (kg)
P2i =  [-1.075e11, 0, 1]       ;%planet location vector
V2i =  [0, 35260, 0]              ;%initial orbital velocity of venus (m/s)

% earth
m3  =  5.976e24              ;%mass of body 3 (kg)
P3i =  [-1.471e11, 0, 0]       ;%planet location vector
V3i =  [0, 30280, 0]              ;%initial orbital velocity of star 3 (m/s)

% mars
m4  =  6.42e23               ;%mass of body 3 (kg)
P4i =  [-2.067e11, 0, 0]       ;%planet location vector
V4i =  [0, 26500, 0]              ;%initial orbital velocity of star 3 (m/s)

% jupiter
m5  =  1.8986e27               ;%mass of body 3 (kg)
P5i =  [-7.409e11, 0, 0]       ;%planet location vector
V5i =  [0, 13720, 0]              ;%initial orbital velocity of star 3 (m/s)

% saturn
m6 =  568.36e24                      ;%mass (kg)
P6i = [-1.352e12, 0, 0]                ;%planet location vector (m)
V6i = [0 , 10180, 0]                ;%initial orbital velocity (m/s)

% Uranus
m7 =  86.816e24                      ;%mass (kg)
P7i = [-2.7413e12, 0, 0]                ;%planet location vector (m)
V7i = [0 , 7110, 0]                ;%initial orbital velocity (m/s)

% Neptune
m8 =  102.42e24                      ;%mass (kg)
P8i = [ -4.44445e12, 0, 0]                ;%planet location vector (m)
V8i = [0 , 5500, 0]                ;%initial orbital velocity (m/s)


initial_values = [P0i, V0i, P1i, V1i, P2i, V2i, P3i, V3i, P4i, V4i P5i, V5i, P6i, V6i, P7i, V7i, P8i, V8i];
masses = [m0, m1, m2, m3, m4 m5 m6, m7, m8];

%% Random Generate all Initial Values for bodies

% for n = 1:Number_of_Bodies
%     initial_values(4*n-3) = randPosX;
%     initial_values(4*n-2) = randPosY;
%     initial_values(4*n-1) = randVelX;
%     initial_values(4*n) = randVelY;
%     masses(n) = randMass;
% end

%% universe and time variables
G = 6.67e-11 ; %Nm^2/kg^2
tyears = 1;
tdays = tyears*365;
tseconds = tdays * 24 *60 * 60;
step = tseconds/1000;
time_span = (0 : step : tseconds);

%% ode 45
%options = odeset ('Events', @events);
[T, Output] = ode45 (@jupitergoesvroom, time_span, initial_values );

%% plot
hold on
%  plot (Output(:,1), Output(:,2), 'y', 'linewidth', 4)
%  plot (Output(:,5), Output(:,6), 'r', 'linewidth', 4)
%  plot (Output(:,9), Output(:,10), 'g', 'linewidth', 4)
%  plot (Output(:,13), Output(:,14), 'r', 'linewidth', 4)
%  plot (Output(:,17), Output(:,18), 'g', 'linewidth', 4)
res = Output;
animate_func (T, Output)
'simulation finished'

%% ODE 45 call Function
    function res = jupitergoesvroom (~, W)
        %W comes from ODE
        
        %unpack positions
        %create vector of positions
        %and initialize it so MATLAB doesn't die
        Position_Vector = zeros (3 * Number_of_Bodies, 1);
        for i = 1:Number_of_Bodies
            Position_Vector (3*i-2) = W(6*i-5); % x position of star n
            Position_Vector (3*i-1) = W(6*i-4); % y position of star n
            Position_Vector (3*i)   = W(6*i-3); % yzposition of star n
        end
        
        % Generate Acceleration values from positions
        Accelerations = acceleration (Position_Vector);
        
        %Rpack into a column vector
        dW = zeros (6 * Number_of_Bodies, 1); %initialize dW vector
        for m = 1 : Number_of_Bodies
            dW(6*m-5) = W(6*m-2);               % change in position x
            dW(6*m-4) = W(6*m-1);               % change in position y
            dW(6*m-3) = W(6*m);                 % change in position y
            
            dW(6*m-2) = Accelerations (3*m-2);  % change in velocity x
            dW(6*m-1) = Accelerations (3*m-1);  % change in velocity y
            dW(6*m)   = Accelerations (3*m);    % change in velocity z
        end
        
        res = (dW) ; %transpose into column vector
    end

%% Acceleration Function
    function res = acceleration (Pos_Vec)
        %% calculate forces on each body from all other bodies
        %gravitational_force = gravity_force_func (m1, P1, m2, P2);
        %initialize force vector so MATLAB doesn't die
        Force = zeros (Number_of_Bodies * 3, 1);
        
        for i = 1:Number_of_Bodies
            %repeat the following force calculations for each body
            %find position of body being acted on
            P1x = Pos_Vec (3*i-2);
            P1y = Pos_Vec (3*i-1);
            P1z = Pos_Vec (3*i);
            P1 = [P1x; P1y; P1z];
            
            %initialize two vectors so MATLAB doesn't die
            forcex = zeros (Number_of_Bodies - 1, 1);
            forcey = zeros (Number_of_Bodies - 1, 1);
            forcez = zeros (Number_of_Bodies - 1, 1);
            
            for j = 1:Number_of_Bodies
                %make sure it doesn't calculate the force the body exerts on
                %itself
                if (j == i)
                    continue
                end
                %calculate all the forces on a given body
                %find position of body causing the force
                P2x = Pos_Vec (3*j-2);
                P2y = Pos_Vec (3*j-1);
                P2z = Pos_Vec (3*j);
                P2 = [P2x; P2y; P2z];
                
                gravforce = GravForce (masses(i), P1 , masses(j), P2);
                %create vector of all the forces in each direction
                forcex(j) = gravforce(1);
                forcey(j) = gravforce(2);
                forcez(j) = gravforce(3);
            end
            
            %sum all values in each direction, then pack into vector
            Force(3*i-2)    =  sum(forcex);
            Force(3*i-1)    =  sum(forcey);
            Force(3*i)      =  sum(forcez);
        end
        
        %% derive acceleration from forces (divide by mass)
        %initialize Acceleration vector so MATLAB doesn't die
        Acceleration = zeros(3 * Number_of_Bodies, 1);
        for k = 1:Number_of_Bodies
            Acceleration(3*k-2) = Force(3*k-2) / masses(k);
            Acceleration(3*k-1) = Force(3*k-1) / masses(k);
            Acceleration(3*k)   = Force(3*k)   / masses(k);
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

%% animation
    function animate_func(T,M)
        %% unpack positions
        % animate the positions of the planets, assuming that the
        % columns of M are x1, y1, x2, y2.
        xPositions = zeros(length(M(:,1)), Number_of_Bodies);
        yPositions = zeros(length(M(:,2)), Number_of_Bodies);
        zPositions = zeros(length(M(:,3)), Number_of_Bodies);
        
        for i= 1:Number_of_Bodies
            xPositions(:,i) = M (:, (6*i-5));
            yPositions(:,i) = M (:, (6*i-4));
            zPositions(:,i) = M (:, (6*i-3));
        end
        
        %% find minimum and maximums
        %minmax = [min(min([xPositions])), max(max([xPositions])), min(min([yPositions])), max(max([yPositions]))];
        minmax = [-3e11, 3e11, -3e11, 3e11, -3e11, 3e11];
        %% do the animations
        for l=1:length(T)
            clf;
            hold on;
            draw_func(xPositions, yPositions, zPositions, l);
            drawnow;
        end
        
        %% plot each body
        function draw_func(xPos, yPos, zPos, i) 
            Colors = planetcolors();
            MarkerSizes = GenerateMarkerSizes ();
            for j= 1:Number_of_Bodies
               scatter3(xPos(i,j), yPos(i,j), zPos(i,j), MarkerSizes(j), Colors (j), 'filled'); %, 'ro', 'MarkerSize', MarkerSizes , 'MarkerEdgeColor', Colors (j,:));
            view (-30,10)
            axis (minmax)
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

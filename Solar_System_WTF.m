% Now with animation!
function res = Solar_System_WTF()
%% initialize values
Number_of_Bodies = 5;
%MatrixWidth = 5;
% for column=1:MatrixWidth
%     for row=1:Number_of_Bodies
%     initial_values(row, column) =
%     end
% end

%% Define Initial Values for each Body

%Matrix of 4 wide and number of bodies down
masses = zeros(Number_of_Bodies, 1);             % 1 value for each star
initial_values = zeros(4 * Number_of_Bodies, 1); % 4 values for each body

% initial_values = [ 0, 0, 0, -131, 778500000000, 0, 0, 13100, -788500000000, 0, 0, -13100,-788500000000, -788500000000, 0, -13100];
% masses = [2e30, 1.8986e27 , 1.8986e27, 1.8986e27];

%sun
m0  =  2e30                    ;%mass of Sun (kg)
P0i =  [0, 0]                  ;%star location vector
V0i =  [0, 0]                  ;%initial orbital velocity of star 2 (m/s)

% mercury
m1  =  3.3e23                    ;%mass of Sun (kg)
P1i =  [-4.6e7, 0]                  ;%star location vector
V1i =  [0, -474000]                  ;%initial orbital velocity of star 2 (m/s)

% venus
m2  =  4.87e24               ;%mass of body 2 (kg)
P2i =  [1.075e8, 0]       ;%planet location vector
V2i =  [0, -350000]              ;%initial orbital velocity of star 2 (m/s)

% earth
m3  =  5.976e24              ;%mass of body 3 (kg)
P3i =  [-1.471e8, 0]       ;%planet location vector
V3i =  [0, -298000]              ;%initial orbital velocity of star 3 (m/s)

% mars
m4  =  6.42e23               ;%mass of body 3 (kg)
P4i =  [2.067e8, 0]       ;%planet location vector
V4i =  [0, -241000]              ;%initial orbital velocity of star 3 (m/s)

% jupiter
m5  =  1.8986e27               ;%mass of body 3 (kg)
P5i =  [-7.409e8, 0]       ;%planet location vector
V5i =  [0, -13100]              ;%initial orbital velocity of star 3 (m/s)



initial_values = [P0i, V0i, P1i, V1i, P2i, V2i, P3i, V3i, P4i, V4i];%, P5i, V5i];% P6i, V6i, P7i, V7i, P8i, V8i];
masses = [m0, m1, m2, m3, m4];%, m5];%, m6, m7, m8];

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
tdays = 365;
tseconds = tdays * 24 *60 * 60;
step = tseconds/10;
time_span = [0, tseconds]; %(0 : step : tseconds);

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
res = Output
animate_func (T, Output)


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
        %minmax = [min(min([xPositions])), max(max([xPositions])), min(min([yPositions])), max(max([yPositions]))];
        minmax = [-3e8, 3e8, -3e8, 3e8];
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
            for j= 1:Number_of_Bodies
               plot(xPos(i,j), yPos(i,j), 'r.', 'MarkerSize', 20);
            end
        end
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

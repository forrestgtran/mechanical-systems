%this code calls 'rocket1' in a for loop with different initial velocities
%for the space ship

%x velocities
xvelocities = [0, 1, 2, 5, 10, 20, 50, 1e2, 2e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6];
yvelocities = [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0, ...
               -1, -2, -5, -10, -20, -50, -100, -200, -500];

Output = zeros (length(xvelocities), length(yvelocities));            

parfor i = 1 : length(xvelocities)
    disp(i)
    parfor j = 1 : length(yvelocities)
        try
        Output(i,j) = rocket1 (xvelocities(i), yvelocities(j));
        catch
            Output(i,j) = 0;
            continue
        end
    end
end

pcolor(xvelocities, yvelocities, Output)
shading interp
xlabel ('X Velocity (m/s)')
ylabel ('Y Velocity (m/s)')
colorbar
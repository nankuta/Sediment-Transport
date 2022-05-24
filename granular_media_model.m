n_particles = 1000;

%set properties of particles and environment
g = .02;
pho_particle = 10;
pho_air = 1;
m_particle = .1;
nu = 1.5*10^-5;
sigma_p = pho_particle/pho_air;
d = ((6*m_particle)/(pi*pho_particle))^(1/3);
reynolds_particle = 1*d/nu;

%position
x = randn(n_particles, 1);
y = randn(n_particles, 1);
z = abs(randn(n_particles, 1));

%velocity
u = randn(n_particles, 1);
v = randn(n_particles, 1);
w = abs(randn(n_particles, 1));

%acceleration
du = zeros(n_particles, 1);
dv = zeros(n_particles, 1);
dw = zeros(n_particles, 1);

%value indicating if given particle is in motion
trajectories = ones(n_particles, 1);

%time parameters
t_final = 20000;
n_steps = 50000;
dt = t_final/n_steps;

%field of view
x_bounds = [-2, 2];
y_bounds = [-2, 2];
z_bounds = [-2, 2];

%create random terrain using gradient noise
coefficients = randn(2+ x_bounds(2)-x_bounds(1), 2+ y_bounds(2)-y_bounds(1), 2);
for a = 1:length(coefficients)
    for b = 1:length(coefficients)
        mag = norm([coefficients(a,b, 1), coefficients(a,b, 2)]);
        coefficients(a,b, 1) = coefficients(a,b, 1)/mag;
        coefficients(a,b, 2) = coefficients(a,b, 2)/mag;
    end
end

%number of points rendered on given axis of terrain
npoints = 100;

%generate heightmap
x_terrain = linspace(1, 1+x_bounds(2)-x_bounds(1), npoints);
y_terrain = linspace(1, 1+y_bounds(2)-x_bounds(1), npoints);
heightmap = gradheightmap(x_terrain, y_terrain, coefficients);

%set field to reflect area and boundary parameters
x_terrain = linspace(x_bounds(1), 1+x_bounds(2), npoints);
y_terrain = linspace(y_bounds(1), 1+y_bounds(2), npoints);

%magnitude of intial velocity for stochastically emitted particles
liftoffmag = 25;

for t = 1:n_steps
    clf
    for i = 1:n_particles
        
        %stochastically put stationary particles in motion
        trajectories = trajectories + floor(abs(randn(n_particles, 1)));
        
        if (trajectories(i) == 0)
            continue
        end
        
        %intialize velocity if particle was previously stationary
        if (u(i) == 0 && v(i) == 0 && w(i) == 0)
            u(i) = liftoffmag*randn;
            v(i) = liftoffmag*randn;
            w(i) = liftoffmag*abs(randn);
        end
        
        x_prev = x(i);
        y_prev = y(i);
        z_prev = z(i);
        
        %update accleration
        [du(i), dv(i), dw(i)] = d_vel(x(i), y(i), z(i), u(i), v(i), w(i), g, sigma_p, d, nu);
        
        %update velocity
        u(i) = dt*du(i);
        v(i) = dt*dv(i);
        w(i) = dt*dw(i);
        
        %update position
        x(i) = x(i) + dt*u(i);
        y(i) = y(i) + dt*v(i);
        z(i) = z(i) + dt*w(i);
      
        %particles bounce upon hitting boundary
        if(x(i) > x_bounds(2))
            x(i) = x_bounds(2)-.01;
            u(i) = -u(i);
        end
        if(y(i) > y_bounds(2))
            y(i) = y_bounds(2)-.01;
            v(i) = -v(i);
        end
        if(x(i) < x_bounds(1))
            x(i) = x_bounds(1)+.01;
            u(i) = -u(i);
        end
        if(y(i) < y_bounds(1))
            y(i) = y_bounds(1)+.01;
            v(i) = -v(i);
        end
        
        %find closest value x,y in heightmap
        [val,idxx]=min(abs(x_terrain-x(i)));
        [val, idxy]=min(abs(y_terrain-y(i)));
        
        %deposition
        %particle becomes stationary upon colliding with terrain
        if(z(i) <= heightmap(idxy, idxx))
            u(i) = 0;
            v(i) = 0;
            w(i) = 0;
            x(i) = x_prev;
            y(i) = y_prev;
            z(i) = z_prev;
            %reset movement status
            trajectories(i) = -100;
        end
    end
    
    scatter3(x, y, z, "red"); hold on
    surf(x_terrain, y_terrain, heightmap); hold on
    
    xlim(x_bounds);
    ylim(y_bounds);
    zlim(z_bounds);
    
    if (mod(t, 10) == 0)
        drawnow
    end
end

%calculate acceleration in a shear flow
function [du, dv, dw] = d_vel(x, y, z, u, v, w, g, sigma_p, d, nu)

    %set various parameters
    shear_mag = .2;
    r = norm([x, y]);
    reynolds_particle = shear_mag*d/(nu*r^2);
    C_d = (24/reynolds_particle)*(1+.15*(reynolds_particle)^.687);

    c = -.75*C_d/(sigma_p*d);
    du = c*shear_mag*((u-y)/(r^2));
    dv = c*shear_mag*((v+x)/(r^2));
    dw = w*shear_mag*(x/(r^2)) - g;
end

%interpolation function used in generating height map
function val = interp(a0, a1, w)
    %val = (a1 - a0)*(3 - w*2)*w*w + a0;
    val = (a1 - a0) * ((w * (w * 6.0 - 15.0) + 10.0) * w * w * w) + a0;
end

%create a height map using a gradient noise algorithm
function height = gradheightmap(X, Y, coefficients)

    height = zeros(length(X), length(Y));
    for r = 1:length(X)
        for c = 1:length(Y)
            x = X(r);
            y = Y(c);
            x_floor = floor(x);
            y_floor = floor(y);
            
            vi = [x_floor-x, y_floor-y];
            gradi = [coefficients(x_floor,y_floor, 1), coefficients(x_floor,y_floor, 2)];
            n0 = dot(gradi, vi/sqrt(2));
            vi = [1+x_floor-x, y_floor-y];
            gradi = [coefficients(x_floor+1, y_floor, 1), coefficients(x_floor+1, y_floor, 2)];
            n1 = dot(gradi, vi/sqrt(2));
            ix0 = interp(n0, n1, x-x_floor);

            vi = [x_floor-x, 1+y_floor-y];
            gradi = [coefficients(x_floor, y_floor+1, 1), coefficients(x_floor, y_floor+1, 2)];
            n0 = dot(gradi, vi/sqrt(2));
            vi = [1+x_floor-x, 1+y_floor-y];
            gradi = [coefficients(x_floor+1, y_floor+1, 1), coefficients(x_floor+1, y_floor+1, 2)];
            n1 = dot(gradi, vi/sqrt(2));
            ix1 = interp(n0, n1, x-x_floor);
            height(r,c) = 2*interp(ix0, ix1, y-y_floor);
        end
    end
end

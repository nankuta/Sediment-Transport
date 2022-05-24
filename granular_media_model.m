n_particles = 5000;
g = .1;

%set properties of particles
pho_particle = 10;
pho_air = 1;
m_particle = .1;
nu = 1.5*10^-5;

sigma_p = pho_particle/pho_air;
d = ((6*m_particle)/(pi*pho_particle))^(1/3);
reynolds_particle = 1*d/nu;

%position
x = abs(randn(n_particles, 1));
y = abs(randn(n_particles, 1));
z = 1+abs(randn(n_particles, 1));

%velocity
u = randn(n_particles, 1);
v = randn(n_particles, 1);
w = randn(n_particles, 1);

%acceleration
du = zeros(n_particles, 1);
dv = zeros(n_particles, 1);
dw = zeros(n_particles, 1);

t_final = 20000;
n_steps = 10000;
dt = t_final/n_steps;

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

npoints = 50;
linscale = npoints/(x_bounds(2)-x_bounds(1));
linshift = -x_bounds(1);

x_terrain = linspace(1, 1+x_bounds(2)-x_bounds(1), npoints);
y_terrain = linspace(1, 1+y_bounds(2)-x_bounds(1), npoints);

heightmap = gradheightmap(x_terrain, y_terrain, coefficients);
x_terrain = linspace(x_bounds(1), 1+x_bounds(2), npoints);
y_terrain = linspace(y_bounds(1), 1+y_bounds(2), npoints);

for t = 1:n_steps
    clf
    for i = 1:n_particles
        [du(i), dv(i), dw(i)] = d_vel(x(i), y(i), z(i), u(i), v(i), w(i), g, sigma_p, d, nu, t);
        
        u(i) = dt*du(i);
        v(i) = dt*dv(i);
        w(i) = dt*dw(i);
        
        x(i) = x(i) + dt*u(i);
        y(i) = y(i) + dt*v(i);
        z(i) = z(i) + dt*w(i);
      
        %particles bounce upon hitting boundary
        if(x(i) > x_bounds(2))
            x(i) = x_bounds(2)-.01;
            u(i) = -.5*u(i);
        end
        if(y(i) > y_bounds(2))
            y(i) = y_bounds(2)-.01;
            v(i) = -.5*v(i);
        end
        if(z(i) > z_bounds(2))
            w(i) = z_bounds(1)+(z(i)-z_bounds(2));
        end
        
        if(x(i) < x_bounds(1))
            x(i) = x_bounds(1)+.01;
            u(i) = -.5*u(i);
        end
        if(y(i) < y_bounds(1))
            y(i) = y_bounds(1)+.01;
            v(i) = -.5*v(i);
        end
        
        [val,idxx]=min(abs(x_terrain-x(i)));
        [val, idxy]=min(abs(y_terrain-y(i)));
        
        %keep particles from falling through terrain
        if(z(i) <= heightmap(idxy, idxx))
%             dU = m_particle*g*(heightmap(idxy,idxx)-z(i));
%             dV = (2*dU/m_particle)/(u(i)^2+v(i)^2+w(i)^2);
%             u(i) = u(i)*(1-dV);
%             v(i) = v(i)*(1-dV);
            w(i) = -w(i);
            z(i) = heightmap(idxy, idxx);
        end
    end
    
    scatter3(x, y, z); hold on
    surf(x_terrain, y_terrain, heightmap); hold on
    
    
    xlim(x_bounds);
    ylim(y_bounds);
    zlim(z_bounds);
    
    if (mod(t, 10) == 0)
        drawnow
    end
end

%calculate accleration in a rotational flow
function [du, dv, dw] = d_vel(x, y, z, u, v, w, g, sigma_p, d, nu, t)

    %set various parameters
    shear_mag = .05;
    r = norm([x, y]);
    reynolds_particle = shear_mag*d/(nu*r^2);
    C_d = (24/reynolds_particle)*(1+.15*(reynolds_particle)^.687);

    c = -.75*C_d/(sigma_p*d);
    du = c*shear_mag*((u-y)/(r^2));
    dv = c*shear_mag*((v+x)/(r^2));
    dw = (1+sin(t/100))*w*shear_mag*(x/(r^2)) - g;
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
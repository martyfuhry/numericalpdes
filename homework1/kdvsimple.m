function kdvsimple
    dx=0.1;
    xmin = -8;
    xmax = +8;
    ymin = -20;
    ymax = +20;
    x=(xmin+dx:dx:xmax)';

    k = dx^3;
    nsteps = 100/k;
    mesh2d = [];

    u = 20*exp(-x.^2);

    for ii = 1:nsteps
        k1 = k * kdvequ(u,dx);
        k2 = k * kdvequ(u+k1/2,dx);
        k3 = k * kdvequ(u+k2/2,dx);
        k4 = k * kdvequ(u+k3,dx);
        u = u + k1/6 + k2/3 + k3/3 + k4/6;
        mesh2d(ii,:) = -u;

        if mod(ii,10) == 0
            plot(x,u,'b-','LineWidth',2);
            axis([xmin,xmax,ymin,ymax])
            drawnow;
        end
    end

    function dudt=kdvequ(u,dx)
        u = [u(end-1:end); u; u(1:2)];
        dudt = -(u(4:end-1)-u(2:end-3))/(2*dx) - (u(5:end)-2*u(4:end-1)+2*u(2:end-3)-u(1:end-4))/(2*dx^3);

function cell = MakeNewCell(mothercell, isdaughter, sister, NN, ks)

    x = mothercell.position;
    u = mothercell.orientation;
    L = mothercell.length;
    r = mothercell.radius; 
    c1 = mothercell.connectedto;
    connect = [];

    %gavg = 0.224; %0.0862; %exponential growth rate
    %gavg = 0.2199; %linear growth rate in micrometers/hr
    g = 0.224; %gavg; % + 0.01*randn; %randomly assigned growth rate from N(mu = a, sigma = b)
    %rotation is clockwise: positive theta rotates down
    R = @(th) [cosd(th) -sind(th); sind(th) cosd(th)]; %rotation matrix
    L = (L - (pi*r/2))/2;

    if(isdaughter == 1) %this is the 'right' new daughter cell
        xnew = x + (r + L/2)*u;
        isLRd = 2;
        if (NN == 3 || NN == 4) %force M shape
            unew = u*R(-15);
        elseif (NN == 5 || NN == 8)
            unew = u*R(10);
        else
            unew = u; %*R(randn);
        end
        for i = 1:size(c1, 2)
            if (c1(2, i) == 2) %if mother cell is connected to another cell at right, keep this connection
                connect = [connect, c1(:, i)];
            end
        end
        connect = [connect, [sister 1 0 ks 1]'];
    else %this is the 'left' new daughter cell (will replace its mother)
        xnew = x - (r + L/2)*u;
        isLRd = 1;
        if (NN == 3 || NN == 4)
            unew = u*R(15);
        elseif (NN == 5 || NN == 8)
            unew = u*R(-10);
        else
            unew = u; %*R(randn);
        end
        for i = 1:size(c1,2)
            if (c1(2, i) == 1) %if mother cell is connected to another cell at left, keep this connection
                connect = [connect, c1(:, i)];
            end
        end
        connect = [connect, [sister 2 0 ks 1]'];
    end
    
    if unew(1) < 0
    unew = unew*[-1 0; 0 -1];
    end
    unew = unew./norm(unew);
     
    cell = struct('position', xnew, 'orientation', unew, 'length', L, 'radius', r, 'p', NaN, 'q', NaN, 'mass', NaN, 'moment', NaN, 'growrate', g, 'v0', mothercell.v0, 'w0', mothercell.w0, 'contactpts', [0; 0; 0], 'connectedto', connect, 'isLRd', isLRd, 'oldL', mothercell.oldL);

end

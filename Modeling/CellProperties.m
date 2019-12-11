function cell = CellProperties(cell1)
%this function takes a cell's center, orientation, and length and fills out
%the rest of the fields such as p, q, mass, moment

r = cell1.radius; %radius of cell
rho = 1; %cell density
x = cell1.position;
u = cell1.orientation;
L = cell1.length; %length
L0 = cell1.oldL;
isLRd = cell1.isLRd;

if L > 1.1
     isLRd = 0;
end

u = u./norm(u); %normalize orientation vector
if isLRd == 1 %left daughter cell
    x = x + ((L - L0)/2)*u;
    p = x - (L/2)*u;
    q = x + (L/2)*u;
elseif isLRd == 2 %right daughter cell
    x = x - ((L - L0)/2)*u;
    p = x - (L/2)*u;
    q = x + (L/2)*u;
elseif isLRd == 0 %grow evenly
    p = x - (L/2)*u; %left endpoint of body axis
    q = x + (L/2)*u; %right endpoint of body axis
end

mass = rho*pi*[(4/3)*r^3 + L*r^2]; %cell mass
moment = rho*pi*[(1/12)*r^2*L^3 + (3/4)*L*r^4 + (16/30)*r^5 + (1/12)*r^3*L^2]; %moment of inertia
g = cell1.growrate;

cell = struct('position', x, 'orientation', u, 'length', L, 'radius', r, 'p', p, 'q', q, 'mass', mass, 'moment', moment, 'growrate', g, 'v0', cell1.v0, 'w0', cell1.w0, 'contactpts', cell1.contactpts, 'connectedto', cell1.connectedto, 'isLRd', isLRd, 'oldL', cell1.oldL);

end
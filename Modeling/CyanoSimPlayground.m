%clear all

function CellArray = CyanoSimPlayground(k, bt, br, ks)
%force parameters: k cell-cell interaction, bt translational friction, br
%rotational friction, ks spring constant (force of cell-cell connections)
%% 

R = @(th) [cosd(th) -sind(th); sind(th) cosd(th)];
Ccount = 0;

%% give initial condition
%define the first cell horizontally aligned, still and centered at (0, 0)
CellArray.position = [0 0]; %[-1.5 0];
CellArray.orientation = [1 0];
CellArray.length = 1;
CellArray.radius = 1;
CellArray.growrate = 0.224;
CellArray.v0 = [0 0];
CellArray.w0 = [0 0 0];
CellArray.contactpts = zeros(3, 1);
CellArray.connectedto = []; %[2; 0; 0; ks; 1]; %[(index of linked cell) (0 = connected from center, 1 = connected at left, 2 = connected at right) (time connected) (spring constant) (1 = can be broken)]';
CellArray.isLRd = 0; %1; %for recently divided cells, 1 = left daughter, 2 = right daughter, 0 = not recently divided (L >= (1/2)*Ldiv)
CellArray.oldL = 1;
CellArray = CellProperties(CellArray(1));

% CellArray(2).position = [2.2 0];
% CellArray(2).orientation = [0 -1]*R(15);
% CellArray(2).length = 1;
% CellArray(2).radius = 1;
% CellArray(2).growrate = 0; %0.224;
% CellArray(2).v0 = [0 0];
% CellArray(2).w0 = [0 0 0];
% CellArray(2).contactpts = zeros(3, 1);
% CellArray(2).connectedto = [1; 0; 0; ks; 1]; %[1; 1; 0; ks; 1]; %[(index of linked cell) (1 = connected at left, 2 = connected at right) (time connected) (spring constant)]';
% CellArray(2).isLRd = 2;
% CellArray(2).oldL = 1;
% CellArray(2) = CellProperties(CellArray(2));

%% declare constants
r = 1; %cell radius
l0 = 0.01; %0.2; %equilibrium length of spring
Ldiv = 2.5; %division length

tc = 0.01;
dt = 0.25;
tfinal = 25; %32;

nframe = tfinal/dt;
mov(1:nframe) = struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren')
cframe = 1;
A = zeros(tfinal*(1/dt), 17);
A(:, 1) = 1:(tfinal*(1/dt))';

%%
while(tc < tfinal)
    %update cell growth and division
    
% if cframe == 107
%     CellArray(10).orientation = CellArray(10).orientation*R(10);
%     CellArray(2).position(2) = CellArray(2).position(2) + 0.1; 
%     CellArray(13).position(2) = CellArray(13).position(2) + 0.1; 
%     CellArray(12).position(2) = CellArray(12).position(2) - 0.1; 
%     CellArray(11).position(2) = CellArray(11).position(2) - 0.1; 
% end
% if cframe == 80
%     CellArray(8).orientation = CellArray(8).orientation*R(-10);
%     CellArray(2).orientation = CellArray(2).orientation*R(10);
% end

    N = numel(CellArray); %only grow cells present at the beginning of this time step
    for n = 1:N
        f = @(t, y) CellArray(n).growrate*y;
        [tout, yout] = ode45(f, [tc (tc + dt)], CellArray(n).length);
        L = yout(end);
        CellArray(n).oldL = CellArray(n).length;
        CellArray(n).length = L;
        CellArray(n) = CellProperties(CellArray(n));
        
%% update cell-cell connections        
        if(~isempty(CellArray(n).connectedto))
            c1 = CellArray(n).connectedto;
            Csize = size(c1, 2);
            for ii = 1:Csize
                c1(3, ii) = c1(3, ii) + 1; %update time connected
                if c1(5, ii) == 1
                    c1(4, ii) = ks*(1 - c1(3, ii)/55); %spring constant, decreases linearly with time
                end
                if (CellArray(n).length > Ldiv) %this cell is about to divide
                    ncell = c1(1, ii);
                    side = c1(2, ii);
                    if side == 1
                        c2 = CellArray(ncell).connectedto;
                        col = find(c2(1, :) == n); 
                        c2(1, col) = n;
                        CellArray(ncell).connectedto = c2;
                    elseif side == 2
                        c2 = CellArray(ncell).connectedto;
                        col = find(c2(1, :) == n);
                        c2(1, col) = (numel(CellArray) + 1);
                        CellArray(ncell).connectedto = c2;
                    end
                end
            end
            for ii = fliplr(1:Csize)
                %if(c1(3, ii) == 25) %shortly after division
                if(abs(c1(4, ii)) < 1) %if spring constant is low
                    c1(:, ii) = []; %delete the connection
                end
            end
            CellArray(n).connectedto = c1;
        end
        if(CellArray(n).length > Ldiv)
            
 %% perform cell division
            NN = (numel(CellArray) + 1);
            CellArray(n).oldL = 0.5;
            CellArray(NN) = MakeNewCell(CellArray(n), 1, n, NN, ks); %add new daughter cell
            CellArray(NN) = CellProperties(CellArray(NN));
            CellArray(n) = MakeNewCell(CellArray(n), 0, NN, NN, ks); %replace mother with next daughter cell
            CellArray(n) = CellProperties(CellArray(n));
        end
    end
        
%% mechanics, construct force matrix (cell-cell interactions)
    N = numel(CellArray);
%     k %strength of cell-cell interaction
%     bt %translational friction
%     br %rotational friction
    for kk = 1:N
        %clear contact points
        CellArray(kk).contactpts = zeros(3, 1);
    end
    Fij = zeros(N, N, 2); %Fij = force on cell i exerted by cell j
    C = zeros(N, N, 2); %C_ij = dist from center of cell i to contact point in cell i of force by cell j
    TS = zeros(N, 1); %torque from spring force
    
%% find force on cell i from cell j spring connection    
    for i = 1:N-1
       for j = (i + 1):N
           if(~isempty(CellArray(j).connectedto))
               c2 = CellArray(j).connectedto;
               for jj = 1:size(c2, 2)
                   if(c2(1, jj) == i)
                       if(c2(2, jj) == 1) %connected at left on j
                           xx1 = CellArray(i).q + CellArray(i).orientation; %right end of i
                           xx2 = CellArray(j).p - CellArray(j).orientation; %left end of j
                           n1 = xx2 - xx1; %vector from i to j 
                           if (norm(n1) > l0)
                               theta1 = acos(n1(1)/norm(n1));
                               theta2 = asin(n1(2)/norm(n1));
                               d1 = [n1(1) - l0*cos(theta1); n1(2) - l0*sin(theta2)]; %displacement of spring
                               Fij(i, j, 1) = c2(4, jj)*d1(1);
                               Fij(i, j, 2) = c2(4, jj)*d1(2);
                               Cspr1 = xx1 - CellArray(i).position; %CellArray(i).q - CellArray(i).position;
                               Cspr2 = CellArray(j).position - xx2; %CellArray(j).position - CellArray(j).p; 
                               TS(i) = Cspr1(1)*Fij(i, j, 2) - Cspr1(2)*Fij(i, j, 1);
                               TS(j) = Cspr2(1)*(Fij(i, j, 2)) - Cspr2(2)*(Fij(i, j, 1));
                           end
                       elseif(c2(2, jj) == 2); %connected at right on j
                           xx1 = CellArray(j).q + CellArray(j).orientation; %right end of j
                           xx2 = CellArray(i).p - CellArray(i).orientation; %left end of i                         
                           n1 = xx1 - xx2; %vector from i to j 
                           if (norm(n1) > l0)
                               theta1 = acos(n1(1)/norm(n1));
                               theta2 = asin(n1(2)/norm(n1));
                               d1 = [n1(1) - l0*cos(theta1); n1(2) - l0*sin(theta2)]; %displacement of spring
                               Fij(i, j, 1) = c2(4, jj)*d1(1);
                               Fij(i, j, 2) = c2(4, jj)*d1(2);
                               Cspr1 = xx1 - CellArray(j).position;
                               Cspr2 = CellArray(i).position - xx2;
                               TS(j) = Cspr1(1)*(Fij(i, j, 2)) - Cspr1(2)*(Fij(i, j, 1));
                               TS(i) = Cspr2(1)*(Fij(i, j, 2)) - Cspr2(2)*(Fij(i, j, 1));
                           end
                       elseif(c2(2, jj) == 0); %connected at center
                           if CellArray(j).position(1) > CellArray(i).position(1) %j is right of i
                               ss1 = sign(CellArray(i).orientation(2));
                               ss2 = sign(CellArray(j).orientation(2));
                               if (ss1 == 1 && ss2 == -1 || ss1 == -1 && ss2 == 1)
                                   if CellArray(i).p - CellArray(j).q < CellArray(i).q - CellArray(j).p
                                       xx1 = CellArray(i).q + CellArray(i).orientation*R(ss1*90);
                                       xx2 = CellArray(j).p + CellArray(j).orientation*R(ss1*90);
                                   else
                                       xx1 = CellArray(i).p + CellArray(i).orientation*R(ss1*90);
                                       xx2 = CellArray(j).q + CellArray(j).orientation*R(ss1*90);
                                   end
                               else
                                   if CellArray(i).p - CellArray(j).p < CellArray(i).q - CellArray(j).q
                                       xx1 = CellArray(i).q + CellArray(i).orientation*R(ss1*90);
                                       xx2 = CellArray(j).q + CellArray(j).orientation*R(ss1*90);
                                   else
                                       xx1 = CellArray(i).p + CellArray(i).orientation*R(ss1*90);
                                       xx2 = CellArray(j).p + CellArray(j).orientation*R(ss1*90);
                                   end
                               end
                               %xx1 = CellArray(i).position + CellArray(i).orientation*R(ss1*90); %right side of i
                               %xx2 = CellArray(j).position + CellArray(j).orientation*R(ss2*(-90)); %left side of j
                               n1 = xx2 - xx1; %vector from i to j 
                               if (norm(n1) > l0)
                                   theta1 = acos(n1(1)/norm(n1));
                                   theta2 = asin(n1(2)/norm(n1));
                                   d1 = [n1(1) - l0*cos(theta1); n1(2) - l0*sin(theta2)]; %displacement of spring
                                   Fij(i, j, 1) = c2(4, jj)*d1(1);
                                   Fij(i, j, 2) = c2(4, jj)*d1(2);
                                   Cspr1 = CellArray(i).orientation*R(ss1*90); %CellArray(i).q - CellArray(i).position;%xx1 - CellArray(i).position;
                                   Cspr2 = -CellArray(j).orientation*R(ss2*(-90)); %CellArray(j).position - CellArray(j).p; %CellArray(j).position - xx2;  
                                   TS(i) = Cspr1(1)*Fij(i, j, 2) - Cspr1(2)*Fij(i, j, 1);
                                   TS(j) = Cspr2(1)*(Fij(i, j, 2)) - Cspr2(2)*(Fij(i, j, 1));
                               end
                           else %j is left of i 
                               ss1 = sign(CellArray(i).orientation(2));
                               ss2 = sign(CellArray(j).orientation(2));
                               if (ss1 == 1 && ss2 == -1 || ss1 == -1 && ss2 == 1)
                                   if CellArray(j).p - CellArray(i).q < CellArray(j).q - CellArray(i).p
                                       xx1 = CellArray(j).q + CellArray(j).orientation*R(ss1*90);
                                       xx2 = CellArray(i).p + CellArray(i).orientation*R(ss1*90);
                                   else
                                       xx1 = CellArray(i).p + CellArray(i).orientation*R(ss1*90);
                                       xx2 = CellArray(j).q + CellArray(j).orientation*R(ss1*90);
                                   end
                               else
                                   if CellArray(j).p - CellArray(i).p < CellArray(j).q - CellArray(i).q
                                       xx1 = CellArray(j).q + CellArray(j).orientation*R(ss1*90);
                                       xx2 = CellArray(i).q + CellArray(i).orientation*R(ss1*90);
                                   else
                                       xx1 = CellArray(j).p + CellArray(j).orientation*R(ss1*90);
                                       xx2 = CellArray(i).p + CellArray(i).orientation*R(ss1*90);
                                   end
                               end
%                                xx1 = CellArray(i).position + CellArray(i).orientation*R(ss1*(-90)); %left side of i
%                                xx2 = CellArray(j).position + CellArray(j).orientation*R(ss2*90); %right side of j
                               n1 = xx2 - xx1; %vector from i to j 
                               if (norm(n1) > l0)
                                   theta1 = acos(n1(1)/norm(n1));
                                   theta2 = asin(n1(2)/norm(n1));
                                   d1 = [n1(1) - l0*cos(theta1); n1(2) - l0*sin(theta2)]; %displacement of spring
                                   Fij(i, j, 1) = c2(4, jj)*d1(1);
                                   Fij(i, j, 2) = c2(4, jj)*d1(2);
                                   %Cspr1 = CellArray(i).q - CellArray(i).position; %xx1 - CellArray(i).position; 
                                   %Cspr2 = CellArray(j).p - CellArray(j).position; %CellArray(j).position - xx2; 
                                   Cspr1 = CellArray(i).orientation*R(ss1*(-90)); %xx1 - CellArray(i).position; 
                                   Cspr2 = -CellArray(j).orientation*R(ss2*90); %CellArray(j).position - xx2; 
                                   TS(i) = Cspr1(1)*Fij(i, j, 2) - Cspr1(2)*Fij(i, j, 1);
                                   TS(j) = Cspr2(1)*(Fij(i, j, 2)) - Cspr2(2)*(Fij(i, j, 1));
                               end
                           end
                       end
                   end
               end
           end
           
%% find cell-cell forces, contact points
           [s, n, t1, t2] = DistBetween2Segment(CellArray(i).p, CellArray(i).q, CellArray(j).p, CellArray(j).q);
           n = n./norm(n);
           if(2*r > s)
               %include force on cell i by cell j from cell-cell interaction
               Dij = 2*r - s;
               Fij(i, j, 1) = Fij(i, j, 1) + k*Dij*n(1); %x component
               Fij(i, j, 2) = Fij(i, j, 2) + k*Dij*n(2); %y component
               dij = t1 - CellArray(i).position;
               C(i, j, 1) = dij(1);
               C(i, j, 2) = dij(2);
               d2ij = t2 - CellArray(j).position;
               C(j, i, 1) = d2ij(1);
               C(j, i, 2) = d2ij(2);
               M = [(t1(1) + t2(1))/2; (t1(2) + t2(2))/2; sqrt(Fij(i, j, 1)^2 + Fij(i, j, 2)^2)];
               CellArray(i).contactpts = [CellArray(i).contactpts, M];
               CellArray(j).contactpts = [CellArray(j).contactpts, M];
           end
       end
    end   
    Fij(:, :, 1) = Fij(:, :, 1) - Fij(:, :, 1)';
    Fij(:, :, 2) = Fij(:, :, 2) - Fij(:, :, 2)';
    Fi = zeros(N, 2);
    for i = 1:N
        Fi(i, :) = [sum(Fij(i, :, 1)) sum(Fij(i, :, 2))]; %net cell-cell force on cell i as [x y]
    end
    
%% integrate ODEs for position, angular velocity, and orientation
    for i = 1:N
        tspan = [tc (tc + dt)];
        y0 = [CellArray(i).position, CellArray(i).v0, CellArray(i).w0, CellArray(i).orientation, 0];
        m = CellArray(i).mass;
        mo = CellArray(i).moment;
        tau = 0;
        for j = 1:N
            a = C(i, j, 1)*Fij(i, j, 2) - C(i, j, 2)*Fij(i, j, 1);
            tau = tau + a;
        end
        torq = [0 0 (tau + TS(i))];
        [t, y] = ode45(@(t, y) newtonODEs(t, y, m, mo, bt, br, Fi(i, :)', torq'), tspan, y0);
        CellArray(i).position = y(end, 1:2);
        CellArray(i).v0 = y(end, 3:4);
        CellArray(i).w0 = y(end, 5:7);
        CellArray(i).orientation = y(end, 8:9)./norm(y(end, 8:9));
    end

   %% update cell growth rate by dist from edge
%     v = [CellArray.position];
%     v = v(1:2:end);
%     minx = min(v);
%     maxx = max(v);
%     for i = 1:numel(CellArray)
%         p = CellArray(i).position(1);
%         if p > 0
%             xdist = maxx - p;
%         else
%             xdist = minx - p;
%         end
%         xdist = floor(abs(xdist));
%         if xdist == 1
%             CellArray(i).growrate = 0.223;
%         elseif xdist == 2
%             CellArray(i).growrate = 0.222;
%         elseif xdist == 3
%             CellArray(i).growrate = 0.221;
%         elseif xdist > 3
%             CellArray(i).growrate = 0.22;
%         end
%     end

% Make growth rates proportional to force acting on a cell
for i = 1:numel(CellArray)
    s1 = CellArray(i).contactpts;
    ls = size(s1, 2);
    fc = 0;
    if (ls > 1)
        for m = 2:ls
            if (s1(3, m) >= 15)
                fc = fc + 5;
            elseif((s1(3, m) < 15) && (s1(3, m) >= 10))
                fc = fc + 4;
            elseif((s1(3, m) < 10) && (s1(3, m) >= 5))
                fc = fc + 3;
            else
                fc = fc + 2;
            end
        end
    end
    if fc >= 10
        CellArray(i).growrate = 0.22; %yellow
    elseif (fc < 10 && fc >= 7)
        CellArray(i).growrate = 0.221; %orange
    elseif (fc < 7 && fc >= 5)
        CellArray(i).growrate = 0.222; %pink
    elseif (fc < 5 && fc >= 2)
        CellArray(i).growrate = 0.223; %purple
    else
        CellArray(i).growrate = 0.224; %blue
    end
    A(cframe, i) = fc;
end
                
    

%% connect cells to force alignment
%     if (N >= 4 && conncount == 1)
%         for j = 1:N-1
%             %[ind, dist] = getClosestCell(CellArray, j);
%                 for i = (j+1):N
%                     otr = 0;
%                     disp(i), disp(j)
%                     cnto = CellArray(j).connectedto
%                     for kbc = 1:size(cnto, 2)
%                         ab = cnto(1, kbc);
%                         if ab == i
%                             otr = 1;
%                             disp('otr'), disp(otr)
%                         end
%                     end
%                     if (abs(CellArray(j).position(1)) < 2.2 && CellArray(i).position(1) > CellArray(j).position(1) && otr == 0)
%                         if (CellArray(i).position(1) < (CellArray(j).position(1) + 2.2))
%                             if ((CellArray(i).position(2) < CellArray(j).position(2) + .5) || (CellArray(i).position(2) > CellArray(j).position(2) - .5))
%                                 CellArray(i).orientation = CellArray(i).orientation*R(2);
%                                 CellArray(j).orientation = CellArray(j).orientation*R(-2);
%                                 CellArray(i).connectedto = [j; 0; 5; 2*ks; 0];
%                                 CellArray(j).connectedto = [i; 0; 5; 2*ks; 0];
%                             end
%                         end
%                     elseif (abs(CellArray(j).position(1)) < 2.2 && CellArray(i).position(1) < CellArray(j).position(1) && otr == 0)
%                         if (CellArray(i).position(1) < CellArray(j).position(1) + 2.2)
%                             if ((CellArray(i).position(2) < CellArray(j).position(2) + .5) || (CellArray(i).position(2) > CellArray(j).position(2) - .5))
%                                 CellArray(i).orientation = CellArray(i).orientation*R(-2);
%                                 CellArray(j).orientation = CellArray(j).orientation*R(2);
%                                 CellArray(i).connectedto = [j; 0; 5; 2*ks; 0];
%                                 CellArray(j).connectedto = [i; 0; 5; 2*ks; 0];
%                             end
%                         end
%                     end
%                 end
%             end
%         end
        


    ks2 = 2*ks;
    if (numel(CellArray) == 4 && Ccount == 0)
        CellArray(2).connectedto = [3, 4; 0, 2; 0, 3; ks2, 9.45; 0, 1];
        CellArray(3).connectedto = [2, 1; 0, 1; 0, 3; ks2, 9.45; 0, 1];
        Ccount = 1;
    elseif(numel(CellArray) == 8 && Ccount == 1)
        CellArray(5).connectedto = [1; 1; 14; 7.45; 1];
        CellArray(7).connectedto = [3; 0; 0; ks2; 0]; %[2, 3; 0, 2; 0, 14; ks2, 4; 0, 0]; %[2; 0; 0; ks2; 0]; %[CellArray(7).connectedto, [2; 0; 0; ks2; 0]]; 
        CellArray(2).connectedto = [8; 0; 0; ks2; 0]; %[CellArray(2).connectedto, [7; 0; 0; ks2; 0]]; 
        CellArray(3).connectedto = [7; 0; 0; ks2; 0]; %[6, 7; 0, 1; 0, 14; ks2, 4; 0, 0]; %[6; 0; 0; ks2; 0];%[CellArray(3).connectedto, [6; 0; 0; ks2; 0]]; 
        CellArray(8).connectedto = [2; 0; 0; ks2; 0];
        %CellArray(6).connectedto = [3, 2, 4; 0, 1, 2; 0, 14, 43; ks2, 3, 2.1; 0, 0, 1]; %[CellArray(6).connectedto, [3; 0; 0; ks2; 0]];
        Ccount = 2;
%         CellArray(5).connectedto = [1; 1; 14; 7.45; 1];
%         CellArray(8).connectedto = [2; 0; 0; ks2; 0]; %[2, 3; 0, 2; 0, 14; ks2, 4; 0, 0]; %[2; 0; 0; ks2; 0]; %[CellArray(7).connectedto, [2; 0; 0; ks2; 0]]; 
%         CellArray(2).connectedto = [8, 7; 0, 2; 0, 14; ks2, ks; 0, 0]; %[CellArray(2).connectedto, [7; 0; 0; ks2; 0]]; 
%         CellArray(3).connectedto = [7; 0; 0; ks2; 0]; %[6, 7; 0, 1; 0, 14; ks2, 4; 0, 0]; %[6; 0; 0; ks2; 0];%[CellArray(3).connectedto, [6; 0; 0; ks2; 0]]; 
%         CellArray(7).connectedto = [3, 2, 4; 0, 1, 2; 0, 14, 43; ks2, 3, 2.1; 0, 0, 1]; %[CellArray(6).connectedto, [3; 0; 0; ks2; 0]];
%         conncount = 2;
    elseif (numel(CellArray) == 16 && cframe == 109) %conncount == 2)
%         CellArray(15).connectedto = [10; 0; 0; ks2; 0]; %[CellArray(15).connectedto, [2; 0; 0; ks2; 0]]; 
%         CellArray(2).connectedto = [7; 0; 0; ks2; 0]; %[CellArray(2).connectedto, [15; 0; 0; ks2; 0]]; 
%         CellArray(7).connectedto = [2; 0; 0; ks2; 0]; %[CellArray(7).connectedto, [10; 0; 0; ks2; 0]]; 
%         CellArray(10).connectedto = [15; 0; 0; ks2; 0]; %[CellArray(10).connectedto, [7; 0; 0; ks2; 0]];
%         CellArray(11).connectedto = [14; 0; 0; ks2; 0]; %[CellArray(11).connectedto, [6; 0; 0; ks2; 0]]; 
%         CellArray(6).connectedto = [3, 12; 0, 0; 0, 0; ks2, ks2; 0, 0]; %[CellArray(6).connectedto, [11; 0; 0; ks2; 0]]; 
%         CellArray(3).connectedto = [6; 0; 0; ks2; 0]; %[CellArray(3).connectedto, [14; 0; 0; ks2; 0]]; 
%         CellArray(14).connectedto = [11; 0; 0; ks2; 0]; %[CellArray(14).connectedto, [3; 0; 0; ks2; 0]];
%         conncount = 3;
        CellArray(10).connectedto = [8; 0; 0; ks2; 0]; %[CellArray(15).connectedto, [2; 0; 0; ks2; 0]]; 
        CellArray(2).connectedto = [13; 0; 0; ks2; 0]; %[CellArray(2).connectedto, [15; 0; 0; ks2; 0]]; 
        CellArray(7).connectedto = [3; 0; 0; ks2; 0]; %[CellArray(7).connectedto, [10; 0; 0; ks2; 0]]; 
        CellArray(8).connectedto = [10; 0; 0; ks2; 0]; %[CellArray(10).connectedto, [7; 0; 0; ks2; 0]];
        CellArray(11).connectedto = [12; 0; 0; ks2; 0]; %[CellArray(11).connectedto, [6; 0; 0; ks2; 0]]; 
        %CellArray(6).connectedto = [3, 12; 0, 0; 0, 0; ks2, ks2; 0, 0]; %[CellArray(6).connectedto, [11; 0; 0; ks2; 0]]; 
        CellArray(12).connectedto = [11, 4; 0, 0; 0, 0; ks2, ks; 0, 0];
        CellArray(3).connectedto = [7; 0; 0; ks2; 0]; %[CellArray(3).connectedto, [14; 0; 0; ks2; 0]]; 
        CellArray(13).connectedto = [2; 0; 0; ks2; 0]; %[CellArray(14).connectedto, [3; 0; 0; ks2; 0]];
        Ccount = 3;
    end
    Ccount    
    
%% plot cells
%if (cframe == 13 || cframe == 19 || cframe == 44 || cframe == 48 || cframe == 68 || cframe == 75 || cframe == 77 || cframe == 108 || cframe == 112)
    %figure
    plotCells(CellArray, k, bt, br, ks, l0, cframe);
%end
    mov(cframe) = getframe(gcf);
    cframe = cframe + 1;

    tc = tc + dt;
   
end
%% write force data to excel
%filename = 'CyanoSimData.xlsx';
%xlswrite(filename, A, 1, 'B2')
%% write the video
vid = VideoWriter('CyanoSim12.avi');
open(vid)
writeVideo(vid, mov)
close(vid)
end
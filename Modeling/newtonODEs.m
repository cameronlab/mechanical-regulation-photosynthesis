function nextX = newtonODEs(t, y, m, mo, bt, br, Fsum, torq)
    nextX = zeros(10, 1);
    nextX(1:2) = y(3:4);
    nextX(3:4) = (1/m)*(Fsum - bt*y(3:4));
    nextX(5:7) = (1/mo)*(torq - br*y(5:7));
    nextX(8:10) = cross(y(5:7), y(8:10));
end

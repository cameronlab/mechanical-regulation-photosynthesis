%From cell 1
cbx12 = cellLine(12).cbxCount(136:192); %only include cbx count from cell 'birth' to division or end
cbx13 = cellLine(13).cbxCount(146:186);
%From cell 2
cbx10 = cellLine(10).cbxCount(35:134);
cbx11 = cellLine(11).cbxCount(35:134);
cbx16 = cellLine(16).cbxCount(134:192);
cbx17 = cellLine(17).cbxCount(134:191);
cbx18 = cellLine(18).cbxCount(135:192);
cbx19 = cellLine(19).cbxCount(135:192);
%From cell 4
cbx8 = cellLine(8).cbxCount(11:94);
cbx9 = cellLine(9).cbxCount(11:128);
cbx14 = cellLine(14).cbxCount(128:192);
cbx15 = cellLine(15).cbxCount(129:181);

%matrix of cbx count values, column 1 = age 0
A = zeros(12, 150);
A(1, 1:length(cbx12)) = cbx12;
A(2, 1:length(cbx13)) = cbx13;
A(3, 1:length(cbx10)) = cbx10;
A(4, 1:length(cbx11)) = cbx11;
A(5, 1:length(cbx16)) = cbx16;
A(6, 1:length(cbx17)) = cbx17;
A(7, 1:length(cbx18)) = cbx18;
A(8, 1:length(cbx19)) = cbx19;
A(9, 1:length(cbx8)) = cbx8;
A(10, 1:length(cbx9)) = cbx9;
A(11, 1:length(cbx14)) = cbx14;
A(12, 1:length(cbx15)) = cbx15;

%iterate through each row and record the index at which (1) the count is
%higher than before and (2) remains constant for more than one frame
growframe = zeros(12, 10); %gives # of frame when cbx was added
grownumber = zeros(12, 10); %gives # of cbx at that frame
for i = 1:12
    n = 1; %counter
    m = 4; %compare to min value
    for j = 1:149
        if (A(i, j) > m) && (A(i, j) < 16) && (A(i, j+1) == A(i, j))
            growframe(i, n) = j;
            grownumber(i, n) = A(i, j);
            n = n+1;
            m = max(grownumber(i, :));
        end
    end
end
growframe
grownumber


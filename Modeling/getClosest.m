function [indy] = getClosest(x, y, contactpt)

    N = size(x, 2);
    A = zeros(N, 3); %[x y z]
    A(:, 1) = x;
    A(:, 2) = y;
    ptCloud = pointCloud(A);
    point = [contactpt, 0];
    [ind, dist] = findNearestNeighbors(ptCloud, point, N);
    indy = ind(1);
end

% function [ind, dist] = getClosest(CellArray, n)
% 
%     N = numel(CellArray);
%     A = zeros(N, 3); %[x y z]
%     for i = 1:N
%         A(i, 1) = CellArray(i).position(1);
%         A(i, 2) = CellArray(i).position(2);
%     end
%     ptCloud = pointCloud(A);
%     point = [CellArray(n).position, 0];
%     [ind, dist] = findNearestNeighbors(ptCloud, point, N);
% end
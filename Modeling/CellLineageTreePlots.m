%Cell Line starting from Cell 1
b1 = [1 2];
d1 = [114; 114; 0].*5;
N1 = {'Cell 12', 'Cell 13', 'Cell 1'};

%Cell Line from Cell 2
b2 = [1 2; 3 4; 5 6];
d2 = [134; 134; 135; 135; 35; 35; 0].*5;
N2 = {'Cell 16', 'Cell 17', 'Cell 18', 'Cell 19', 'Cell 10', 'Cell 11', 'Cell 2'};

%Cell Line from Cell 4
b4 = [1 2; 3 4];
d4 = [128; 128; 11; 11; 0].*5;
N4 = {'Cell 14', 'Cell 15', 'Cell 8', 'Cell 9', 'Cell 4'};

%Plot trees ??How can I plot these on the same axis??
% plot(phytree(b1, d1, N1), 'BranchLabels', 'true', 'TerminalLabels', true)
%plot(phytree(b2, d2, N2),'BranchLabels', 'true', 'TerminalLabels', true)
% plot(phytree(b4, d4, N4),'BranchLabels', 'true', 'TerminalLabels', true)
% title('Cell Lineage'), xlabel('Time (minutes = frame*5)')

%All as one tree
b = [1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16];
d = [128; 128; 134; 134; 135; 135; 114; 114; 11; 11; 35; 35; 0; 0; 0; 0; 0].*5;
N = {'Cell 14', 'Cell 15', 'Cell 16', 'Cell 17', 'Cell 18', 'Cell 19', 'Cell 12', 'Cell 13', 'Cell 8', 'Cell 9', 'Cell 10', 'Cell 11', 'Cell 1', 'Cell 4', 'Cell 2', 'x', 'X'};
plot(phytree(b, d, N),'BranchLabels', 'true', 'TerminalLabels', true)
title('Cell Lineage'), xlabel('Time (minutes = frame*5)')
%ignore x's, they are included as phytrees must be binary
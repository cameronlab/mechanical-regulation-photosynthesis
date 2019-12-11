function g = plotCells(CellArray, k, bt, br, ks, l0, cframe)
    r = CellArray(1).radius; %cell radius
    %rho = 1; %cell density
    colormap gray;
    R = @(th) [cosd(th) -sind(th); sind(th) cosd(th)];
    N = numel(CellArray);
    N1 = 74;
    for i = 1:N
%         if i == 2 %SOLID OBJECT IN FIELD
%             theta = linspace(0, 2*pi, 74);
%             x = 2*cos(theta);
%             y = 2*sin(theta);
%         else
            u = r*CellArray(i).orientation;
            p = CellArray(i).p;
            q = CellArray(i).q;
            R1 = R(90); %rotate CCW
            R2 = R(-5); %rotate CW
            u1 = u*R1;
            u2 = u1;
            P = zeros(2, N1);
            P(:, 1) = q + u1;
            for kk = 2:37
                u2 = u2*R2;
                P(:, kk) = q + u2;
            end
            P(:, 38) = p - u1;
            u3 = -u1;
            for j = 39:74
                u3 = u3*R2;
                P(:, j) = p + u3;
            end
            x = P(1, :);
            y = P(2, :);
        %end
        C = 0.08*ones(N1, 1);
        s1 = CellArray(i).contactpts;
        ls = size(s1, 2);
        fc = 0;
        if(ls > 1)
            for m = 2:ls
                cptind = getClosest(x, y, [s1(1, m) s1(2, m)]);
                xnn = (1:5)';
                if(s1(3, m) >= 15)
                    fc = fc + 5;
                    if cptind < 6
                        C(74) = C(74) + .5;
                        C(6) = 1;
                        C(7:11) = C(7:11) + (1 - xnn./10);
                        C(1:5) = C(1:5) + fliplr(1 - xnn./10);
                    elseif (cptind >= 33 && cptind < 38)
                        C(33) = 1;
                        C(34:38) = C(34:38) + (1 - xnn./10);
                        C(28:32) = C(28:32) + fliplr(1 - xnn./10);
                    elseif (cptind >= 38 && cptind < 43)
                        C(42) = 1;
                        C(43:47) = C(43:47) + (1 - xnn./10);
                        C(37:41) = C(37:41) + fliplr(1 - xnn./10);
                    elseif (cptind >= 70)
                        C(70) = 1;
                        C(71:end) = C(71:end) + (1 - xnn(1:4)./10);
                        C(1) = C(1) + (1 - xnn(5)./10);
                        C(65:69) = C(65:69) + fliplr(1 - xnn./10);
                    else
                        C(cptind) = 1;
                        C((cptind + 1):(cptind + 5)) = C((cptind + 1):(cptind + 5)) + (1 - xnn./10);
                        C((cptind - 5):(cptind - 1)) = C((cptind - 5):(cptind - 1)) + fliplr(1 - xnn./10);
                    end
                elseif((s1(3, m) < 15) && (s1(3, m) >= 10))
                    fc = fc + 4;
                    if cptind < 6
                        C(74) = C(74) + 0.7*0.5;
                        C(6) = 0.7;
                        C(7:11) = C(7:11) + 0.7*(1 - xnn./10);
                        C(1:5) = C(1:5) + fliplr(0.7*(1 - xnn./10));
                    elseif (cptind >= 33 && cptind < 38)
                        C(33) = 0.7;
                        C(34:38) = C(34:38) + 0.7*(1 - xnn./10);
                        C(28:32) = C(28:32) + fliplr(0.7*(1 - xnn./10));
                    elseif (cptind >= 38 && cptind < 43)
                        C(42) = 0.7;
                        C(43:47) = C(43:47) + 0.7*(1 - xnn./10);
                        C(37:41) = C(37:41) + fliplr(0.7*(1 - xnn./10));
                    elseif (cptind >= 70)
                        C(70) = 0.7;
                        C(71:end) = C(71:end) + 0.7*(1 - xnn(1:4)./10);
                        C(1) = C(1) + 0.7*(1 - xnn(5)./10);
                        C(65:69) = C(65:69) + fliplr(0.7*(1 - xnn./10));
                    else
                        C(cptind) = 0.7;
                        C((cptind + 1):(cptind + 5)) = C((cptind + 1):(cptind + 5)) + 0.7*(1 - xnn./10);
                        C((cptind - 5):(cptind - 1)) = C((cptind - 5):(cptind - 1)) + fliplr(0.7*(1 - xnn./10));
                    end
                elseif((s1(3, m) < 10) && (s1(3, m) >= 5))
                    fc = fc + 3;
                    if cptind < 6
                        C(74) = C(74) + 0.5*0.5;
                        C(6) = 0.5;
                        C(7:11) = C(7:11) + 0.5*(1 - xnn./10);
                        C(1:5) = C(1:5) + fliplr(0.5*(1 - xnn./10));
                    elseif (cptind >= 33 && cptind < 38)
                        C(33) = 0.5;
                        C(34:38) = C(34:38) + 0.5*(1 - xnn./10);
                        C(28:32) = C(28:32) + fliplr(0.5*(1 - xnn./10));
                    elseif (cptind >= 38 && cptind < 43)
                        C(42) = 0.5;
                        C(43:47) = C(43:47) + 0.5*(1 - xnn./10);
                        C(37:41) = C(37:41) + fliplr(0.5*(1 - xnn./10));
                    elseif (cptind >= 70)
                        C(70) = 0.5;
                        C(71:end) = C(71:end) + 0.5*(1 - xnn(1:4)./10);
                        C(1) = C(1) + 0.5*(1 - xnn(5)./10);
                        C(65:69) = C(65:69) + fliplr(0.5*(1 - xnn./10));
                    else
                        C(cptind) = 0.5;
                        C((cptind + 1):(cptind + 5)) = C((cptind + 1):(cptind + 5)) + 0.5*(1 - xnn./10);
                        C((cptind - 5):(cptind - 1)) = C((cptind - 5):(cptind - 1)) + fliplr(0.5*(1 - xnn./10));
                    end
                else
                    fc = fc + 2;
                    if cptind < 6
                        C(74) = C(74) + 0.3*0.5;
                        C(6) = 0.3;
                        C(7:11) = C(7:11) + 0.3*(1 - xnn./10);
                        C(1:5) = C(1:5) + fliplr(0.3*(1 - xnn./10));
                    elseif (cptind >= 33 && cptind < 38)
                        C(33) = 0.3;
                        C(34:38) = C(34:38) + 0.3*(1 - xnn./10);
                        C(28:32) = C(28:32) + fliplr(0.3*(1 - xnn./10));
                    elseif (cptind >= 38 && cptind < 43)
                        C(42) = 0.3;
                        C(43:47) = C(43:47) + 0.3*(1 - xnn./10);
                        C(37:41) = C(37:41) + fliplr(0.3*(1 - xnn./10));
                    elseif (cptind >= 70)
                        C(70) = 0.3;
                        C(71:end) = C(71:end) + 0.3*(1 - xnn(1:4)./10);
                        C(1) = C(1) + 0.3*(1 - xnn(5)./10);
                        C(65:69) = C(65:69) + fliplr(0.3*(1 - xnn./10));
                    else
                        C(cptind) = 0.3;
                        C((cptind + 1):(cptind + 5)) = C((cptind + 1):(cptind + 5)) + 0.3*(1 - xnn./10);
                        C((cptind - 5):(cptind - 1)) = C((cptind - 5):(cptind - 1)) + fliplr(0.3*(1 - xnn./10));
                    end
                end
            end
        end
        for j = 1:N1
            abc = C(j);
            if abc > 1
                C(j) = 1;
            end
        end
        h = fill(x, y, C, 'EdgeColor', 'interp', 'LineWidth', 4);
        %set(h,'facealpha', 1)
        if fc >= 10
            set(h, 'FaceColor', [255 255 0]/256) %yellow
        elseif (fc < 10 && fc >= 7)
            set(h, 'FaceColor', [255 127 36]/256) %orange
        elseif (fc < 7 && fc >= 5)
            set(h, 'FaceColor', [238 59 59]/256) %pink
        elseif (fc < 5 && fc >= 2)
            set(h, 'FaceColor', [127 0 255]/256) %purple
        else
            set(h, 'FaceColor', [0 0 0]/256) %blue
        end
        hold on
        %txt1 = num2str(i);
        %text(CellArray(i).position(1),CellArray(i).position(2),txt1)
        %plot(s1(1, :), s1(2, :), 'ok')
    end
    axis([-15 15 -15 15])
    set(gca,'Yticklabel',[]) 
    set(gca,'Xticklabel',[])
    set(gca, 'Color', [0.9 0.9 0.9])
    title(['k = ' (num2str(k)) ', bt = ' num2str(bt) ', br = ' num2str(br) ', ks = ' num2str(ks) ', Spring eq. = ' num2str(l0) ', t = ' (num2str(cframe)) '*dt'])
    caxis([0 1])
    %colorMap = [[0 0 0]/256; [127 0 255]/256; [238 59 59]/256; [255 127 36]/256; [255 255 0]/256];
    %colormap(colorMap)
    colorbar
    hold off 
end


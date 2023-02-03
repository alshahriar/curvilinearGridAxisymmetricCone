function [x] = incremental_dx_dy2(xStart,xEnd,nPx,direction,smallest_dx)

dxV = zeros(nPx - 1,1);
dx = abs(xEnd-xStart)/(nPx-1);

if(smallest_dx>=dx)
    disp(["smallest dx should be more smaller:",sprintf("smaller than %f , cunrrent value %f",dx,smallest_dx)])
    return;
end

smaller_dx = (dx - smallest_dx)/(nPx-2);

x = zeros(nPx,1);

if (direction == 1)
    for i = 1:nPx-1
        dxV(i) = smallest_dx + 2*(nPx-i-1)*smaller_dx;
    end
    x(1) = xStart;
    for i = 1:nPx-1
        x(i+1) = x(i) + dxV(i);
    end
elseif (direction == 0)
    smaller_dx = 0;
    for i = 1:nPx-1
        sum1 = 0;
        for j =  2: nPx-i
            sum1 = sum1 + smaller_dx;
        end
        dx1 = dx - (i-1)*smaller_dx;
        dxV(i) = dx1 + sum1;
        %    oldIncrement = oldIncrement + (nPx-i)*smaller_dx;
    end
    x(1) = xStart;
    for i = 1:nPx-1
        x(i+1) = x(i) + dxV(i);
    end
elseif(direction == -1)
    for i = 1:nPx-1
        dxV(i) = smallest_dx + 2*(i-1)*smaller_dx;
    end
    
    x(1) = xStart;
    for i = 1:nPx-1
        x(i+1) = x(i) + dxV(i);
    end
end

end
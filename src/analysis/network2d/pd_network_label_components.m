function [labels, wrapsX, wrapsY] = pd_network_label_components(phaseMask, validMask, boundary)
%PD_NETWORK_LABEL_COMPONENTS Label 4-connected components and winding loops.

    [rows, columns] = size(phaseMask);
    labels = zeros(rows, columns);
    phaseCount = nnz(phaseMask);
    if phaseCount == 0
        wrapsX = false(0, 1);
        wrapsY = false(0, 1);
        return;
    end
    queue = zeros(phaseCount, 1);
    wrapsX = false(phaseCount, 1);
    wrapsY = false(phaseCount, 1);
    offsetX = nan(rows, columns);
    offsetY = nan(rows, columns);
    component = 0;
    phaseIndices = find(phaseMask);

    for index = 1:numel(phaseIndices)
        seed = phaseIndices(index);
        if labels(seed) ~= 0, continue; end
        component = component + 1;
        head = 1;
        tail = 1;
        queue(1) = seed;
        labels(seed) = component;
        offsetX(seed) = 0;
        offsetY(seed) = 0;
        componentWrapsX = false;
        componentWrapsY = false;
        while head <= tail
            current = queue(head);
            head = head + 1;
            [row, column] = ind2sub([rows, columns], current);
            for direction = 1:4
                [nextRow, nextColumn, exists, shiftX, shiftY] = ...
                    pd_network_neighbor(row, column, direction, rows, columns, boundary);
                if ~exists || ~validMask(nextRow, nextColumn) || ~phaseMask(nextRow, nextColumn)
                    continue;
                end
                nextOffsetX = offsetX(current) + shiftX;
                nextOffsetY = offsetY(current) + shiftY;
                nextIndex = sub2ind([rows, columns], nextRow, nextColumn);
                if labels(nextIndex) == 0
                    tail = tail + 1;
                    queue(tail) = nextIndex;
                    labels(nextIndex) = component;
                    offsetX(nextIndex) = nextOffsetX;
                    offsetY(nextIndex) = nextOffsetY;
                else
                    componentWrapsX = componentWrapsX || nextOffsetX ~= offsetX(nextIndex);
                    componentWrapsY = componentWrapsY || nextOffsetY ~= offsetY(nextIndex);
                end
            end
        end
        wrapsX(component) = componentWrapsX;
        wrapsY(component) = componentWrapsY;
    end
    wrapsX = wrapsX(1:component);
    wrapsY = wrapsY(1:component);
end

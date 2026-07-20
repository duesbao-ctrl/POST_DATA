function topology = pd_network_compute_topology(phaseMask, boundary, beta0, hasMissingCells)
%PD_NETWORK_COMPUTE_TOPOLOGY Betti numbers and Euler characteristic.
% Uses a 4-neighbor digital cell complex for open and wrapped boundaries.

    phaseMask = logical(phaseMask);
    topology = struct('beta0', beta0, 'beta1', 0, 'beta2', 0, ...
        'chi', 0, 'note', conventionNote(boundary, hasMissingCells));
    if ~any(phaseMask(:)), return; end

    vertices = nnz(phaseMask);
    edges = countAdjacencyEdges(phaseMask, boundary);
    faces = countPlaquettes(phaseMask, boundary);
    topology.chi = vertices - edges + faces;
    if strcmp(boundary, 'periodic-xy') && all(phaseMask(:))
        topology.beta2 = 1;
    end
    topology.beta1 = max(0, beta0 + topology.beta2 - topology.chi);
end

function count = countAdjacencyEdges(mask, boundary)
    [rows, columns] = size(mask);
    count = 0;
    for row = 1:rows
        for column = 1:columns
            if ~mask(row, column), continue; end
            [nextRow, nextColumn, exists] = neighbor(row, column, 2, rows, columns, boundary);
            if exists && mask(nextRow, nextColumn), count = count + 1; end
            [nextRow, nextColumn, exists] = neighbor(row, column, 4, rows, columns, boundary);
            if exists && mask(nextRow, nextColumn), count = count + 1; end
        end
    end
end

function count = countPlaquettes(mask, boundary)
    [rows, columns] = size(mask);
    count = 0;
    for row = 1:rows
        for column = 1:columns
            if ~mask(row, column), continue; end
            [sameRow, rightColumn, hasRight] = neighbor(row, column, 2, rows, columns, boundary);
            [downRow, sameColumn, hasDown] = neighbor(row, column, 4, rows, columns, boundary);
            if hasRight && hasDown && mask(sameRow, rightColumn) && ...
                    mask(downRow, sameColumn) && mask(downRow, rightColumn)
                count = count + 1;
            end
        end
    end
end

function [nextRow, nextColumn, exists] = neighbor(row, column, direction, rows, columns, boundary)
    nextRow = row;
    nextColumn = column;
    exists = true;
    switch direction
        case 2
            nextColumn = column + 1;
            if nextColumn > columns
            if pd_network_is_periodic_axis(boundary, 'x'), nextColumn = 1; else, exists = false; end
            end
        case 4
            nextRow = row + 1;
            if nextRow > rows
            if pd_network_is_periodic_axis(boundary, 'y'), nextRow = 1; else, exists = false; end
            end
    end
end

function note = conventionNote(boundary, hasMissingCells)
    if strcmp(boundary, 'open')
        note = ['Topology uses the 4-neighbor digital cell complex. ', ...
            'Under open boundaries beta1 equals the ordinary hole count.'];
    else
        note = ['Topology uses the wrapped 4-neighbor digital cell complex. ', ...
            'Under periodic boundaries beta1 counts all non-trivial 1D loops, ', ...
            'including ordinary enclosed holes and periodic wrapping loops.'];
    end
    if hasMissingCells
        note = [note, ' Missing chunk cells are treated as absent observed cells and may create loops around data gaps.'];
    end
end

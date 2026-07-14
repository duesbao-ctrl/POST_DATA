function [nextRow, nextColumn, exists, shiftX, shiftY] = ...
        pd_network_neighbor(row, column, direction, rows, columns, boundary)
%PD_NETWORK_NEIGHBOR Resolve one 4-neighbor step with periodic displacement.

    nextRow = row;
    nextColumn = column;
    exists = true;
    shiftX = 0;
    shiftY = 0;
    switch direction
        case 1
            nextColumn = column - 1;
            if nextColumn < 1
                if pd_network_is_periodic_axis(boundary, 'x')
                    nextColumn = columns;
                    shiftX = -1;
                else
                    exists = false;
                end
            end
        case 2
            nextColumn = column + 1;
            if nextColumn > columns
                if pd_network_is_periodic_axis(boundary, 'x')
                    nextColumn = 1;
                    shiftX = 1;
                else
                    exists = false;
                end
            end
        case 3
            nextRow = row - 1;
            if nextRow < 1
                if pd_network_is_periodic_axis(boundary, 'y')
                    nextRow = rows;
                    shiftY = -1;
                else
                    exists = false;
                end
            end
        case 4
            nextRow = row + 1;
            if nextRow > rows
                if pd_network_is_periodic_axis(boundary, 'y')
                    nextRow = 1;
                    shiftY = 1;
                else
                    exists = false;
                end
            end
        otherwise
            exists = false;
    end
end

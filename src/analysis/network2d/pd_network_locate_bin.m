function binId = pd_network_locate_bin(coord, edges)
%PD_NETWORK_LOCATE_BIN Locate a coordinate, including the final edge.

    binId = discretize(coord, edges);
    if isnan(binId)
        tol = max(1e-12, 1e-12 * max(abs(edges(end)), 1));
        if abs(coord - edges(end)) <= tol
            binId = numel(edges) - 1;
        else
            binId = 0;
        end
    end
end

function [x, y] = pd_network_scale_coordinates(x, y, coordScale)
%SCALECOORDINATES Apply the configured physical coordinate scale.

    x = x .* coordScale;
    y = y .* coordScale;
end

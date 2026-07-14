function [poreMask, matrixMask] = pd_network_classify_phases(ncountGrid, validMask, thresholdN)
%CLASSIFYPHASES Classify valid cells using the documented Ncount threshold.

    poreMask = validMask & (ncountGrid < thresholdN);
    matrixMask = validMask & ~poreMask;
end

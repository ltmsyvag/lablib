def peaks2binary(nWinPnts, analogData, height=1):
    """
    convert analog data series (peaks) into binary data (1001020...), which is returned.
    bin size is nWinPnts, the tailing points less than nWinPnts is thrown away
    """
    idPeaks, _ = find_peaks(analogData, height=height)
    peakPredicates = np.array([False]*len(analogData))
    for id in idPeaks: peakPredicates[id] = True
    binary = peakPredicates[:-(len(analogData)%nWinPnts)].reshape((-1, nWinPnts))
    binary = binary.sum(axis=1)
    return binary

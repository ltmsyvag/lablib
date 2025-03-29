#%%
from scipy.signal import find_peaks
import numpy as np
def peaks2binary(nWinPnts, analogData, height=1):
    """
    convert analog data series (peaks) into binary data (1001020...), which is returned.
    bin size is nWinPnts, the tailing points less than nWinPnts is thrown away
    """
    idPeaks, _ = find_peaks(analogData, height=height)
    peakPredicates = np.array([False]*len(analogData))
    for id in idPeaks: peakPredicates[id] = True
    remainder = len(analogData)%nWinPnts
    binary = peakPredicates[:-remainder if remainder else None].reshape((-1, nWinPnts))
    binary = binary.sum(axis=1)
    return binary
def peaks2binary2(nWinPnts, analogData, height=1):
    """
    same as peaks2binary, but the binary has the same length as the analogData
    note that nWinPnts does nothing, it's just for the sake of compatibility with peaks2binary, 
    so that fast switching between the two functions is possible
    """
    idPeaks, _ = find_peaks(analogData, height=height)
    peakPredicates = np.array([False]*len(analogData))
    for id in idPeaks: peakPredicates[id] = True
    binary = peakPredicates
    return binary

def color_right_yax(c: type, color : str)->None:
    """
    decorator 专用函数, 将 Axes class 对象右侧 yax 涂成颜色 color
    """
    c.tick_params(colors = color) # tick color 
    c.spines["right"].set_color(color) # edge color
    c.yaxis.label.set_color(color) # label color
def add_Axes_cls_methods(c: type)-> type:
    """
    decorator for plt.Axes type
    """
    c.color_right_yax = color_right_yax
    return c


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.Axes = add_Axes_cls_methods(plt.Axes) # decorate
    
    # fig, ax = plt.subplots()
    # axx = ax.twinx()
    # ax.plot([1,2] , label= 'line1')
    # line, = axx.plot([2,1], color = "r" , label= 'line2') # single item unpacking
    # lcolor = line.get_color() 
    # ax.legend(loc = (0,0))
    # axx.legend(loc = (1,1))
    # axx.set_ylabel("right")
    # axx.color_right_yax("r")

    import inspect
    print(inspect.getsource(plt.Axes.color_right_yax))

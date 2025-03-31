#%%
from scipy.signal import find_peaks
import numpy as np
from matplotlib.axes import Axes
from scipy.fft import fftshift, ifftshift

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


def extend_Axes_methods(c: type[Axes])-> type[Axes]: # 所有的类的 type 都是 `type`, 但是由于下面的函数 type annotation 中要写具体的 class `Axes`, 因此, 反正都要 import 这个 `Axes` 对象, 就在 decorator 输入的 type annotation 中也 explicitly 写出 `type[Axes]` 吧. 否则直接写 `type` 也无不可
    """
    为 matplotlib Axes 实例(e.g. ax) 追加方法
    """
    def color_right_yax(self: Axes, color: str)->None:
        """
        decorator 专用函数, 将 Axes 对象右侧 yax 涂成颜色 color
        """
        self.tick_params(colors = color) # tick color 
        self.spines["right"].set_color(color) # edge color
        self.yaxis.label.set_color(color) # label color
    c.color_right_yax = color_right_yax # 追加一个实例方法
    return c
    
def fdata_keep_n_lowfreq_pnts(fdata, nPositive_freq_pnts_kept: int)->np.ndarray:
    """
    一个简单的频域高频成分截断 filter, 可以用于对任何数据序列的 smoothing (不需要是时域数据)
    fft(data) 后得到的 fdata 频率序列有两种情况:
    0 1 2 3 -4 -3 -2 -1
    0 1 2 3 -3 -2 -1
    fftshift(fdata) 后有两种情况:
    -4 -3 -2 -1 0 1 2 3
       -3 -2 -1 0 1 2 3
    其中 [1 2 3] 被称为正频率成分, 不包含 0 (DC), 以上两种情况中, fdata_keep_n_lowfreq_pnts 均为 3
    如果保留两个正频率成分, 那么上述两种情况下均得到 (其他成分设为 0):
    -2 -1 0 1 2
    如果保留 0 个正频率成分, 那么得到:
    0 
    如果保留所有的点那么, 两种情况下最终均会返回原始成分 (那奎斯特频率 -4 如果有, 会保留):
    -4 -3 -2 -1 0 1 2 3
       -3 -2 -1 0 1 2 3
    """
    assert isinstance(nPositive_freq_pnts_kept, int) and (nPositive_freq_pnts_kept >=0), "需要保留的点数是非负整数"
    nPositive_freq_pnts = int(len(fdata)/2-0.1) # 永远返回正确的 positive frequency components 数量, 见 docstring
    assert nPositive_freq_pnts_kept<=nPositive_freq_pnts, "正频率成分数量上限 ≈ signal 点数的一半!"
    sfdata = fftshift(fdata)
    nPositive_freq_pnts_thrown = nPositive_freq_pnts - nPositive_freq_pnts_kept
    if nPositive_freq_pnts_thrown: # when this number is 0, do nothing
        sfdata[-nPositive_freq_pnts_thrown:] = 0 # set positive high freq components to zero
    if nPositive_freq_pnts_thrown: # 如果根本没有丢弃正频率成分, 那么也不必要动负频率成分
        nNegative_freq_pnts_thrown = nPositive_freq_pnts_thrown
        if len(fdata)%2 == 0: 
            nNegative_freq_pnts_thrown += 1
        sfdata[:nNegative_freq_pnts_thrown] = 0 # set negative high freq components to zero
    fdata_filtered = ifftshift(sfdata)
    return fdata_filtered
if __name__ == "__main__":
    
    import matplotlib.pyplot as plt
    plt.Axes = extend_Axes_methods(plt.Axes) # decorate

    fig, ax = plt.subplots()
    axx = ax.twinx()
    # ax.plot([1,2] , label= 'line1')
    # line, = axx.plot([2,1], color = "r" , label= 'line2') # single item unpacking
    # lcolor = line.get_color() 
    # ax.legend(loc = (0,0))
    # axx.legend(loc = (1,1))
    # axx.set_ylabel("right")
    axx.color_right_yax("r")

    import inspect
    print(inspect.getsource(plt.Axes.color_right_yax))

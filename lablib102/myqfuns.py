from scipy.constants import h, k
import numpy as np
def n_th(freq, T):
    """
    return thermal photons. freq in Hz, T in K
    """
    return 1/(np.exp(h*freq/k/T)-1)

def numberToBase(n, b):
    "labelKet 的辅助函数, 源自 https://stackoverflow.com/questions/2267362/how-to-convert-an-integer-to-a-string-in-any-base"
    if n == 0:
        return [0]
    digits = []
    while n:
        digits.append(int(n % b))
        n //= b
    return digits[::-1]

def labelKet(N, ket):
    assert ket.isket, "只能 label ket 对象"
    nParts = len(ket.dims[0])
    print("amplitude\tlabel")
    for id, amplitude in enumerate(ket):
        label = numberToBase(id, N)
        label = [0]*(nParts-len(label)) + label
        print(amplitude, label, sep = "\t")

from typing import Sequence
from qutip import Qobj, projection
def Omega_couple(N: int, Omega, ij: Sequence[int]) -> Qobj:
    """
    `N` is space size, `Omega` is the usual Rabi frequency Ω
    """
    i,j = ij
    assert i != j, "Omega coupling cannot be diagonal"
    return (projection(N,i,j)+projection(N,j,i))*Omega

from qutip import Qobj, projection
import numpy as np
def make_c_ops(N: int, decayDict: dict=None, dephaseDict: dict=None) -> list:
    c_ops = []
    if decayDict:
        for ij, Gammaij in decayDict.items():
            if "u" in ij: # if it is an unphysical channel, 需要 "u" 是为了防止重名, e.g. "51" 和 "51u" 均可以出现
                i,j = (int(e) for e in ij[:2])
            else:
                i,j = (int(e) for e in ij)
            L = np.sqrt(Gammaij)*projection(N,j,i)
            c_ops.append(L)
    if dephaseDict:
        for ij, gammaij in dephaseDict.items():
            i,j = (int(e) for e in ij)
            L = np.sqrt(gammaij/2)*(projection(N,i,i) - projection(N,j,j))
            c_ops.append(L)
    return c_ops

import re
from arc import *

def make_decayRateDict(termList: list[str], 
                       additional_channels = [], 
                       unphysicalChannels = [], 
                       temperature=400):
    """
    - additional channels 中用 string 写明非禁戒的自发辐射 decay channel, 
        比如 `"51"`, `"40"` (数字顺序无所谓 "40", 和 "04" 等价)
        若写入禁戒的 channels 是无效的, 不会报错
    - unphysicalChannels 中用 string-val tuple 写明需要添加的 decay channel 及其 rate, 
        无视禁戒法则, 数字顺序重要, `("51",3)` 和 `("15",3)` 是不同的.
        它们分别描述 5 -> 1 和 1 -> 5 两种不同的 decay
    """
    pattern = r"([1-9]?[0-9])[ ]*([A-Z])[ ]*([1-9])/([2])"
    nljList, decayRate, udStr = [], dict() , ""
    atom = Rubidium85()
    letter_l_dict = {"S" : 0, "P" : 1, "D" : 2, "F" : 3, "G" : 4, "H" : 5, "I" : 6}
    
    for term in termList:
        match = re.fullmatch(pattern, term)
        assert match, "term 格式输入错误(例: 5S1/2)"
        jmanif = int(match.group(1)), letter_l_dict[match.group(2)], int(match.group(3))/int(match.group(4))
        nljList.append(jmanif)

    for id, (jmanif1, jmanif2) in enumerate(zip(nljList, nljList[1:])):
        freq12 = atom.getTransitionFrequency(*jmanif1, *jmanif2)
        if freq12>0: thisUD = "u"
        else: thisUD = "d"
        udStr+=thisUD
        decayRate[f"{id}{id+1}"] = atom.getTransitionRate(*jmanif1, *jmanif2, temperature=temperature)
        decayRate[f"{id+1}{id}"] = atom.getTransitionRate(*jmanif2, *jmanif1, temperature=temperature)
    
    for ij in additional_channels:
        i,j = (int(e) for e in ij)
        decayRate[ij] = atom.getTransitionRate(*nljList[i], *nljList[j], temperature=temperature)
        ji = ij[1]+ij[0]
        decayRate[ji] = atom.getTransitionRate(*nljList[j], *nljList[i], temperature=temperature)
    for ij, rate in unphysicalChannels:
        i,j = (int(e) for e in ij)
        decayRate[ij+"u"] = rate*1e6 # 加 "u"(unphysical) 防止 e.g. "51u" 和 "51" 重名; rate 转化为 Hz
    for key, val in decayRate.items():
        decayRate[key] = val*1e-6 # 转化为 MHz
        # decayRate[key] = round(val*1e-6, 6) # might make calculation faster by eliminating too many significant digits
    return decayRate, udStr

from qutip import steadystate, qdiags
def steadystateMWM(termList: list[str],
                   Ωlist_no2pi: list, 
                   DeltaList_no2pi: list,
                   nonNNchannels=[], 
                   unphysicalChannels=[],
                   laserDephList = [], # laser dephasings
                   temperature=400):
    """
    接收能级 term symbols, 拉比频率(MHz, 无 2π), 和 detuning 频率 (MHz, 无 2π), 
    返回 Bloch 模型稳态密度矩阵.
    本函数考虑 Rb87 自发辐射速率和所涉及有限能级(e.g. 6WM 就只有6个能级)之间的 300 K 黑体辐射速率.
    如果有非近邻能级之间的自发辐射, 需要手动添加. 
    nonNNchannels 的意思是 non-NearestNeighbor-channels.
    用 dephList 添加 dephasing, 其中元素是 tuple, 格式是 ("51", 3) (能级 5,1 之间的 激光的有 3 MHz 的 dephasing)
    """
    Ωlist = [2*np.pi*e for e in Ωlist_no2pi]
    DeltaList = [2*np.pi*e for e in DeltaList_no2pi]
    decayRate, udStr = make_decayRateDict(
        termList, nonNNchannels, unphysicalChannels, temperature)
    N, Ω, DeltaSumList =len(termList), dict(), [0]

    for i,j, Ωij in zip(range(N),range(1,N), Ωlist): Ω[f"{i}{j}"] = Ωij

    DeltaSum = 0
    for ud, Deltaij in zip(udStr, DeltaList):
        if ud == "u": DeltaSum += Deltaij
        else: DeltaSum -= Deltaij
        DeltaSumList.append(DeltaSum)

    H = qdiags(DeltaSumList)
    for ij, Ωij in Ω.items():
        i,j = (int(e) for e in ij)
        H += Omega_couple(N, Ωij/2, [i,j])
    # make dephasing dict
    dephDict = dict()
    for ij, rate in laserDephList:
        dephDict[ij] = rate
    c_ops = make_c_ops(N, decayRate, dephDict)
    return steadystate(H, c_ops)

def steadystateLDY(Ωlist_no2pi: list, DeltaList_no2pi: list, temperature = 400):
    return steadystateMWM(
        ["5S1/2", "5P1/2", "34S1/2", "34P3/2", "33D5/2", "5P3/2"],
        Ωlist_no2pi=Ωlist_no2pi, 
        DeltaList_no2pi=DeltaList_no2pi, 
        nonNNchannels=["50", "30", "25", "41"], 
        temperature=temperature)

def steadystateZXL(Ωlist_no2pi: list, DeltaList_no2pi: list, temperature=400):
    return steadystateMWM(
        ["5S1/2", "5P1/2", "25S1/2", "24P3/2", "23D5/2", "5P3/2"],
        Ωlist_no2pi=Ωlist_no2pi, 
        DeltaList_no2pi=DeltaList_no2pi, 
        nonNNchannels=["50", "30", "41","25"], 
        temperature=temperature)

if __name__ == "__main__":
    termList = ["5S1/2", "5P1/2", "34S1/2", "34P3/2", "33D5/2", "5P3/2"]
    termList = ["5S1/2", "5P1/2", "25S1/2", "24P3/2", "23D5/2", "5P3/2"]
    decayRate, udStr = make_decayRateDict(termList, ["50", "30"])
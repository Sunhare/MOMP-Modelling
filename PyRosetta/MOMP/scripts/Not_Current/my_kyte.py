# import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

def find_contiguous_colors(colors):
    # finds the continuous segments of colors and returns those segments
    segs = []
    curr_seg = []
    prev_color = ''
    for c in colors:
        if c == prev_color or prev_color == '':
            curr_seg.append(c)
        else:
            segs.append(curr_seg)
            curr_seg = []
            curr_seg.append(c)
        prev_color = c
    segs.append(curr_seg) # the final one
    return segs
 
def plot_multicolored_lines(x,y,colors):
    segments = find_contiguous_colors(colors)
    plt.figure()
    start= 0
    for seg in segments:
        end = start + len(seg)
        l, = plt.gca().plot(x[start:end],y[start:end],lw=2,c=seg[0]) 
        start = end

def in_range(i, max_len):
	if i > 0 and i < max_len:
		return True
	return False

D = {'A': 1.8, 'C': 2.5, 'D': -3.50, 'E': -3.50,
	 'F': 2.8, 'G':-0.4, 'H': -3.20, 'I': 4.50,
	 'K': -3.90, 'L': 3.80, 'M': 1.90, 'N': -3.50,
	 'P': -1.60, 'Q': -3.50, 'R': -4.50, 'S':-0.80,
	 'T': -0.7, 'V': 4.20, 'W': -0.90, 'Y': -1.30}

# myseq = 'ACDEFGHI'
myseq = 'LPVGNPAEPSLMIDGILWEGFGGDPCDPCTTWCDAISLRLGYYGDFVFDRVLKTDVNKQFEMGPVPTTTDTDAAADITTSTPRENPAYGKHMQDAEMFTNAAYMALNIWDRFDVFCTLGATSGYLKGNSASFNLVGLFGDGVANAANAIATVAADSLPNVSLSQAVVELYTDTAFAWSVGARAALWECGCATLGASFQYAQSKPKVEELNVLCNAAQFTINKPKGYVGKEFPLALTAGTDSATDTKDASIDYHEWQASLALSYRLNMFTPYIGVKWSRASFDADTIRIAQPKLAEAILDVTTWNPTIAGAGTIADGTGAAATANGLADTLQIVSLQLNKMKSRKSCGLAIGTTIVDADKYAVTVETRLIDERAAHVNAQFRF'
# myseq = 'EVKLSGDARMGVMYNGDDWNFSSRSRVLFTMSGTTDSGLEFGASFKAHESVGAETGEDGTVFLSGAFGKIEMGDALGASEALFGDLYEVGYTDLDDRGGNDIPYLTGDERLTAEDNPVLLYTYSAGAFSVAASMSDGKVGETSEDDAQEMAVAAAYTFGNYTVGLGYEKIDSPDTALMADMEQLELAAIAKFGATNVKAYYADGELDRDFARAVFDLTPVAAAATAVDHKAYGLSVDSTFGATTVGGYVQVLDIDTIDDVTYYGLGASYDLGGGASIVGGIADNDLPNSDMVADLGVKFKF'
x = [i for i in range(len(myseq))]


#Window size of 2 means skip to the right 4 and skip left 4
def BB_kyte_doolittle(seq, window):
	global D
	n = len(seq)
	# hphob = [0 for i in range(n)]
	hphob = np.zeros(n)
	for i in range(n):
		hphob[i] += D[seq[i]]
		for j in range(window+1):

			ind_pos = i+(j*2)
			ind_neg = i-(j*2)
			
			if in_range(ind_pos, n) and j != 0:#Don't double count
				hphob[i] += D[seq[ind_pos]]
			if in_range(ind_neg, n) and j != 0: 
				hphob[i] += D[seq[ind_neg]]
	return hphob

hphob = np.zeros(len(myseq))
# for z in range(2):
# 	hphob += BB_kyte_doolittle(myseq, z+1)

# hphob0 = BB_kyte_doolittle(myseq, 0)
# hphob1 = BB_kyte_doolittle(myseq, 1)
hphob = BB_kyte_doolittle(myseq, 1)
# hphob3 = BB_kyte_doolittle(myseq, 3)

colors = ['blue']*len(myseq) #Color everything blue
#Variable Domain Colors for C. Muridarum
vd1_start=86;vd1_end=106;vd2_start=162;vd2_end=187;vd3_start=251;vd3_end=264;vd4_start=315;vd4_end=349;
colors[vd1_start:vd1_end] = ['orange']*(vd1_end-vd1_start)
colors[vd2_start:vd2_end] = ['green']*(vd2_end-vd2_start)
colors[vd3_start:vd3_end] = ['cyan']*(vd3_end-vd3_start)
colors[vd4_start:vd4_end] = ['red']*(vd4_end-vd4_start)
 
plot_multicolored_lines(x,hphob,colors)
# plt.plot(x, hphob0, c='red')
# plt.plot(x, hphob1, c='blue')
# plt.plot(x, hphob2, c='green')
# plt.plot(x, hphob, c='blue')
plt.axhline(c='black')
plt.xlabel("Sequence Position")
plt.ylabel("Hydrophobicity (Kyte Doolittle)")
plt.title("Kyte Doolittle Plot for C. Muridarum")
plt.show()


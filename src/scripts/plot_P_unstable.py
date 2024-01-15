import glob
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import paths
plt.style.use(paths.scripts / "matplotlibrc")

def read_data(file_path, N_skip=None, convert = True):
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i==0:
                col = line.split()
                break
    binary_file = str(file_path[:-4]) + ".npy"
    if os.path.isfile(binary_file):
        # try reading binary
        src = np.load(binary_file)
    else:
        print("binary data version does not exist, it may take a while...")
        src = np.genfromtxt(file_path, skip_header=1)
        if convert:
            # save binary version for speedup
            print("saving binary output file...")
            np.save(binary_file, src)
    if N_skip:
        return src[::N_skip, :], col
    else:
        return src, col


def select_fixed_merger_mass_loss(src, col):
    f_dm = src[:, col.index("f_dm")]
    values = np.unique(f_dm)
    list_src = []
    for v in values:
        ind = f_dm == v
        list_src.append(src[ind, :])
    return list_src, values


def get_binary_name_from_file(file_path):
    name = file_path.split('/')[-1].split("_")[0]
    return name

def plot_P_unstable(file_path, ax, N_skip=None):
    print(file_path)
    bin_name = get_binary_name_from_file(file_path)
    # ax.text(0.01, 0.5, bin_name, fontsize=30, transform=ax.transAxes, va="center", ha="left")
    ax.set_title(bin_name, fontsize=30)
    src, col = read_data(file_path, N_skip=N_skip)
    # get list of f_dm tried
    list_src, f_dm_values = select_fixed_merger_mass_loss(src, col)
    linestyles = ["-", "--"]
    for i, src_f_dm in enumerate(list_src):
        P_unstable = src_f_dm[:, col.index("P_unstable")]
        # ax.hist(P_unstable, range=(0,1), color=colors[i], histtype='step', bins=100, density=True)
        ax.hist(P_unstable, range=(0,1), ls=linestyles[i], # color=color,
                histtype='step', bins=100, cumulative=-1, density=True,
                label=r"$f_\Delta="+f"{f_dm_values[i]:.1f}"+r"$", lw=3)
    # ax.set_yscale('log')

def make_Punstable_plot(N_skip=None, fig_name=None):
    fig = plt.figure(figsize=(10,10))
    gs = gridspec.GridSpec(130, 100)
    ax1 = fig.add_subplot(gs[:30, :])
    ax2 = fig.add_subplot(gs[50:80, :])
    ax3 = fig.add_subplot(gs[100:, :])
    axes = [ax1,ax2,ax3]
    files = glob.glob(str(paths.data)+'/triples/*_triples.txt')
    print(files)
    for i, file_path in enumerate(files):
        ax = axes[i]
        plot_P_unstable(file_path, ax, N_skip=N_skip)
    for ax in axes:
        ax.set_xlim(0,1)
    ax3.set_xlabel(r"$P$(triple unstable at ZAMS)",fontsize=25)
    ax2.set_ylabel(r"Fraction "+r"$P\geq x$",fontsize=25)
    ax3.legend(loc="lower left", ncol=2, columnspacing=0.5)
    if fig_name:
        plt.savefig(fig_name)
        # save static too
        #plt.savefig(paths.static/"P_unstable.pdf")
        #print(fig_name)
    else:
        plt.show()


if __name__ == "__main__":
    make_Punstable_plot(N_skip=None, fig_name=paths.figures/"P_unstable.pdf")

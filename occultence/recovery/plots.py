from ..imports import *
from matplotlib.patches import Patch

def twod_plot(x, y, z, x_bins, y_bins, zlabel, statistic, vlims=[], nsig='%0.2f', ylog=False, xlog=False, smooth=False,
              addtext=False, svname=""):
    # cmap = plt.get_cmap("Blues")
    cmap = plt.get_cmap("viridis")
    # fig, (ax2) = plt.subplots(1, 1, figsize=set_size(width, fraction=1.2, ratio=0.8))
    fig, (ax2) = plt.subplots(1, 1, figsize=set_size(width, fraction=1.7, ratio=0.5))
    fig.canvas.draw()
    ret = binned_statistic_2d(x, y, z * 100, statistic=statistic, bins=[x_bins, y_bins])
    ret2 = binned_statistic_2d(x, y, z * 100, statistic="count", bins=[x_bins, y_bins])
    ret3 = binned_statistic_2d(x, y, z * 100, statistic="sum", bins=[x_bins, y_bins])
    err = ret2.statistic.T.copy()

    # im = ax2.pcolormesh(spt_bins, energy_bins, sens_z, cmap=cmap)
    if not smooth:
        if len(vlims) == 0:
            im = ax2.pcolormesh(ret.x_edge, ret.y_edge, ret.statistic.T, cmap=cmap)
        else:
            im = ax2.pcolormesh(ret.x_edge, ret.y_edge, ret.statistic.T, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
        cbar = fig.colorbar(im, ax=ax2)
        plt.ylabel("Injected Radius (R$_{\oplus}$)")
        plt.xlabel("Injected Period (days)")
        cbar.set_label(zlabel, rotation=90, labelpad=20)
        ax2.set_facecolor('gray')

        if addtext:
            for i in range(len(x_bins) - 1):
                for j in range(len(y_bins) - 1):
                    N = ret2.statistic.T[j][i]
                    n_success = ret3.statistic.T[j][i]
                    error_in_frac = (1 / N) * np.sqrt(n_success * (1 + (n_success / N)))
                    err[j][i] = error_in_frac
                    if not np.isnan(ret.statistic.T[j][i]):
                        # print(spt_bins[i],energy_bins[j],sens_z.values[j][i])
                        txt = plt.text(x_bins[i] + 0.5 * (x_bins[i + 1] - x_bins[i]),
                                       y_bins[j] + 0.5 * (y_bins[j + 1] - y_bins[j]),
                                       nsig % ret.statistic.T[j][i],
                                       color="k", ha="center", va="center", fontweight="bold", fontsize=10)
                        txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w', alpha=0.7)])

                        # txt.set_path_effects([PathEffects.withSimplePatchShadow(linewidth=4, foreground='w')])
        # plt.plot(1.51087, 1.144, marker='x', color='k', label="TRAPPIST-1b")
        plt.plot(2.7299024, 1.320, marker='*', color='r', markersize=10, label="SPECULOOS-2b")
        plt.plot(8.457463, 1.366, marker='*', color='#07b2e6', markersize=10, label="SPECULOOS-2c")
        if xlog:
            plt.xscale('log')
        if ylog:
            plt.yscale('log')

    else:
        from scipy.ndimage import gaussian_filter
        z_smooth = gaussian_filter(ret.statistic.T, sigma=1.5)
        if len(vlims) == 0:
            im = ax2.pcolormesh(ret.x_edge, ret.y_edge, z_smooth, cmap=cmap)
        else:
            im = ax2.pcolormesh(ret.x_edge, ret.y_edge, z_smooth, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
        cbar = fig.colorbar(im, ax=ax2)
        plt.ylabel("Radius of Injected Planet (R$_{Earth}$)")
        plt.xlabel("Period (d)")
        cbar.set_label(zlabel, rotation=270, labelpad=20)
        ax2.set_facecolor('gray')
        # plt.plot(1.51087, 1.144, marker='x', color='k', label="TRAPPIST-1b")
        if xlog:
            plt.xscale('log')
        if ylog:
            plt.yscale('log')

    if svname == "":
        plt.show()
    else:
        plt.savefig(svname + ".pdf", bbox_inches='tight')
    plt.close()


def twod_plot_with_errors(x, y, z, x_bins, y_bins):
    cmap = plt.get_cmap("Blues")
    # fig, (ax2) = plt.subplots(1, 1, figsize=set_size(width, fraction=1.2, ratio=0.8))
    fig, (ax2) = plt.subplots(1, 1, figsize=(10, 8))
    fig.canvas.draw()
    ret_mean = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])
    ret = binned_statistic_2d(x, y, z, statistic='count', bins=[x_bins, y_bins])
    # ret2 = binned_statistic_2d(x, y, z, statistic='std', bins=[x_bins, y_bins])
    ret3 = binned_statistic_2d(x, y, z, statistic='sum', bins=[x_bins, y_bins])
    # err = [std / np.sqrt(n) for std, n in zip(ret2.statistic.T, ret.statistic.T)]
    # im = ax2.pcolormesh(spt_bins, energy_bins, sens_z, cmap=cmap)
    im = ax2.pcolormesh(ret_mean.x_edge, ret_mean.y_edge, ret_mean.statistic.T, cmap=cmap)
    cbar = fig.colorbar(im, ax=ax2)
    plt.ylabel("Radius of Injected Planet (R$_{Earth}$)")
    plt.xlabel("Period (d)")
    cbar.set_label("Fraction of Recovered Planets", rotation=270, labelpad=20)
    ax2.set_facecolor('gray')

    err = ret.statistic.T.copy()
    numtests = ret.statistic.T.copy()
    numsuccess = ret3.statistic.T.copy()

    for i in range(len(x_bins) - 1):
        for j in range(len(y_bins) - 1):
            if not np.isnan(err[j][i]):
                N = ret.statistic.T[j][i]
                n_success = ret3.statistic.T[j][i]
                error_in_frac = (1 / N) * np.sqrt(n_success * (1 + (n_success / N)))
                err[j][i] = error_in_frac
                # print(spt_bins[i],energy_bins[j],sens_z.values[j][i])
                txt = plt.text(x_bins[i] + 0.5 * (x_bins[i + 1] - x_bins[i]),
                               y_bins[j] + 0.5 * (y_bins[j + 1] - y_bins[j]),
                               "%.3f+-%.3f" % (ret_mean.statistic.T[j][i], err[j][i]),
                               color="k", ha="center", va="center", fontweight="bold", fontsize=9)
                txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w', alpha=0.7)])
                # txt.set_path_effects([PathEffects.withSimplePatchShadow(linewidth=4, foreground='w')])

    plt.xscale('log')
    plt.show()
    plt.close()


def plot_expected_planets(x, y, z, x_bins, y_bins, numstars, yscale, xscale, svname=""):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    cmap = plt.get_cmap("Blues")
    fig, (ax2) = plt.subplots(1, 1, figsize=(6, 6))
    fig.canvas.draw()
    ret = binned_statistic_2d(x, y, z, statistic='mean', bins=[x_bins, y_bins])
    im = ax2.pcolormesh(ret.x_edge, ret.y_edge, ret.statistic.T * numstars, cmap=cmap)
    cbar = fig.colorbar(im, ax=ax2)
    plt.ylabel("Radius of Injected Planet (R$_{Earth}$)")
    plt.xlabel("Period (d)")
    cbar.set_label("Number of Expected Planets", rotation=270, labelpad=20)
    ax2.set_facecolor('gray')

    for i in range(len(x_bins) - 1):
        for j in range(len(y_bins) - 1):
            if not np.isnan(ret.statistic.T[j][i]):
                # print(spt_bins[i],energy_bins[j],sens_z.values[j][i])
                txt = plt.text(x_bins[i] + 0.5 * (x_bins[i + 1] - x_bins[i]),
                               y_bins[j] + 0.5 * (y_bins[j + 1] - y_bins[j]),
                               "%.2f" % (ret.statistic.T[j][i] * numstars),
                               color="k", ha="center", va="center", fontweight="bold", fontsize=10)
                txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w', alpha=0.7)])
                # txt.set_path_effects([PathEffects.withSimplePatchShadow(linewidth=4, foreground='w')])
    if xscale == "log":
        plt.xscale('log')
    if yscale == "log":
        plt.yscale('log')
    if svname == "":
        plt.show()
    else:
        plt.savefig(svname + ".pdf", bbox_inches='tight')
    plt.close()


def plot_transitparams(x, y, z, xlabel="Radius", ylabel="SNR", zlabel="Detected?", ylims=[], yscale='uniform',
                       xscale='uniform', add_points={}, svname="", figsize=(6,6)):
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    markers = ['*', "x", "D", "^", 'o']
    plt.figure(figsize=figsize)
    plt.scatter(x[z == 0], y[z == 0], s=1.2, color="#2694c7", alpha=0.7)
    plt.scatter(x[z == 1], y[z == 1], s=1.2, color="#f08f18", alpha=0.7)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if yscale == "log":
        plt.yscale('log')
    if xscale == "log":
        plt.xscale('log')
    if len(ylims) > 0:
        plt.ylim(ylims[0], ylims[1])
    if len(add_points) > 0:
        mc = 0
        for k, v in add_points.items():
            count = 0
            for vi in v:
                if count == 0:
                    plt.plot(vi[0], vi[1], label=k, marker=markers[mc], color='k', linestyle="None")
                else:
                    plt.plot(vi[0], vi[1], marker=markers[mc], color='k', linestyle="None")
                count = count + 1
            mc = mc + 1
    lgd = plt.legend()
    add_patch(lgd)
    if svname == "":
        plt.show()
    else:
        plt.savefig(svname + ".png", bbox_inches='tight', dpi=800)
    plt.close()

def add_patch(legend):
    ax = legend.axes

    handles, labels = ax.get_legend_handles_labels()
    handles.append(Patch(facecolor="#2694c7", edgecolor="#2694c7"))
    labels.append("Not Recovered")
    handles.append(Patch(facecolor="#f08f18", edgecolor="#f08f18"))
    labels.append("Recovered")

    legend._legend_box = None
    legend._init_legend_box(handles, labels)
    legend._set_loc(legend._loc)
    legend.set_title(legend.get_title().get_text())
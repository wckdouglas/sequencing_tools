import re
from builtins import range, zip
from collections import OrderedDict

import numpy as np
from pandas import Series
from sequencing_tools.utils import SeqUtilsError


def douglas_palette():
    """
    Automatic set color if seaborn is installed, otherwise return list of colors
    Example::

        colors = douglas_palette()

        ax=plt.subplot();
        for i in range(14):
            ax.plot(np.arange(10),np.arange(10) + i,label=i, color = colors[i])
        ax.legend()
    """
    colors = [
        "#B98476",
        "#6297B6",
        "#BF7F8E",
        "#8E955A",
        "#6F94B9",
        "#579F79",
        "#3C9F9D",
        "#EDAEB5",
        "#B0BDEA",
        "#D8BA90",
        "#7BCBD5",
        "#87CDA9",
        "#B1C68D",
        "#E7C039",
    ]
    return colors


def maximum_palette():
    """
    modified from:
    https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/

    Example::

        colors = maximum_palette()

        ax=plt.subplot();
        for i in range(22):
            ax.plot(np.arange(10),np.arange(10) + i,label=i, color = colors[i])
        ax.legend()
    """
    colors = [
        "#f58231",
        "#e6194b",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#911eb4",
        "#03A8FB",
        "#F8BF6C",
        "#CAF5CB",
        "#fabebe",
        "#008080",
        "#e6beff",
        "#9a6324",
        "#fffac8",
        "#800000",
        "#aaffc3",
        "#808000",
        "#ffd8b1",
        "#000075",
        "#808080",
        "#ffffff",
        "#000000",
    ]
    return colors


def simpsons_palette():
    """
    A palette from ggsci R package
    https://github.com/road2stat/ggsci/blob/master/data-raw/data-generator.R

    Example::

        colors = simpsons_palette()

        ax=plt.subplot();
        for i in range(16):
            ax.plot(np.arange(10),np.arange(10) + i,label=i, color = colors[i])
        ax.legend()
    """
    colors = [
        "#FED439",
        "#709AE1",
        "#8A9197",
        "#D2AF81",
        "#FD7446",
        "#D5E4A2",
        "#197EC0",
        "#F05C3B",
        "#46732E",
        "#71D0F5",
        "#370335",
        "#075149",
        "#C80813",
        "#91331F",
        "#1A9993",
        "#FD8CC1",
    ]
    return colors


def okabeito_palette():
    """
    Color palette proposed by Okabe and Ito
    copy from colorboindr R package
    https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R

    Example::

        colors = okabeito_palette()

        ax=plt.subplot();
        for i in range(8):
            ax.plot(np.arange(10),np.arange(10) + i,label=i, color = colors[i])
        ax.legend()
    """
    colors = [
        "#E69F00",
        "#56B4E9",
        "#009E73",
        "#F0E442",
        "#0072B2",
        "#D55E00",
        "#CC79A7",
        "#999999",
    ]
    return colors


def cor_plot(plot_df, fig, diagonal_line=True, method="pearson", **kwargs):
    """
    given a data frame with columns storing data for each sample,
    output a matplotlib figure object as correlation plots.

    usage: cor_plot(plot_df, fig=fig, diagonal_line = True, method = 'pearson')

    Args:
        plot_df: a pandas dataframe
        fig: matplotlib figure object (optional)
        diagonal_line: Boolean controlling if a diagonal line should be drawn
        method: correlation method, [spearman or pearson]

    Returns:
        matplotlib.pyplot.figure: matplotlib figure object
    """

    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("white")
    assert method in ["spearman", "pearson"], "Wrong correlation method"
    box_size = len(plot_df.columns)
    for row in range(box_size):
        for col in range(box_size):
            ax = fig.add_subplot(box_size, box_size, row * box_size + col + 1)

            #### Right side plots
            if col < row:
                ax.scatter(plot_df.iloc[:, col], plot_df.iloc[:, row], **kwargs)

                if diagonal_line:
                    plot_data = plot_df.iloc[:, [col, row]]
                    maxima = plot_data.max()
                    minima = plot_data.min()
                    ax.plot([minima, maxima], [minima, maxima], color="red")

            #### Diagonal plots
            elif col == row:
                sns.distplot(plot_df.iloc[:, row], ax=ax, color="salmon", hist=False)

            #### Right side plot ->>> correlation value
            else:
                correlation = (
                    plot_df.iloc[:, [col, row]]
                    .corr(method=method)
                    .reset_index()
                    .iloc[1, 1]
                )
                ax.text(0.3, 0.5, "%.3f" % correlation, fontsize=20)

            if col != 0:
                ax.yaxis.set_visible(False)
            else:
                label = plot_df.columns[row]
                ax.set_ylabel(label)

            if row != box_size - 1:
                ax.xaxis.set_visible(False)
            else:
                label = plot_df.columns[col]
                ax.set_xlabel(label)
    sns.despine()


### COLOR ENCODER ##########


def assert_color_vector(categorical_vector, color_vector):
    categories = categorical_vector.unique()
    if len(categories) > len(color_vector):
        raise SeqUtilsError(
            "Not enough colors!! {} colors for {} categories".format(
                len(color_vector), len(categories)
            )
        )
    return categories.tolist()


class ColorEncoder:
    """
    color-encoding a categoric vector

    Example::

        categorical_vector = ['a','b','c','a']
        colors = obakeito_palette()
        ce = color_encoder()
        ce.fit(categorical_vector, colors)
        encoded_colors = ce.transform(new_categorical_vector)

    or::

        ce = color_encoder()
        encoded_colors = ce.fit_transform(categorical_vector, colors)

    access color encoder::

        encoded_color_map = ce.encoder

    """

    def __init__(self):
        self.x = None
        self.categories = None
        self.encoder = None  #: color encoder dictionary

    def fit(self, x, colors=okabeito_palette()):
        """
        Usage::

            ce = color_encoder()
            ce.fit(categroical_vector, colors)
        """
        self.x = Series(x)
        self.categories = assert_color_vector(self.x, colors)
        self.encoder = {c: col for c, col in zip(self.categories, colors)}

    def transform(self, xs):
        """
        Usage::

            ce = color_encoder()
            ce.fit(categroical_vector, colors)
            encoded_colors = ce.transform(new_categorical_vector)
        """
        if not self.encoder:
            raise ValueError("Call color_encoder.fit() first!!")

        union_set = set(self.categories).union(set(xs))
        if len(union_set) != len(self.categories):
            unseen = union_set - set(self.categories)
            unseen = ", ".join(list(unseen))
            raise SeqUtilsError("Contain unseen data!!: %s" % unseen)

        return Series(xs).map(self.encoder)

    def fit_transform(self, xs, colors=okabeito_palette()):
        """
        ce = color_encoder()
        encoded_colors = ce.fit_transform(categorical_vector, colors)
        """
        self.fit(xs, colors=colors)
        colors = self.transform(xs)
        return colors

    def show_legend(self, ax=None, sort=False, **kwargs):
        """
        Adding matplotlib legend

        """
        import matplotlib.patches as mpatches

        if sort:
            self.encoder = OrderedDict(sorted(self.encoder.items(), key=lambda x: x[0]))
        pat = [
            mpatches.Patch(color=col, label=lab) for lab, col in self.encoder.items()
        ]
        lgd = ax.legend(handles=pat, **kwargs)
        return lgd


color_encoder = ColorEncoder  # back compatibility


def mixed_sort(list_of_elements):
    """
    https://arcpy.wordpress.com/2012/05/11/sorting-alphanumeric-strings-in-python/

    Args:
       a list of

    Returns:
      list: ordered list

    Example::

        a = ['A1','A10','A2','A20']
        mixed_sort(a)
        > ['A1','A2', 'A10', 'A20']

    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", key)]
    return sorted(list_of_elements, key=alphanum_key)


def plot_upset(
    fig, upset_df, ylab="Intersected count", matrix_to_plot_ratio=0.4, fontsize=15
):
    """
    Args:
        fig: matplotlib figure object
        upset_df: pandas dataframe, contain column: sample matrix (binary), "count", 'index'
        ylab: y label
        matrix_to_plot_ratio: ratio between interaction matrix and batplot
        fontsize: fontsize

    input example::

        index   HeLa    K562    UHRR    Plasma  count
        3       0       1       1       1       1
        4       1       0       0       0       2
        5       1       1       0       0       3
        1       0       1       0       0       11
        6       1       1       0       1       12
        7       1       1       1       1       12
        2       0       1       0       1       13
        0       0       0       0       1       30


    """
    bar_ax = fig.add_axes([0, matrix_to_plot_ratio, 1, 1])
    heat_ax = fig.add_axes([0, 0, 1, matrix_to_plot_ratio])

    # plot upset bar
    upset_df.plot.bar("index", "count", width=0.9, ax=bar_ax)
    bar_ax.xaxis.set_visible(False)
    bar_ax.set_xlim(0, upset_df.shape[0] - 0.5)
    bar_ax.legend().set_visible(False)

    ymax = round(upset_df["count"].max(), -1)
    ticks = np.linspace(0, ymax, 5)
    ticks = np.round(ticks, -1)[1:]
    ticks = np.array(ticks, dtype="int")
    # bar_ax.set_yticks(ticks)
    # bar_ax.set_yticklabels(ticks, fontsize=fontsize)
    [bar_ax.spines[s].set_visible(False) for s in ["top", "right"]]
    bar_ax.set_ylabel(ylab, fontsize=fontsize)
    bar_ax.set_xlim(-0.5, upset_df.shape[0] - 0.5)

    matrix = upset_df.drop(["index", "count"], axis=1)
    sample_number = matrix.shape[1]
    heat_ax.imshow(matrix.transpose(), aspect="auto", cmap="binary", alpha=0.4)
    heat_ax.set_yticks(range(sample_number))
    heat_ax.set_yticklabels(matrix.columns, fontsize=fontsize)
    heat_ax.xaxis.set_visible(False)
    heat_ax.hlines(
        y=np.arange(sample_number) + 0.5, xmin=-0.5, xmax=upset_df.shape[0] - 0.5
    )
    heat_ax.vlines(
        x=np.arange(upset_df.shape[0]) + 0.5, ymax=sample_number + 0.5, ymin=-0.5
    )
    heat_ax.set_ylim(-0.5, sample_number - 0.5)

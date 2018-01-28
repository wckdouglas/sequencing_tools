import matplotlib.pyplot as plt
import seaborn as sns
from pandas import Series
from builtins import zip

def douglas_palette():
    '''
    Automatic set color if seaborn is installed, otherwise return list of colors
    usage: douglas_palette()

    ax=plt.subplot();
    for i in range(14):
        ax.plot(np.arange(10),np.arange(10) + i,label=i)
    ax.legend()
    '''
    colors = ['#B98476', '#6297B6', '#BF7F8E', '#8E955A', '#6F94B9', '#579F79', '#3C9F9D', 
            '#EDAEB5', '#B0BDEA', '#D8BA90', '#7BCBD5', '#87CDA9', '#B1C68D', '#E7C039']
    return colors

def okabeito_palette():
    '''
    Color palette proposed by Okabe and Ito
    copy from colorboindr R package 
    https://github.com/clauswilke/colorblindr/blob/master/R/palettes.R
    '''
    colors = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999"]
    return colors


def cor_plot(plot_df, fig=plt.figure(figsize=(17,17)), 
            diagonal_line = True, method = 'pearson'):
    '''
    given a data frame with columns storing data for each sample,
    output a matplotlib figure object as correlation plots.

    usage: cor_plot(plot_df, fig=plt.figure(figsize=(17,17)), diagonal_line = True, method = 'pearson')
    
    input:
    * plot_df: a pandas dataframe
    * fig: matplotlib figure object (optional)
    * diagonal_line: Boolean controlling if a diagonal line should be drawn
    * method: correlation method, [spearman or pearson]

    output:
    * fig: matplotlib figure object
    '''

    sns.set_style('white')
    assert method in ['spearman', 'pearson'], 'Wrong correlation method'
    box_size = len(plot_df.columns) 
    for row in range(box_size):
        for col in range(box_size):
            ax = fig.add_subplot(box_size, box_size, row * box_size + col + 1)

            #### Right side plots
            if col < row:
                ax.scatter(plot_df.iloc[:,col], plot_df.iloc[:,row])

                if diagonal_line:
                    plot_data = plot_df.iloc[:, [col,row]]
                    maxima = plot_data.max()
                    minima = plot_data.min()
                    ax.plot([minima,maxima],[minima,maxima], color='red')

            #### Diagonal plots
            elif col == row:
                sns.distplot(plot_df.iloc[:,row], 
                         ax = ax, color = 'salmon',
                        hist=False)

            #### Right side plot ->>> correlation value
            else:
                correlation = plot_df\
                    .iloc[:,[col,row]]\
                    .corr(method=method)\
                    .reset_index()\
                    .iloc[1,1]
                ax.text(0.3, 0.5, '%.3f' %correlation, fontsize=20)
            
            if col != 0:
                ax.yaxis.set_visible(False)
            else:
                label = plot_df.columns[row]
                ax.set_ylabel(label)
                
            if row != box_size -1:
                ax.xaxis.set_visible(False)
            else:
                label = plot_df.columns[col]
                ax.set_xlabel(label)        
    sns.despine()
    return fig



### COLOR ENCODER ##########

def assert_color_vector(categorical_vector, color_vector):
    categories = set(categorical_vector)
    assert len(categories) <= len(color_vector), 'Not enough colors!! %i colors for %i categories' %(len(color_vector),len(categories))
    return categories

class color_encoder():
    '''
    color-encoding a categoric vector

    Example:

    colors = obakeito_palette()
    ce = color_encoder()
    ce.fit(categroical_vector, colors)
    encoded_colors = ce.transform(new_categorical_vector) 

    #or

    ce = color_encoder()
    encoded_colors = ce.fit_transform(categorical_vector, colors)

    #access color encoder
    encoded_color_map = ce.encoder
    '''
    def __init__(self):
        self.x  = None
        self.categories = None
        self.encoder = None


    def fit(self, x, colors = okabeito_palette()):
        self.x = Series(x)
        self.categories = assert_color_vector(self.x, colors)
        self.encoder = {c:col for c, col in zip(self.categories, colors)}

    def transform(self, xs):
        if not self.encoder:
            raise ValueError('Call color_encoder.fit() first!!')

        if not self.categories.union(set(xs)) == self.categories:
            raise ValueError('Contain unseen data!!')

        return xs.map(encoder)

    def fit_transform(self, xs, colors = okabeito_palette()):
        self.x = Series(xs)
        self.categories = assert_color_vector(self.x, colors)
        self.encoder = {c:col for c, col in zip(self.categories, colors)}
        colors = self.x.map(self.encoder)
        return colors

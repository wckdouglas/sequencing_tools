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
    try:
        import seaborn as sns
        sns.set_style('white')
        palette = sns.color_palette(colors)
        sns.set_palette(palette)
        return colors
    except:
        return colors

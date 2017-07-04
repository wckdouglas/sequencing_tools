def douglas_palette():
    '''
    Automatic set color if seaborn is installed, otherwise return list of colors
    usage: douglas_palette()

    ax=plt.subplot();
    for i in range(17):
        ax.plot(np.arange(10),np.arange(10) + i,label=i)
    ax.legend()
    '''
    colors =['salmon',
         'khaki',
         'lightsteelblue',
         'plum',
         'darkseagreen',
         'darkorange',
         'palevioletred',
         'tan',
         'dimgray',
         'firebrick',
         'goldenrod',
         'navy',
         'mediumslateblue',
         'darkolivegreen',
        'darkgoldenrod',
        'crimson',
        'black']
    try:
        import seaborn as sns
        palette = sns.color_palette(colors)
        sns.set_palette(palette)
        return 0
    except:
        return colors

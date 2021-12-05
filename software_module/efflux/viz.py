# Import plotting utilities
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl

# Seaborn, useful for graphics
import seaborn as sns

import numpy as np


# Default RP plotting style
def pboc_style_mpl():
    """
    Formats matplotlib plotting enviroment to that used in 
    Physical Biology of the Cell, 2nd edition.
    """
    rc = {'lines.linewidth': 1.25,
          'axes.labelsize': 8,
          'axes.titlesize': 9,
          'axes.facecolor': '#E3DCD0',
          'xtick.labelsize': 7,
          'ytick.labelsize': 7,
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.color': '#ffffff',
          'legend.fontsize': 8,
          'figure.dpi': 300,
          'savefig.dpi': 300}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('xtick.major', pad=-1)
    plt.rc('ytick.major', pad=-1)
    plt.rc('mathtext', fontset='stixsans', sf='sansserif')
    plt.rc('figure', figsize=[3.5, 2.5])
    plt.rc('svg', fonttype='none')
    plt.rc('legend', title_fontsize='8', frameon=True, 
           facecolor='#E3DCD0', framealpha=1)
    sns.set_style('darkgrid', rc=rc)
    sns.set_palette("colorblind", color_codes=True)
    sns.set_context('notebook', rc=rc)
 
def pboc_style_bokeh():
    '''
    Formats bokeh plotting enviroment to that used in 
    Physical Biology of the Cell, 2nd edition.
    '''
    theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#E3DCD0',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'axis_line_color': "white",
            'major_tick_in': 7,
            'major_tick_line_width': 2.5,
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': 'white',
            },
            'Legend': {
                'background_fill_color': '#E3DCD0',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

    return theme_json

def pboc_color_selector():
    """
    Select the color palette of your choice.
    Parameters
    ----------
    style: str "mut" or "pboc"
        A string identifier for the style. "mut" gives colors for single and double mutants.
        "pboc" returns the PBoC2e color palette.
    Returns
    -------
    colors: dict
        Dictionary of colors. If "dna", "double", or "inducer" is the selected style,
        keys will be the mutants in upper case. Double mutant keys will be DNA-IND. For
        pboc, the keys will be the typical color descriptors. 
    """
    
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
             'pale_green': '#DCECCB', 'yellow': '#EAC264',
             'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
             'blue': '#738FC1', 'light_blue': '#A9BFE3',
             'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
             'pale_red': '#F1D4C9', 'purple': '#AB85AC',
             'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}
    return colors

def pboc_single(p):
    '''
    Formats a single bokeh plot to that used in 
    Physical Biology of the Cell, 2nd edition.
    Parameters
    ----------
    p : bokeh Figure
    '''
    # Figure
    p.background_fill_color = '#E3DCD0'
    p.outline_line_color = '#FFFFFF'

    # Axis
    p.axis.axis_line_color = "white"
    p.axis.major_tick_in = 7
    p.axis.major_tick_line_width = 2.5
    p.axis.major_tick_line_color = "white"
    p.axis.minor_tick_line_color = "white"
    p.axis.axis_label_text_font = 'Helvetica'
    p.axis.axis_label_text_font_style = 'normal'

    # Grid
    p.grid.grid_line_color = 'white'

    # Legend
    if p.legend:
        p.legend.background_fill_color = '#E3DCD0'
        p.legend.border_line_color = '#FFFFFF'
        p.legend.border_line_width = 1.5
        p.legend.background_fill_alpha = 0.5

    # Title
    if p.title:
        p.title.background_fill_color = '#FFEDC0'
        p.title.text_font_style = 'normal'
        p.title.align = 'center'
        p.title.text_font = 'Helvetica'
        p.title.offset = 2


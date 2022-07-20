YELLOW = '#b58900'
ORANGE = '#cb4b16'
RED = '#dc322f'
MAGENTA = '#d33682'
VIOLET = '#6c71c4'
BLUE = '#268bd2'
CYAN = '#2aa198'
GREEN = '#859900'
WHITE = '#ffffff'
BLACK = '#000000'
GREY = '#808080'
TURQ = '#0ffac9'

DEFAULT_SEQUENCE = [BLUE, RED, GREEN, YELLOW, VIOLET, CYAN, MAGENTA, ORANGE]

layout_white = dict(
    paper_bgcolor=WHITE,
    plot_bgcolor=WHITE,
    font_color=BLACK,
    xaxis_gridcolor=BLACK,
    yaxis_gridcolor=BLACK,

    # height=800, width=800,
    showlegend=False,
    font_family="Droid Sans",
)

scene3d_white = dict(
    xaxis_gridcolor=BLACK,
    xaxis_showbackground=False,
    yaxis_gridcolor=BLACK,
    yaxis_showbackground=False,
    zaxis_gridcolor=BLACK,
    zaxis_showbackground=False,
    aspectmode='cube',
)

layout3d_white = dict(
    scene=scene3d_white,
    margin=dict(l=0, r=0, b=0, t=0),
    legend=dict(
        orientation='h',
        itemsizing='constant',
        font_color=BLACK,
    ),
    **layout_white,
)

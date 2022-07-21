BASE03 = '#002b36'
BASE02 = '#073642'
BASE01 = '#586e75'
BASE00 = '#657b83'
BASE0 = '#839496'
BASE1 = '#93a1a1'
BASE2 = '#eee8d5'
BASE3 = '#fdf6e3'
YELLOW = '#b58900'
ORANGE = '#cb4b16'
RED = '#dc322f'
MAGENTA = '#d33682'
VIOLET = '#6c71c4'
BLUE = '#268bd2'
CYAN = '#2aa198'
GREEN = '#859900'

# discrete colors
DEFAULT_SEQUENCE = [BLUE, RED, GREEN, YELLOW, VIOLET, CYAN, MAGENTA, ORANGE]

layout_dark = dict(
  paper_bgcolor=BASE03,
  plot_bgcolor=BASE03,
  font_color=BASE1,
  xaxis_gridcolor=BASE0,
  yaxis_gridcolor=BASE0,
  font_family="Droid Sans",
)

scene3d_dark = dict(
    xaxis_gridcolor=BASE00,
    xaxis_showbackground=False,
    yaxis_gridcolor=BASE00,
    yaxis_showbackground=False,
    zaxis_gridcolor=BASE00,
    zaxis_showbackground=False,
    aspectmode='cube',
)

layout3d_dark = dict(
    # height=768, width=768,
    scene=scene3d_dark,
    margin=dict(l=0, r=0, b=0, t=0),
    legend=dict(
        orientation='h',
        itemsizing='constant',
        font_color=BASE2,
    ),
    **layout_dark,
)

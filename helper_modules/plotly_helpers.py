
import plotly.express as px
import plotly.graph_objects as go

def format_plotly(fig, height = 600, width = 800):
    fig.update_xaxes(ticks='outside', 
                     showline=True, 
                     linewidth=2., 
                     title_font=dict(size=22),
                     gridcolor='#E1D8C3',
                     linecolor='#273139', 
                     mirror = True)
    fig.update_yaxes(ticks='outside', 
                     showline=True, 
                     title_font=dict(size=22),
                     gridcolor='#E1D8C3',
                     linewidth=2., 
                     linecolor='#273139', 
                     mirror = True)
    fig.update_layout(
        height=height,
        width = width,
        paper_bgcolor = 'white',
        plot_bgcolor = 'white',
        template='plotly_white',
        hoverlabel=dict(
            bgcolor = 'white',
            font_size=14
        ),
        title={
        'y':0.9,
        'x':0.5,
        'xanchor': 'center',
        'yanchor': 'top'},
        xaxis = dict(
            title={
                #'text': 'First Dimension',
                'standoff': 0,
                'font': {
                    'size': 17
                } 
            },
            zeroline=True, 
            zerolinecolor='#BDB5A3', 
            zerolinewidth=2,
        ),
        yaxis = dict(
            title={
                #'text': 'Second Dimension',
                'standoff': 0,
                'font': {
                    'size': 17
                }
            },
            zeroline=True, 
            zerolinecolor='#BDB5A3', 
            zerolinewidth=2
        ),
        legend=dict(
            font=dict(
                color='#333333',
                size=12.5),
            orientation="v",
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.03,
            bgcolor="#eeeeee"
        ),
        dragmode='pan',
#         margin=dict(l=30, 
#                     r=30, 
#                     t=5, 
#                     b=30),
#         modebar=dict(orientation='h', 
#                      activecolor='#1f76b1')
    )
    return fig
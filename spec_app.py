#!/usr/bin/env python
import io
from base64 import b64encode

import dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd

# buffer = io.StringIO()
# df = px.data.iris()
# fig = px.scatter(
#     df, x="sepal_width", y="sepal_length", 
#     color="species")
# fig.write_html(buffer)

# html_bytes = buffer.getvalue().encode()
# encoded = b64encode(html_bytes).decode()

# app = dash.Dash(__name__)
# app.layout = html.Div([
#     dcc.Graph(id="graph", figure=fig),
#     html.A(
#         html.Button("Download HTML"), 
#         id="download",
#         href="data:text/html;base64," + encoded,
#         download="plotly_graph.html"
#     )
# ])

app = dash.Dash(__name__)

# assume you have a "long-form" data frame
# see https://plotly.com/python/px-arguments/ for more options
df = pd.DataFrame({
    "Fruit": ["Testing", "Oranges", "Bananas", "Apples", "Oranges", "Bananas"],
    "Amount": [4, 1, 2, 2, 4, 5],
    "City": ["SF", "SF", "SF", "Montreal", "Montreal", "Montreal"]
})

fig = px.bar(df, x="Fruit", y="Amount", color="City", barmode="group")

app.layout = html.Div(children=[
    html.H1(children='Hello Dash'),

    html.Div(children='''
        Dash: A web application framework for your data.
    '''),

    dcc.Graph(
        id='example-graph',
        figure=fig
    )
])

if __name__ == '__main__':
    app.run_server(debug=True)





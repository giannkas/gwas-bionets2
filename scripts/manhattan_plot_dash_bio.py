from dash import Dash, dcc, html, Input, Output
import dash_bio as dashbio
import pandas as pd


app = Dash(__name__)

df = pd.read_csv('results/psoriasis/magma/magma_output_50kb_window_snp-array/pso.scores.genes.out_converted.tsv', sep="\t")

figure = dashbio.ManhattanPlot(
	dataframe=df,
	chrm='Chr',
	bp='Start',
	p='Pvalue',
	snp='GENE',
	gene='Gene',
	title='Manhattan plot for MAGMA gene P-values',
	suggestiveline_value=False,
	showlegend=False,
	showgrid=False,
)
figure.write_image("manhattan_plot.pdf", width=1200, height=600)

app.layout = html.Div([
    'Threshold value',
    dcc.Slider(
        id='slider',
        min=1,
        max=10,
        marks={
            i: {'label': str(i)} for i in range(10)
        },
        value=6
    ),
    html.Br(),
    html.Div(
        dcc.Graph(
            id='graph',
            figure=dashbio.ManhattanPlot(
                dataframe=df,
								chrm='Chr',
								bp='Start',
								p='Pvalue',
								snp='GENE',
								gene='Gene',
								title='Manhattan plot for MAGMA gene P-values'
            )
        )
    )
])

@app.callback(
    Output('graph', 'figure'),
    Input('slider', 'value')
)
def update_manhattanplot(threshold):

    return dashbio.ManhattanPlot(
        dataframe=df,
        chrm='Chr',
				bp='Start',
				p='Pvalue',
				snp='GENE',
				gene='Gene',
				title='Manhattan plot for MAGMA gene P-values',
        genomewideline_value=threshold
    )


if __name__ == '__main__':
    app.run(debug=True)
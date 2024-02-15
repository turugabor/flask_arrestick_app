import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go


class Plotter:
    def __init__(self):
        pass

    def plot(self, seq, probability, confidence, name):
        """
        Plot the arreSTick sequence probability.

        Args:
            seq (str): The amino acid sequence.
            probability (list): The probability values.
            confidence (list): The confidence values.

        Returns:
            plotly.graph_objects.Figure: The plotly figure.
        """
        data = self._prepare_data(seq, probability, confidence)
        return self._plot(data, name)

    def _prepare_data(self, seq, probability, confidence):
        """
        Prepare the data for plotting.

        Args:
            seq (str): The amino acid sequence.
            probability (list): The probability values.
            confidence (list): The confidence values.

        Returns:
            pandas.DataFrame: The prepared data.
        """
        data = pd.DataFrame([probability, list(seq), np.arange(len(seq)) + 1],
                            index=["Probability", "Amino acid", "Amino acid position"]).transpose()
        data["region"] = ""
        
        for i, row in data.iterrows():
            data.loc[i, "region"] = f"{i+1}-{i+16} {seq[i:i+15]}"
        
        data['alphafold'] = [1 if conf > 70 else 0 for conf in confidence]
        data["structure"] = ['>70' if conf > 70 else '≤70' for conf in confidence]
        
        return data

    def _plot(self, data, name):
        """
        Plot the arreSTick sequence probability.

        Args:
            data (pandas.DataFrame): The prepared data.

        Returns:
            plotly.graph_objects.Figure: The plotly figure.
        """
        ymin = 0
        ymax = 1
        
        subfig = make_subplots(specs=[[{"secondary_y": True}]])
        
        plot = px.line(data_frame=data, y="Probability", x="Amino acid position", hover_data=["region"],
                       color_discrete_sequence=["black"])
        
        plot.update_layout(yaxis_range=[ymin, 1])
        plot.update_traces(customdata=data.region, hovertemplate='Position: %{x} <br>Probability: %{y:.2f} <br>Region: %{customdata} ')
        
        heatmap_range = np.arange(ymin, ymax)
        
        trace1 = go.Heatmap(
            z=np.concatenate([data.alphafold.values.reshape(1, -1)]),
            x=data["Amino acid position"],
            y=[0, 1],
            colorscale=["#ef8354", "#2d3142"],
            customdata=np.concatenate([data.structure.values.reshape(1, -1)] * heatmap_range.shape[0]),
            hovertemplate='Alphafold confidence: %{customdata}',
            name="",
            showscale=False,
            opacity=0.6
        )

        
        
        
        subfig.add_trace(trace1)
        subfig.add_traces(plot.data)
        
        subfig.layout.xaxis.title = "Amino acid position"
        subfig.layout.yaxis.title = "arreSTick sequence probability"
        subfig.layout.yaxis2.update(visible=False, showticklabels=False)
        
        subfig.add_shape(
            type="rect",
            x0=0, y0=1.05, x1=20, y1=1.10,
            line=dict(color="#ef8354"), fillcolor="#ef8354", opacity=0.6
        )
        
        subfig.add_annotation(
            x=25, y=1.08, xanchor="left",
            text="Alphafold confidence ≤ 70",
            showarrow=False
        )
        
        subfig.add_shape(
            type="rect",
            x0=0, y0=1.12, x1=20, y1=1.17,
            line=dict(color="#2d3142"), fillcolor="#2d3142", opacity=0.6
        )
        
        subfig.add_annotation(
            x=25, y=1.15, xanchor="left",
            text="Alphafold confidence > 70",
            showarrow=False
        )
        
        subfig.update_layout(paper_bgcolor='white', plot_bgcolor='white', title=f"Protein sequence prediction of {name}")
         
        return subfig

from flask import Flask
from flask import render_template, session, url_for, redirect, request
from flask_restful import Resource, Api
import plotly
from sequence import Convolution
from plot import Plotter
import json
import os

from logging.config import dictConfig


dictConfig({
    'version': 1,
    'formatters': {'default': {
        'format': '[%(asctime)s] %(levelname)s in %(module)s: %(message)s',
    }},
    'handlers': {'wsgi': {
        'class': 'logging.StreamHandler',
        'stream': 'ext://flask.logging.wsgi_errors_stream',
        'formatter': 'default'
    }},
    'root': {
        'level': 'INFO',
        'handlers': ['wsgi']
    }
})

app = Flask(__name__)
app.secret_key = os.urandom(16)
api = Api(app)

conv = Convolution()
plotter = Plotter()

app.logger.info(f"session")

@app.route('/')
@app.route('/home')
def home():
    session.permanent = False

    if 'protein' in session:
        protein = session['protein']
        app.logger.info(f"protein: {protein}")

        graphJSON = create_graph(protein)
        message = ''

        if graphJSON != None:
            return render_template("index.html",
                                graphJSON=graphJSON,
                                message=message, 
                                title="arreSTick", )
        else:
            message = "Protein not found or sequence is not valid"
            return render_template("index.html", message=message, title="arreSTick")
        
    else:
        return render_template("index.html", message='', title="arreSTick")
    
@app.route('/select', methods=["GET", "POST"])
def select():
    app.logger.info(f"select")
    if request.method == 'POST':
        protein = request.form.get('protein')

        session['protein'] = protein

        return redirect(url_for('home'), code=302)
    
def create_graph(protein):
    protein = protein.upper().replace(" ", "").replace("\n", "").replace("\r", "")
    sequence, probability, confidence = conv.get_arrestick_data(protein)

    if sequence == None:
        session["message"] = "Protein not found"
        return None

    else:
        fig = plotter.plot(sequence, probability, confidence)  
        graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
        return graphJSON
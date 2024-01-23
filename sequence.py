import numpy as np
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import urllib.request
import biotite.structure.io.pdbx as pdbx
import io
import yaml
from yaml import CLoader as Loader


class Convolution:
    def __init__(self):
        self.entries = self.get_dicts()
        self.aa_dict, self.kernel, self.bias, self.sigmoid_weight, self.sigmoid_bias = self.get_params()

    def get_dicts(self):
        """Load data for protein label dictionaries (entry - entry name conversions)"""
        uniprot = pd.read_csv('static/uniprot_entries.tsv', sep='\t')
        entries = uniprot.set_index("Entry Name")["Entry"].to_dict()
        return entries

    def get_params(self):
        """Load model parameters from a YAML file"""
        with open("static/model_params", "r") as f:
            params = yaml.load(f, Loader)
            
        aa_dict = params["model_1"]["aa_dict"]
        kernel = np.array(params["model_1"]["kernel_weights"])
        bias = params["model_1"]["conv_bias"]
        sigmoid_weight = params["model_1"]["sigmoid_weight"]
        sigmoid_bias = params["model_1"]["sigmoid_bias"]
        
        return aa_dict, kernel, bias, sigmoid_weight, sigmoid_bias

    def sigmoid(self, x):
        """Apply the sigmoid function to the input x"""
        if x is not None:
            return 1 / (1 + np.exp(-self.sigmoid_weight * x - self.sigmoid_bias))
        else:
            return None

    def convolve(self, seq):
        """Perform convolution on the input sequence using the given kernel, amino acid dictionary, and bias"""
        seq_translated = np.array(list(map(lambda x: self.aa_dict[x], seq)))
        convoluted = np.convolve(seq_translated, self.kernel[::-1], mode='valid') + self.bias
        return convoluted
    
    def get_probability(self, convoluted):
        return self.sigmoid(convoluted)

    def get_alphafold_data(self, entry_name):
        """Download alphafold .cif file for the input protein and extract sequence and model confidence array

        Args:
            entry_name (str): Uniprot entry name of the protein

        Returns:
            tuple of sequence and confidence array: The protein sequence from the alphafold website with the model confidence array. 
            If the entry is not available on the server, "x" instead of the sequence and an array of [0] as the confidence is returned
        """
        if "_" in entry_name:
            entry = self.entries[entry_name]
        else:
            entry = entry_name

        connection = urllib.request.urlopen(f"https://alphafold.ebi.ac.uk/files/AF-{entry}-F1-model_v4.cif")
        databytes = connection.read()
        connection.close()
        cif_txt = databytes.decode("utf8")

        f = io.StringIO(cif_txt)
        cif = pdbx.PDBxFile.read(f)
        
        confidence = pd.DataFrame(cif.get_category("ma_qa_metric_local")).metric_value.astype(float).values
        sequence = cif.get_category("entity_poly")["pdbx_seq_one_letter_code"]

        return sequence, confidence
    
    def get_arrestick_data(self, protein):
        """Get the arreSTick data for the input protein sequence or entry name"""

        try:

            if protein in self.entries.values():
                entry = protein
                sequence, confidence = self.get_alphafold_data(entry)
            elif protein in self.entries.keys():
                sequence, confidence = self.get_alphafold_data(self.entries[protein])
            else:
                # check if protein is a valid amino acid sequence
                if all(aa in self.aa_dict.keys() for aa in protein):
                    sequence = protein
                    confidence = np.zeros(len(sequence))
                else:
                    sequence = None
            
            if sequence is None:
                return None, None, None
            else:
                convoluted = self.convolve(sequence)
                probability = self.get_probability(convoluted)        
                return sequence, probability, confidence
            
        except ValueError:
            return None, None, None
        
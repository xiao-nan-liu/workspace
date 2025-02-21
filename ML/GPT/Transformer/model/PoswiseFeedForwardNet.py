import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn

class PoswiseFeedForwardNet(nn.Module):
    def __init__(self, d_model, d_ff):
        super(PoswiseFeedForwardNet, self).__init__()
        self.fc = nn.Sequential(
        nn.Linear(d_model, d_ff, bias = False),
        nn.ReLU(),
        nn.Linear(d_ff, d_model, bias = False)
        )
        self.layernorm = nn.LayerNorm(d_model)
    
    def forward(self, inputs):
        ###inputs [batch_size, seq_len, d_model]
        residual = inputs
        output = self.fc(inputs)
        return self.layernorm(output+residual)
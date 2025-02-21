import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn
class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_pos, device):
        super(PositionalEncoding, self).__init__()
        self.device = device
        self.pos_embedding = nn.Embedding(max_pos, d_model)

    def forward(self, inputs):
        seq_len = inputs.size(1)
        pos = torch.arange(seq_len, dtype=torch.long, device=self.device)
        # [seq_len] -> [batch_size, seq_len]
        pos = pos.unsqueeze(0).expand_as(inputs)
        return self.pos_embedding(pos)
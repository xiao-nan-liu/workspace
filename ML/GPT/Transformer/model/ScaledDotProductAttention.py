import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn
class ScaledDotProductAttention(nn.Module):
    def __init__(self, d_k):
        super(ScaledDotProductAttention, self).__init__()
        self.d_k = d_k
    def forward(self, q, k, v, attention_mask):
        ####Tensor dimension####################
        #### q [batch_size, n_heads, len_q, d_k]
        #### k [batch_size, n_heads, len_k, d_k]
        #### v [batch_size, n_heads, len_v, d_v]
        ####attension_mask [batch_size, n_heads, seq_len, seq_len]
        ####
        ####calculate score q and k
        score = torch.matmul(q, k.transpose(-2, -1))/np.sqrt(self.d_k)
        ##original code 
        ####attn = torch.matmul(q / self.temperature, k.transpose(2, 3))
        ####extend the length and fill a small number
        score.masked_fill_(attention_mask, -1e10)
        ####softmax with last dimension
        attn = nn.Softmax(dim=-1)(score)     
        ####After attention [batch_size, n_heads, len_q, d_v]
        context = torch.matmul(attn, v)
        return context, attn
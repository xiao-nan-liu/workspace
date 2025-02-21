import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn
import ScaledDotProductAttention
class MultiHeadAttention(nn.Module):
    def __init__(self, d_model, n_heads, d_k, d_v):
        super(MultiHeadAttention, self).__init__()
        self.d_model = d_model
        self.n_heads = n_heads
        self.d_k = d_k
        self.d_v = d_v
        self.w_q = nn.Linear(d_model, d_k * n_heads, bias = False)
        self.w_k = nn.Linear(d_model, d_k * n_heads, bias = False)
        self.w_v = nn.Linear(d_model, d_v * n_heads, bias = False)
        self.fc = nn.Linear(n_heads * d_v, d_model, bias=False)
        self.layernorm = nn.LayerNorm(d_model)
        
    def forward(self, q, k, v, attention_mask):
        ####Tensor dimension####################
        #### q [batch_size, n_heads, len_q, d_k]
        #### k [batch_size, n_heads, len_k, d_k]
        #### v [batch_size, n_heads, len_v, d_v]
        ####attension_mask [batch_size, n_heads, seq_len, seq_len]
        residual, batch_size = q, q.size(0)
        ####calculate q, k, v####
        q = self.w_q(q).view(batch_size, -1, self.n_heads, self.d_k).transpose(1,2)
        k = self.w_k(k).view(batch_size, -1, self.n_heads, self.d_k).transpose(1,2)
        v = self.w_v(v).view(batch_size, -1, self.n_heads, self.d_v).transpose(1,2)
        ####calculate attention mask####
        attention_mask = attention_mask.unsqueeze(1).repeat(1, self.n_heads, 1, 1)
        ####calculate score####
        context, attn = ScaledDotProductAttention(self.d_k)(q, k, v, attention_mask)
        # context [batch_size, len_q, n_heads * d_v]
        context = context.transpose(1,2).reshape(batch_size, -1, self.n_heads*self.d_v)
        output = self.fc(context)
        return self.layernorm(output+residual), attn
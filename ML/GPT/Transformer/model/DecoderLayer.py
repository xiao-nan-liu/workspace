import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn
import MultiHeadAttention
import PoswiseFeedForwardNet
class DecoderLayer(nn.Module):
    def __init__(self, d_model, n_heads, d_ff, d_k, d_v):
        super(DecoderLayer, self).__init__()
        # Multi-head attention
        self.attention = MultiHeadAttention(d_model, n_heads, d_k, d_v)
        # feedforward layer
        self.pos_ffn = PoswiseFeedForwardNet(d_model, d_ff)
        
    def forward(self, inputs, attention_mask):
        ##
        # inputs: [batch_size, seq_len, d_model]
        # attention_mask: [batch_size, seq_len, seq_len]
        ##
        # outputs: [batch_size, seq_len, d_model]
        # self_attn: [batch_size, n_heads, seq_len, seq_len]
        outputs, self_attn = self.attention(inputs, inputs, inputs, attention_mask)
        # [batch_size, seq_len, d_model]
        outputs = self.pos_ffn(outputs)
        return outputs, self_attn
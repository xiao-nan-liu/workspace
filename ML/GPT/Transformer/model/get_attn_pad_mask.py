import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn

def get_attn_pad_mask(attention_mask):
    batch_size, len_seq = attention_mask.size()
    attention_mask = attention_mask.data.eq(0).unsqueeze(1)
    return attention_mask.expand(batch_size, len_seq, len_seq)
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import torch
from torch import nn

def get_attn_subsequence_mask(seq, device):
    # attention score dimention [batch_size, n_heads, len_seq, len_seq]
    # mask generation [batch_size, len_seq, len_seq]
    attn_shape = [seq.size(0), seq.size(1), seq.size(1)]
    # generate an upper left matrix
    subsequence_mask = np.triu(np.ones(attn_shape), k=1)
    subsequence_mask = torch.from_numpy(subsequence_mask).byte()
    subsequence_mask = subsequence_mask.to(device)
    return subsequence_mask
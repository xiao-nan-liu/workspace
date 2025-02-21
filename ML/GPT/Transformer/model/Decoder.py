class Decoder(nn.Module):
    def __init__(self, d_model, n_heads, d_ff, d_k, d_v, vocab_size, max_pos, n_layers, device):
        super(Decoder, self).__init__()
        self.device = device
        # Tokenization
        self.embedding = nn.Embedding(vocab_size, d_model)
        # positional encodinng
        self.pos_encoding = PositionalEncoding(d_model, max_pos, device)
        self.layers = nn.ModuleList([DecoderLayer(d_model, n_heads, d_ff, d_k, d_v) for _ in range(n_layers)])

    def forward(self, inputs, attention_mask):
        ##
        # inputs: [batch_size, seq_len]
        ##
        # [batch_size, seq_len, d_model]
        outputs = self.embedding(inputs) + self.pos_encoding(inputs)
        # upper triangle，prevent seeing future info， [batch_size, seq_len, seq_len]
        subsequence_mask = get_attn_subsequence_mask(inputs, self.device)
        if attention_mask is not None:
            # padding [batch_size, seq_len, seq_len]
            attention_mask = get_attn_pad_mask(attention_mask)
            # [batch_size, seq_len, seq_len]
            attention_mask = torch.gt((attention_mask + subsequence_mask), 0)
        else:
            attention_mask = subsequence_mask.bool()
        # calculate each layer
        self_attns = []
        for layer in self.layers:
            # outputs: [batch_size, seq_len, d_model],
            # self_attn: [batch_size, n_heads, seq_len, seq_len],
            outputs, self_attn = layer(outputs, attention_mask)
            self_attns.append(self_attn)
        return outputs, self_attns
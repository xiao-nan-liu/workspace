###predict####

import torch

from model import GPTModel
from tokenizer import Tokenizer


def generate(model, tokenizer, text, max_length, device):
    input, att_mask = tokenizer.encode(text)
    input = torch.tensor(input, dtype=torch.long, device=device).unsqueeze(0)
    stop = False
    input_len = len(input[0])
    while not stop:
        if len(input[0]) - input_len > max_length:
            next_symbol = tokenizer.sep_token
            input = torch.cat(
                [input.detach(), torch.tensor([[next_symbol]], dtype=input.dtype, device=device)], -1)
            break
        projected, self_attns = model(input)
        prob = projected.squeeze(0).max(dim=-1, keepdim=False)[1]
        next_word = prob.data[-1]
        next_symbol = next_word
        if next_symbol == tokenizer.sep_token:
            stop = True
        input = torch.cat(
            [input.detach(), torch.tensor([[next_symbol]], dtype=input.dtype, device=device)], -1)
    decode = tokenizer.decode(input[0].tolist())
    decode = decode[len(text):]
    return "".join(decode)


def main():
    model_path = "/kaggle/working/best.pt"
    vocab_path = "/kaggle/input/vocabulary/vocab.json"  # vocabulary
    max_length = 128  # maximum length
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # tokenizer loading
    tokenizer = Tokenizer(vocab_path)
    # model parameters
    model_param = {
        "d_model": 768,  # embedding size
        "d_ff": 2048,  # ffn size
        "d_k": 64,  # K
        "d_v": 64,  # V
        "n_layers": 6,  # number of decoder layers
        "n_heads": 8,  # number of attention heads
        "max_pos": 1800,  # length of positional encoding
        "device": device,  # device
        "vocab_size": tokenizer.get_vocab_size(),  # vocabulary size
    }
    model = GPTModel(**model_param)
    model.load_state_dict(torch.load(model_path))
    model.to(device)

    while True:
        text = input("请输入：")
        if not text:
            continue
        if text == "q":
            break
        res = generate(model, tokenizer, text, max_length, device)
        print("AI: ", res)


if __name__ == '__main__':
    main()
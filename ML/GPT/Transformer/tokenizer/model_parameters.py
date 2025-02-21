import torch

from model import GPTModel


def main():
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # 模型参数
    model_param = {
        "d_model": 768,  # embedding
        "d_ff": 2048,  # forward feed layer
        "d_k": 64,  # K 
        "d_v": 64,  # V 
        "n_layers": 6,  # decoding layer
        "n_heads": 8,  # multi-head attenntion
        "max_pos": 1800,  # positional encoding length
        "device": device,  # device
        "vocab_size": 4825  # vocabulary size
    }
    model = GPTModel(**model_param)
    total_params = sum(p.numel() for p in model.parameters())
    print(model)
    print("total_params: ", total_params)


if __name__ == '__main__':
    main()
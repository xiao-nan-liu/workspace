import torch
from torch.utils.data import DataLoader
from torch.utils.tensorboard import SummaryWriter
from tokenizer import Tokenizer
from model import GPTModel
import QADataset
from tqdm import tqdm
import time, sys, os


def train_model(model, train_loader, val_loader, optimizer, criterion,
                device, num_epochs, model_output_dir, writer):
    batch_step = 0
    best_val_loss = float('inf')
    for epoch in range(num_epochs):
        time1 = time.time()
        model.train()
        for index, data in enumerate(tqdm(train_loader, file=sys.stdout, desc="Train Epoch: " + str(epoch))):
            input_ids = data['input_ids'].to(device, dtype=torch.long)
            attention_mask = data['attention_mask'].to(device, dtype=torch.long)
            labels = data['labels'].to(device, dtype=torch.long)
            optimizer.zero_grad()
            outputs, dec_self_attns = model(input_ids, attention_mask)
            loss = criterion(outputs, labels.view(-1))
            loss.backward()
            # gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1)
            optimizer.step()
            writer.add_scalar('Loss/train', loss, batch_step)
            batch_step += 1
            # print loss every 100 times
            if index % 100 == 0 or index == len(train_loader) - 1:
                time2 = time.time()
                tqdm.write(
                    f"{index}, epoch: {epoch} -loss: {str(loss)} ; lr: {optimizer.param_groups[0]['lr']} ;each step's time spent: {(str(float(time2 - time1) / float(index + 0.0001)))}")
        # validation
        model.eval()
        val_loss = validate_model(model, criterion, device, val_loader)
        writer.add_scalar('Loss/val', val_loss, epoch)
        print(f"val loss: {val_loss} , epoch: {epoch}")
        # save the best model
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_path = os.path.join(model_output_dir, "best.pt")
            print("Save Best Model To ", best_model_path, ", epoch: ", epoch)
            torch.save(model.state_dict(), best_model_path)
        # save current model
        last_model_path = os.path.join(model_output_dir, "last.pt")
        print("Save Last Model To ", last_model_path, ", epoch: ", epoch)
        torch.save(model.state_dict(), last_model_path)


def validate_model(model, criterion, device, val_loader):
    running_loss = 0.0
    with torch.no_grad():
        for _, data in enumerate(tqdm(val_loader, file=sys.stdout, desc="Validation Data")):
            input_ids = data['input_ids'].to(device, dtype=torch.long)
            attention_mask = data['attention_mask'].to(device, dtype=torch.long)
            labels = data['labels'].to(device, dtype=torch.long)
            outputs, dec_self_attns = model(input_ids, attention_mask)
            loss = criterion(outputs, labels.view(-1))
            running_loss += loss.item()
    return running_loss / len(val_loader)


def main():
    train_json_path = "train.json"  # Training set
    val_json_path = "val.json"  # Validation set
    vocab_path = "/kaggle/input/vocabulary/vocab.json"  # Vocabulary path
    max_length = 120  # maximum length
    epochs = 15 # epoch
    batch_size = 128  # each batch 128 
    lr = 1e-4  # leaerning rate
    model_output_dir = "/kaggle/working/"  # save model to dir
    logs_dir = "logs"  # log dir
    # device
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    # Tokenizer loading
    tokenizer = Tokenizer(vocab_path)
    # model parameter
    model_param = {
        "d_model": 768,  # embedding
        "d_ff": 2048,  # forward feed net
        "d_k": 64,  # K
        "d_v": 64,  # V
        "n_layers": 6,  # number of decoding layers
        "n_heads": 8,  # number of multi-head attention heads
        "max_pos": 1800,  # length of positional encoding
        "device": device,  # device
        "vocab_size": tokenizer.get_vocab_size(),  # vocabulary size
    }
    model = GPTModel(**model_param)
    print("Start Load Train Data...")
    train_params = {
        "batch_size": batch_size,
        "shuffle": True,
        "num_workers": 4,
    }
    training_set = QADataset(train_json_path, tokenizer, max_length)
    training_loader = DataLoader(training_set, **train_params)
    print("Start Load Validation Data...")
    val_params = {
        "batch_size": batch_size,
        "shuffle": False,
        "num_workers": 4,
    }
    val_set = QADataset(val_json_path, tokenizer, max_length)
    val_loader = DataLoader(val_set, **val_params)
    # log 
    writer = SummaryWriter(logs_dir)
    # optimizer
    optimizer = torch.optim.AdamW(params=model.parameters(), lr=lr)
    # loss function
    criterion = torch.nn.CrossEntropyLoss(ignore_index=0).to(device)
    model = model.to(device)
    # Start training
    print("Start Training...")
    train_model(
        model=model,
        train_loader=training_loader,
        val_loader=val_loader,
        optimizer=optimizer,
        criterion=criterion,
        device=device,
        num_epochs=epochs,
        model_output_dir=model_output_dir,
        writer=writer
    )
    writer.close()


if __name__ == '__main__':
    main()
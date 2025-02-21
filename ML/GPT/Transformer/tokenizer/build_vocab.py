import json

def build_vocab(file_path):
    # Read all texts
    texts = []
    with open(file_path, 'r', encoding='utf-8') as r:
        for line in r:
            if not line:
                continue
            line = json.loads(line)
            question = line["question"]
            answer = line["answer"]
            texts.append(question)
            texts.append(answer)
    # Split Token
    words = set()
    for t in texts:
        if not t:
            continue
        for word in t.strip():
            words.add(word)
    words = list(words)
    words.sort()
    # Special Token
    # pad、unk UNKNOWN、sep END
    word2id = {"<pad>": 0, "<unk>": 1, "<sep>": 2}
    # build up vocabulary list
    word2id.update({word: i + len(word2id) for i, word in enumerate(words)})
    id2word = list(word2id.keys())
    vocab = {"word2id": word2id, "id2word": id2word}
    vocab = json.dumps(vocab, ensure_ascii=False)
    with open('/kaggle/input/vocabulary/vocab.json', 'r', encoding='utf-8') as w:
        temp = w.read()
        vocab += temp
        #w.write(vocab)
    print(f"finish. words: {len(id2word)}")

if __name__ == '__main__':
    build_vocab("/kaggle/input/chinese-gpt-data/train.jsonl")
with open('meth_out.txt', encoding='gbk', errors='ignore') as f:
    text = f.read()
with open('meth_clean.txt', 'w', encoding='utf-8') as f:
    f.write(text)

# 关于Biopython的笔记

```python
def organism_name(seq_record):
    'get organism name from seq_record and remove space in scientific name'
    name_string = seq_record.annotations["organism"]
    new_name = "_".join(name_string.split())
    return new_name
```
其中，*seq_record*是SeqIO解析出来的Seq对象，通过annotation（字典）获取其中包含的信息，一般为除feature表外的其它信息，如：organism， keywords等

---
```python
for seq_record in SeqIO.parse(args.genbank, 'genbank'):
    #print(f'Dealing with GenBank record \n>{seq_record.id}')
    for seq_feature in seq_record.features:
        if seq_feature.type == args.feature:
            #with open(args.outfile, 'w') as fo:
            print(f'>{seq_record.name}|{organism_name(seq_record)}{keywords(seq_record)} {seq_feature.qualifiers["note"]}\
            \n{subseq(seq_record, seq_feature)}')
```
其中，由seq_record.features获取所有的features表，得到一个可迭代对象，（iterator？），使用for循环，取单一feature记录中的具体信息，
seq_feature.type提取的是例如“CDS，tRNA，rRNA，misc_feature”等特征，本例中，使用seq_feature.type判断是否是指定的特征，从而打印当前feature的信息，
包括注释信息（note，product，gene等标签信息）

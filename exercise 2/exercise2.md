```python
# pip install pysam
import pysam
```


```python
flags = ["    read paired (0x1)",
"    read mapped in proper pair (0x2)",
"    read unmapped (0x4)",
"    mate unmapped (0x8)",
"    read reverse strand (0x10)",
"    mate reverse strand (0x20)",
"    first in pair (0x40)",
"    second in pair (0x80)",
"    not primary alignment (0x100)",
"    read fails platform/vendor quality checks (0x200)",
"    read is PCR or optical duplicate (0x400)",
"    supplementary alignment (0x800)"]
```


```python
# first read
tmpfilename = "/sbgenomics/project-files/merged-tumor.bam"
infile = pysam.AlignmentFile(tmpfilename, "rb")
for read in infile:
    print(read)
    break
```

    C0HVYACXX120402:7:1207:5722:57044	1187	20	9483248	27	76M	20	9483381	76	TTTTCAAACAGTATCTATGCCTGCCAAATGTGAACATATAAAAAAAAACCAGAATGTGCCATTCTGATTTAAACTG	array('B', [28, 28, 27, 29, 31, 30, 31, 31, 29, 31, 35, 30, 29, 31, 34, 30, 29, 23, 41, 32, 20, 30, 29, 34, 34, 29, 30, 31, 30, 30, 30, 33, 33, 26, 39, 12, 25, 19, 32, 30, 35, 28, 35, 33, 23, 33, 35, 36, 30, 38, 33, 41, 34, 35, 31, 33, 23, 30, 30, 36, 27, 32, 29, 34, 35, 41, 33, 31, 33, 29, 32, 32, 31, 31, 31, 34])	[('XA', 'GL000217.1,-110754,76M,1;'), ('BD', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'), ('MD', '76'), ('RG', '1'), ('BI', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'), ('NM', 0), ('MQ', 27), ('AS', 76), ('XS', 71)]



```python
# flag field 
infile = pysam.AlignmentFile(tmpfilename, "rb")
flag = 0
for read in infile:
    flag = read.flag
    break
    
print(f"Flag field in the first read is: {read.flag}")

# print flag
index = 0
while flag > 0:
    if (flag & 1 == 1):
        print(flags[index])
    index += 1
    flag >>= 1
```

    Flag field in the first read is: 1187
        read paired (0x1)
        read mapped in proper pair (0x2)
        mate reverse strand (0x20)
        second in pair (0x80)
        read is PCR or optical duplicate (0x400)



```python
# count unmapped
infile = pysam.AlignmentFile(tmpfilename, "rb")
cnt = 0
for read in infile:
    cnt += (read.flag & 0x4) >> 2
print(f"Count of unmapped: {cnt}")
```

    Count of unmapped: 17765



```python
# count reads
infile = pysam.AlignmentFile(tmpfilename, "rb")
cnt = 0
for read in infile:
    cnt += 1
print(f"Total number of reads: {cnt}")
```

    Total number of reads: 2921629



```python
# mapping quality 0
infile = pysam.AlignmentFile(tmpfilename, "rb")
cnt = 0
for read in infile:
    cnt += 1 if read.mapping_quality == 0 else 0 
print(f"Number of reads with mapping quality 0: {cnt}")
```

    Number of reads with mapping quality 0: 126628



```python
# avg mapping quality
infile = pysam.AlignmentFile(tmpfilename, "rb")
cnt = 0
avg = 0
for read in infile:
    cnt += 1
    avg += read.mapping_quality

print(f"Average mapping quality: {round(avg/cnt, 3)}")
```

    Average mapping quality: 55.914



```python
# avg mapping quality over 0
infile = pysam.AlignmentFile(tmpfilename, "rb")
cnt = 0
avg = 0
for read in infile:
    if (read.mapping_quality != 0):
        cnt += 1
        avg += read.mapping_quality

print(f"Average mapping quality : {round(avg/cnt, 3)}")
```

    Average mapping quality : 58.447



```python
# now to do it all at once 
infile = pysam.AlignmentFile(tmpfilename, "rb")
size = 0
unmapped_count = 0
unmapped_zero_count = 0
avg_all = 0
avg_filter = 0
for read in infile:
    size += 1
    unmapped_count += (read.flag & 0x4) >> 2
    unmapped_zero_count += 1 if (read.mapping_quality == 0) else 0
    avg_all += read.mapping_quality 
    avg_filter += read.mapping_quality if (read.mapping_quality > 0) else 0

infile.close()
print(f"Number of reads: {size}")
print(f"Number of unmapped reads: {unmapped_count}")
print(f"Number of unmapping reads with quality 0: {unmapped_zero_count}")
print(f"Average mapping quality for all the reads: {round(avg_all/size, 3)}")
print(f"Average mapping quality if reads with 0 map quality are filtered out: {round(avg_filter/(size - unmapped_zero_count), 3)}")
```

    Number of reads: 2921629
    Number of unmapped reads: 17765
    Number of unmapping reads with quality 0: 126628
    Average mapping quality for all the reads: 55.914
    Average mapping quality if reads with 0 map quality are filtered out: 58.447


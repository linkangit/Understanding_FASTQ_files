# Understanding, Reading, and Inferring FASTQ Files: A Detailed Tutorial

## Introduction

FASTQ files are the standard format for storing nucleotide sequence data along with quality scores from high-throughput sequencing platforms. Understanding these files is essential for anyone working in bioinformatics, genomics, or computational biology.

## What is a FASTQ File?

FASTQ is a text-based format that stores both biological sequences (usually DNA or RNA) and their corresponding quality scores. It's an extension of the FASTA format that includes quality information for each base call.

## FASTQ File Structure

Each sequence entry in a FASTQ file consists of exactly **four lines**:

### Line 1: Sequence Identifier
- Begins with `@` character
- Contains sequence identifier and optional description
- Format varies by sequencing platform

**Example:**
```
@SEQ_ID optional_description
```

Common identifier formats:
- **Illumina**: `@INSTRUMENT:RUN_ID:FLOWCELL_ID:LANE:TILE:X:Y READ:FILTERED:CONTROL:INDEX`
- **Example**: `@SRR123456.1 HWI-ST123:100:C0ABCACXX:1:1101:1234:2000 1:N:0:ATCACG`

### Line 2: Sequence
- The actual nucleotide sequence
- Contains A, C, G, T, and sometimes N (unknown base)

**Example:**
```
GATCCGAATCGATCGATCGATCGATCG
```

### Line 3: Separator
- Begins with `+` character
- Optionally repeats the sequence identifier (but usually just `+`)
- Acts as a separator between sequence and quality

**Example:**
```
+
```

### Line 4: Quality Scores
- ASCII-encoded quality scores for each base
- Same length as the sequence line
- Each character represents the quality of the corresponding base

**Example:**
```
IIIIIHHHGGGFFFFDDDCCCBBBAA
```

## Complete FASTQ Entry Example

```
@SRR123456.1 HWI-ST123:100:C0ABCACXX:1:1101:1234:2000 1:N:0:ATCACG
GATCCGAATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIHHHHHHGGGGGFFFFFFF@@@DDDCCCCBBBBAAA@@@###########
```

## Understanding Quality Scores

Quality scores represent the confidence in each base call. They're encoded using the **Phred quality score** system.

### Phred Quality Score

The relationship between quality score (Q) and error probability (P):

```
Q = -10 × log₁₀(P)
```

**Common Quality Score Values:**

| Phred Score | Error Probability | Accuracy | ASCII Character (Phred+33) |
|-------------|-------------------|----------|----------------------------|
| 10 | 1 in 10 | 90% | + |
| 20 | 1 in 100 | 99% | 5 |
| 30 | 1 in 1,000 | 99.9% | ? |
| 40 | 1 in 10,000 | 99.99% | I |

### ASCII Encoding

Quality scores are encoded as ASCII characters to save space. There are two main encoding schemes:

**Phred+33 (Sanger/Illumina 1.8+)** - Most common today:
- ASCII characters: `!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ`
- Quality range: 0-41
- Formula: `Quality = ASCII_value - 33`

**Phred+64 (Illumina 1.3-1.7)** - Older format:
- ASCII characters: `@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_'abcdefgh`
- Quality range: 0-40
- Formula: `Quality = ASCII_value - 64`

### Decoding Quality Scores

To convert an ASCII character to a quality score:

```python
ascii_char = 'I'
quality_score = ord(ascii_char) - 33  # For Phred+33
# quality_score = 40
```

## Reading FASTQ Files Programmatically

### Python Example (Using BioPython)

```python
from Bio import SeqIO

# Reading a FASTQ file
for record in SeqIO.parse("example.fastq", "fastq"):
    print(f"ID: {record.id}")
    print(f"Sequence: {record.seq}")
    print(f"Length: {len(record.seq)}")
    print(f"Quality scores: {record.letter_annotations['phred_quality']}")
    print(f"Average quality: {sum(record.letter_annotations['phred_quality']) / len(record.seq):.2f}")
    print("---")
```

### Python Example (Manual Parsing)

```python
def parse_fastq(filename):
    """Generator to parse FASTQ files"""
    with open(filename, 'r') as f:
        while True:
            # Read four lines at a time
            header = f.readline().strip()
            if not header:  # End of file
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            
            # Extract ID (remove @ symbol)
            seq_id = header[1:].split()[0]
            
            # Convert quality scores
            quality_scores = [ord(q) - 33 for q in quality]
            
            yield {
                'id': seq_id,
                'sequence': sequence,
                'quality': quality,
                'quality_scores': quality_scores
            }

# Usage
for read in parse_fastq("example.fastq"):
    print(f"ID: {read['id']}")
    print(f"Avg Quality: {sum(read['quality_scores']) / len(read['quality_scores']):.2f}")
```

### Python Example (Quality Filtering)

```python
def filter_by_quality(input_fastq, output_fastq, min_quality=20):
    """Filter reads by average quality score"""
    from Bio import SeqIO
    
    good_reads = []
    total_reads = 0
    
    for record in SeqIO.parse(input_fastq, "fastq"):
        total_reads += 1
        avg_quality = sum(record.letter_annotations['phred_quality']) / len(record.seq)
        
        if avg_quality >= min_quality:
            good_reads.append(record)
    
    # Write filtered reads
    SeqIO.write(good_reads, output_fastq, "fastq")
    
    print(f"Total reads: {total_reads}")
    print(f"Passed filter: {len(good_reads)} ({len(good_reads)/total_reads*100:.2f}%)")

# Usage
filter_by_quality("raw.fastq", "filtered.fastq", min_quality=25)
```

## Common Inferences and Quality Control

### 1. Basic Statistics

```python
from Bio import SeqIO
import numpy as np

def fastq_stats(filename):
    """Calculate basic statistics from FASTQ file"""
    lengths = []
    qualities = []
    gc_contents = []
    
    for record in SeqIO.parse(filename, "fastq"):
        # Length
        lengths.append(len(record.seq))
        
        # Average quality
        avg_q = np.mean(record.letter_annotations['phred_quality'])
        qualities.append(avg_q)
        
        # GC content
        gc_count = record.seq.count('G') + record.seq.count('C')
        gc_contents.append(gc_count / len(record.seq) * 100)
    
    print(f"Total reads: {len(lengths)}")
    print(f"Average read length: {np.mean(lengths):.2f} bp")
    print(f"Average quality: {np.mean(qualities):.2f}")
    print(f"Average GC content: {np.mean(gc_contents):.2f}%")
    
    return lengths, qualities, gc_contents
```

### 2. Per-Base Quality Distribution

```python
def per_base_quality(filename, max_reads=10000):
    """Calculate per-base quality across reads"""
    from collections import defaultdict
    
    position_qualities = defaultdict(list)
    
    for i, record in enumerate(SeqIO.parse(filename, "fastq")):
        if i >= max_reads:
            break
        
        for pos, qual in enumerate(record.letter_annotations['phred_quality']):
            position_qualities[pos].append(qual)
    
    # Calculate average quality per position
    avg_qualities = {}
    for pos, quals in position_qualities.items():
        avg_qualities[pos] = np.mean(quals)
    
    return avg_qualities
```

### 3. Detecting Quality Encoding

```python
def detect_quality_encoding(filename, num_reads=1000):
    """Detect if FASTQ uses Phred+33 or Phred+64"""
    min_char = 255
    max_char = 0
    
    with open(filename, 'r') as f:
        for i in range(num_reads * 4):
            line = f.readline()
            if i % 4 == 3:  # Quality line
                for char in line.strip():
                    min_char = min(min_char, ord(char))
                    max_char = max(max_char, ord(char))
    
    if min_char < 59:
        return "Phred+33 (Sanger/Illumina 1.8+)"
    else:
        return "Phred+64 (Illumina 1.3-1.7)"
```

## Common Quality Issues to Check

### 1. Low Quality Bases
- Bases with Phred score < 20 are considered low quality
- Common at the end of reads due to signal degradation

### 2. Adapter Contamination
- Sequencing adapters not properly removed
- Appears as over-represented sequences

### 3. Biased Base Composition
- First few bases often show bias due to priming
- GC content should be relatively uniform

### 4. Duplicate Sequences
- PCR duplicates from library preparation
- High duplication may indicate low library complexity

## Paired-End FASTQ Files

Paired-end sequencing produces two FASTQ files:
- `sample_R1.fastq` - Forward reads
- `sample_R2.fastq` - Reverse reads

The reads are paired by position (line number), not by ID necessarily.

```python
from Bio import SeqIO

# Read paired files simultaneously
forward_reads = SeqIO.parse("sample_R1.fastq", "fastq")
reverse_reads = SeqIO.parse("sample_R2.fastq", "fastq")

for fwd, rev in zip(forward_reads, reverse_reads):
    print(f"Forward: {fwd.id}")
    print(f"Reverse: {rev.id}")
    # Process paired reads together
```

## Useful Command-Line Tools

### 1. View First Few Reads
```bash
head -n 20 example.fastq  # Show first 5 reads (4 lines each)
```

### 2. Count Reads
```bash
wc -l example.fastq  # Divide by 4 for number of reads
```

### 3. FastQC (Quality Control)
```bash
fastqc example.fastq  # Generates comprehensive quality report
```

### 4. Seqtk (FASTQ manipulation)
```bash
# Subsample reads
seqtk sample example.fastq 10000 > subset.fastq

# Convert quality encoding
seqtk seq -V example.fastq > converted.fastq
```

## Best Practices

1. **Always validate file format** - Check that every sequence has exactly 4 lines
2. **Monitor quality scores** - Average quality > 30 is excellent, > 20 is acceptable
3. **Check read length distribution** - Unexpected variation may indicate problems
4. **Watch for adapter sequences** - These should be trimmed before analysis
5. **Compress files** - FASTQ files are large; use gzip compression (`file.fastq.gz`)
6. **Keep raw data** - Always preserve original FASTQ files before processing

## Conclusion

FASTQ files are the foundation of modern sequencing analysis. Understanding their structure, quality encoding, and how to parse them programmatically enables you to perform quality control, filtering, and preprocessing essential for downstream genomic analyses. Whether using specialized tools like BioPython or writing custom parsers, mastering FASTQ files is a crucial skill in bioinformatics.

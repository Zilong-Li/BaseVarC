BaseVarC - SNPs Calling From Low-Pass (<1.0x) WGS Data
===============================================================
**__Current Version: 1.0.0__**

BaseVarC was implemented in C++, aiming at speeding up variants calling from large-scale population, and was used in the [CMDB](https://db.cngb.org/cmdb/) project for calling variants from one million samples


## Installation

```
git clone --recursive https://github.com/Zilong-Li/BaseVarC.git
cd BaseVarC
./configure
make 
```

If everything goes well, you can find BaseVarC program in the `src` directory.

## Command-Line

```
BaseVarC
Contact: Zilong Li [zimusen94@gmail.com]
Usage  : BaseVarC <command> [options]

Commands:
         basetype       Variants Caller
         popmatrix      Create population matrix
         concat         Concat popmatrix
```

### Variants Calling

```
Commands: BaseVarC basetype
Usage   : BaseVarC basetype [options]

Options :
  --input,      -i         BAM/CRAM file list, one file per row
  --output,     -o         Output file prefix
  --reference,  -r         Reference file
  --region,     -s         Samtools-like region <chr:start-end>
  --group,      -g         Population group information <SampleID Group>
  --mapq,       -q <INT>   Mapping quality >= INT [10]
  --thread,     -t <INT>   Number of threads
  --batch,      -b <INT>   Number of samples each batch
  --maf,        -a <FLOAT> Minimum allele count frequency [min(0.001, 100/N, maf)]
  --load,                  Load data only
  --rerun,                 Read previous loaded data and rerun
  --keep_tmp,              Don't remove tmp files when basetype finished
  --verbose,    -v         Set verbose output
```


## Testing

In the tests directory, there is a script which contains a example using test data.

```
cd test/
sh test.sh
```

## Note on performance

RAM, run time and I/O all rest squarely on three parameters: `--region`, `--thread` and `--batch`. Depending on your situation, you can customize these parameters for exploiting your HPC servers.

- `--batch` : BaseVarC converts reads from BAM files into an internal temp format. This parameter control how many samples will be bundled as a batch. RAM is linear with this. Larger number means more RAM but less file pointers(I/O).
- `--region`: The longer the genomic region is given, the more RAM is used. Be aware that reading BAM files repeatedly is overhead. So you should split the chromosome into long region as possible as you can.
- `--thread`: The number of threads to use. RAM and I/O are linear with threads. The more threads are given, the faster BaseVarC is.

## License

BaseVarC and the code in this repo is available under a GPL3 license. For more information please see the [LICENSE](LICENSE)

## Citation

TBD

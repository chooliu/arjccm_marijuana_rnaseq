
Analysis code accompanying "Chronic marijuana use is associated with gene expression in bronchoalveolar lavage" (research letter; under review).

```
File Structure
├─ Data
│   └─ count_table.txt             # RNA-seq counts (gene x sample)
│   └─ participant_metadata.txt    # (sample x variable)
├─ Scripts                         # 22 scripts; data processing & analysis
├─ Results 
│   └─ ResultsTable1.xlsx          # pairwise diff. expression results 
│   └─ ResultsTable2.xlsx          # GO enrichment results
└─ Supplementary_Info.pdf          # detailed methods, Table1, QC, RT-PCR, etc.
```

Participant level metadata and raw RNA-seq data is also publicly available at GEO, accession [GSE155213](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155213).

# Cell Annotation

##### Description

The `cell_annotation_operator` annotates clusters based on a reference table of
marker vs. cell population correpondence.

##### Usage

Input projection|.
---|---
`y-axis`        | numeric, enrichment score (from the marker_enrichment_operator)
`row`           | factor, variable (channel, marker) 
`column`        | factor, group (cluster)

Input parameters|.
---|---
`annotation_level`        | Annotation level

Output relations|.
---|---
`population`        | Cell population
`prob`        | Posterior probability of being part of the population, per column (cluster)

##### See Also

[marker_enrichment_operator](https://github.com/tercen/marker_enrichment_operator)



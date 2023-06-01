## MPEG-7 shape dataset

This folder contains a shape dataset consisting of polygonal chains extracted from binary images of everyday objects comprising a famous benchmark dataset for computer vision known as MPEG-7 (n = 1400). Also included is a curated subset of MPEG-7 that we call the Animals dataset (n = 480), which keeps only the animal classes of MPEG-7.

Polygonal chains are provided in .Rdata format. We use a subset of MPEG-7 converted to the MATLAB .mat format.

* `MPEG7closed.mat`: MPEG-7 polygonal chains in MATLAB variable format, comprising a total of 65 shape classes; each shape class has 20 observations.
* `MPEG7closedLabels.txt`: A detailed description of the data.
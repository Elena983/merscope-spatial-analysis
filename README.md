# merscope-spatial-analysis

Data [from](https://console.cloud.google.com/storage/browser/vz-ffpe-showcase/HumanLiverCancerPatient1?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false)

Cell classification
Analysis of Image-based Spatial Data in Seurat

![iScreen Shoter - Google Chrome - 231109105350](https://github.com/Elena983/merscope-spatial-analysis/assets/68946912/641ae14c-2e67-4c6f-b26a-b0e216a34f9c)

This code was create with Voyager package

@Manual{,
  title = {Voyager: From geospatial to spatial omics},
  author = {Lambda Moses and Kayla Jackson and Lior Pachter},
  year = {2022},
  note = {R package version 1.0.3},
  url = {https://github.com/pachterlab/voyager},
}

Moses L, Jackson K, Pachter L (2022). 
Voyager: From geospatial to spatial omics. 
R package version 1.0.3, 
https://github.com/pachterlab/voyager

analis_sfe

![image](https://user-images.githubusercontent.com/68946912/214247285-c7cff14a-cfd6-4d49-a40f-909d121e0609.png)

Percentage of blank features per cell
![image](https://user-images.githubusercontent.com/68946912/214247503-0cd9a653-f9e7-4b5b-9c6f-d1a613f9b03e.png)

Plot that zoomed in patch to visually inspect cell-cell contiguity
![image](https://user-images.githubusercontent.com/68946912/214247727-190ea134-bef2-4a6d-ac7b-032e1d8d83ee.png)

nCounts and nGenes have sizable positive Moran’s I’s
![image](https://user-images.githubusercontent.com/68946912/214247830-443c217c-b903-4928-9d12-91522e491b45.png)

Plot the top genes with positive spatial autocorrelation
![image](https://user-images.githubusercontent.com/68946912/214247932-50081a1d-9c23-47de-8399-8388550597ad.png)

Plot top gene loadings in each PC 
Many of these genes seem to be related to the endothelium??
![image](https://user-images.githubusercontent.com/68946912/214248198-ba1c8d60-8eac-409e-a2f8-b0ca19e3778e.png)

Plot the first 4 PCs in space
![image](https://user-images.githubusercontent.com/68946912/214248307-b6cb780c-5713-4456-8de6-964273596aa1.png)

Makeshift function to plot a zoomed in patch spatial structure in PC2 and PC3 at a smaller scale
![image](https://user-images.githubusercontent.com/68946912/214248442-535a86c0-e72a-4017-83c0-28f49489f86f.png)

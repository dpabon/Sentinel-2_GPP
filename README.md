# On the potential of Sentinel-2 for estimating Gross Primary Production #

This is the code repository for the scientific manuscript: "On the potential of Sentinel-2 for estimating Gross Primary Production". For further details see: XXXX

1. Example scripts for the post-processing of Sentinel-2 images are presented in the [post-processing](post-processing) folder.

2. The scripts to generate the Sentinel-2 cubes (A netcdf file per site with 4 dimensions: Lat, Lon, Time, Variable) are presented in the [cubes_generation](cube_generation) folder.

3. The scripts to predict GPP using linear regressions and random forest are presented in the [GPP_prediction](GPP_prediction) folder.

4. In the folder [final_models](final_models) we present the random forest models pre-trained to predict GPP. You can use the function [sentinel-2_gpp_function.R](final_models/sentinel-2_gpp_function.R) to upscale GPP using any Sentinel-2 L2A product from [Copernicus](https://scihub.copernicus.eu/). An example of how to use the function is presented in the script [example_upscaling_GPP_Sentinel2.R](final_models/example_upscaling_GPP_Sentinel2.R) 
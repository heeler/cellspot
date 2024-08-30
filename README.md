To install this utility:

`pip install git+https://github.com/heeler/cellspot.git#egg=cellspot`

This  package is intended to find Nuclei using the blue channel of the images. 
It then creates a mask for each nuclei found and grows each mask to create an additional 
buffer around the nuclei. Each mask is then used with a threshold value to detect the presensce of 
inclusion bodies, with an intensity greater than the defined threshold. 


usage: 

`cell_spot [-h] [-r RED_THRESHOLD] [-c CLIP_TO] -o OUTPUT_FOLDER [-p PIXEL_CYCLES] filename`

    Required Args:
    filename - a tif input file with rgb channels
    -o output/directory - the output directory to write json results in 

    Optional Args:
    -r RED_THRESHOLD - the threshold to use for inclusion bodies to count in intensity
    -c CLIP_TO - a tuple of 2 flows, aka, (5.0, 90.0). 
    -p PIXEL_CYCLES - the number of cycles to grow the mask, aka, 6 

To install this utility:

`pip install git+https://github.com/heeler/cellspot.git#egg=cellspot`

This  package is intended to find Nuclei using the blue channel of the images. 
It then creates a mask for each nuclei found and grows each mask to create an additional 
buffer around the nuclei. Each mask is then used with a threshold value to detect the presensce of 
inclusion bodies, with an intensity greater than the defined threshold. 


usage: 

`cell_spot [-h] [-r RED_THRESHOLD] [-c CLIP_TO] -o OUTPUT_FOLDER [-p PIXEL_CYCLES] [-d DEAD_CELLS] filename`

    Required Args:
    filename - a tif input file with rgb channels
    -o output/directory - the output directory to write json results in 

    Optional Args:
    -r RED_THRESHOLD - the threshold to use for inclusion bodies to count in intensity
    -c CLIP_TO - a tuple of 2 flows, aka, (5.0, 90.0). 
    -p PIXEL_CYCLES - the number of cycles to grow the mask, aka, 6 
    -d DEAD_CELLS - a minimum number of pixels for the nuclei to occupy for the cell to be considered alive, dead cells are removed.


Recommendations: 
1) work from the terminal 
2) Use miniconda to create a python 3.11 environment (`conda create -n myenv python=3.11`)
3) Activate your new environment (`conda activate myenv`)
4) install the package and dependencies into the environment (`pip install git+https://github.com/heeler/cellspot.git#egg=cellspot`)

Possibly helpful bash commands:

assumptions 
1) the images are in /Users/myname/images 
2) I want the outputs in /Users/myname/results

```bash
for file in /Users/myname/images/*.tif
do
cell_spot $file -o /Users/myname/results
done
```

also if you wanted to run the command with different red thresholds you could do

```bash
for file in /Users/myname/images/*.tif
do
cell_spot $file -o /Users/myname/results_red_50 -r 50 
done
```

or similarly with clipping values for the blue channel

```bash
for file in /Users/myname/images/*.tif
do
cell_spot $file -o /Users/myname/results_clip_10_to_80 -c (10, 80)
done
```

Also see [cellpose](https://cellpose.org)
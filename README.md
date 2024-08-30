To install this utility:

`pip install git+https://github.com/heeler/cellspot.git#egg=cellspot`

This  package is intended to find Nuclei using the blue channel of the images. 
It then creates a mask for each nuclei found and grows each mask to create an additional 
buffer around the nuclei. Each mask is then used with a threshold value to detect the presensce of 
inclusion bodies, with an intensity greater than the defined threshold. 

script usage:


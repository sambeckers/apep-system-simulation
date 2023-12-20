
# Apep System Simulation
## Name of the Member:

- Vincent van Rie
- Divyansh Srivastava
- Sam Beckers
  
## Minimum Requirements

- **Morphological Reproduction**: The primary objective is to accurately simulate the morphology of the Apep system. The code must be capable of reproducing the observed shape (serpent) and maybe the structure of the system as viewed from various angles. The code will use gravity, stellar evolution, hydrodynamics and radiative transfer codes from the AMUSE framework.
  
- **Evolutionary Animation**: The code should generate a time-lapse animation detailing the visual representation of the system's evolution from its inception to its current state. 

## Additional Objectives

- **Supernova Simulation**: After meeting the minimum requirements, the next goal is to simulate the supernova of the WC8(carbon-burning) Wolf-Rayet star within the system.
  
- **Morphological Impact Analysis**: Post-supernova, the resultant changes in the system’s morphology.

- **X-Ray Flux**: Post-processing element derived from the evolving system.

## Make the movie
*Command Line Code that put the plots together.*
*This might need ffmpeg installation*
ffmpeg -r 2 -f image2 -s 1920x1080 -i hydro_outflow_step_%d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p output_movie.mp4
